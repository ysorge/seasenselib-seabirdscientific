# Copyright (C) 2026 Yves Sorge
# SPDX-License-Identifier: AGPL-3.0-only
# Portions derived from SeaSenseLib; see NOTICE.

"""Beta CNV reader backed by :mod:`seabirdscientific`.

``seabirdscientific`` deliberately provides a small CNV value parser.  This
module keeps that parser as the decoding backend and adds the format handling
SeaSenseLib needs around it: validation, deterministic time coordinates,
coordinates, bad-value handling, and lossless raw-header provenance.
"""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime
import hashlib
from importlib import metadata as importlib_metadata
import logging
from pathlib import Path
import re
import tempfile
from typing import Any, Iterator

import numpy as np
import pandas as pd
import xarray as xr

import seasenselib.parameters as params
from seasenselib.readers.base import AbstractReader
from seasenselib.readers.utils.conductivity_units import infer_conductivity_unit


logger = logging.getLogger(__name__)

_NAME_RE = re.compile(
    r"^#\s*name\s+(?P<index>\d+)\s*=\s*(?P<name>[^:]+):\s*(?P<details>.*)$"
)
_NVALUE_RE = re.compile(r"^#\s*nvalues\s*=\s*(?P<count>\d+)")
_INTERVAL_RE = re.compile(
    r"^#\s*interval\s*=\s*(?P<kind>[^:]+):\s*(?P<value>[-+0-9.eE]+)"
)
_START_TIME_RE = re.compile(
    r"^#\s*start_time\s*=\s*"
    r"(?P<value>[A-Za-z]{3}\s+\d{1,2}\s+\d{4}\s+\d{2}:\d{2}:\d{2})"
)
_BAD_FLAG_RE = re.compile(r"^#\s*bad_flag\s*=\s*(?P<value>\S+)")
_NUMBER_RE = re.compile(
    r"[-+]?(?:(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?|nan|inf)",
    re.IGNORECASE,
)
_HEADER_DATETIME_FORMAT = "%b %d %Y %H:%M:%S"
_TIME_CONSISTENCY_TOLERANCE_SECONDS = 5.0


@dataclass(frozen=True)
class _CnvColumn:
    """One ``# name`` declaration in source-column order."""

    index: int
    name: str
    key: str
    original_label: str
    description: str
    unit: str


@dataclass(frozen=True)
class _CnvHeader:
    """Header evidence required to build a reproducible dataset."""

    raw: str
    columns: tuple[_CnvColumn, ...]
    declared_sample_count: int | None
    actual_sample_count: int
    interval_kind: str | None
    interval_value: float | None
    start_time: datetime | None
    nmea_time: datetime | None
    system_time: datetime | None
    upload_time: datetime | None
    bad_flag: float | None
    latitude: float | None
    longitude: float | None
    encoding: str
    encoding_basis: str
    source_sha256: str
    sanitizations: tuple[dict[str, Any], ...]


@dataclass(frozen=True)
class _TimeCoordinate:
    """Calculated time values plus an explicit provenance record."""

    values: np.ndarray
    source_name: str
    source_type: str
    reference_source: str | None = None
    reference_time: datetime | None = None
    validation_source: str | None = None
    validation_max_difference_seconds: float | None = None


class _CnvSourceParser:
    """Inspect and normalize the CNV text format before numeric decoding.

    This component owns byte decoding, header parsing, row validation and the
    small, recorded rewrites required by the upstream parser. It never changes
    the source file; normalized input is exposed through a temporary file.
    """

    @staticmethod
    def _parse_datetime(value: str) -> datetime | None:
        """Parse a Sea-Bird header timestamp, returning ``None`` if invalid."""
        try:
            return datetime.strptime(value.strip(), _HEADER_DATETIME_FORMAT)
        except ValueError:
            return None

    @staticmethod
    def _find_header_datetime(
        text: str,
        labels: tuple[str, ...],
    ) -> datetime | None:
        """Find the first timestamp associated with one of ``labels``."""
        label_pattern = "|".join(re.escape(label) for label in labels)
        pattern = re.compile(
            rf"^[*#]+\s*(?:{label_pattern})\s*=\s*"
            r"(?P<value>[A-Za-z]{3}\s+\d{1,2}\s+\d{4}\s+\d{2}:\d{2}:\d{2})",
            re.MULTILINE | re.IGNORECASE,
        )
        match = pattern.search(text)
        return (
            _CnvSourceParser._parse_datetime(match.group("value"))
            if match
            else None
        )

    @staticmethod
    def _parse_position(text: str, coordinate: str) -> float | None:
        """Parse a header latitude or longitude into signed decimal degrees."""
        pattern = re.compile(
            rf"^[*#]+\s*(?:NMEA\s+)?{coordinate}\s*(?:=|:)\s*"
            r"(?P<degrees>\d{1,3})\s+(?P<minutes>\d+(?:\.\d+)?)\s*"
            r"(?P<hemisphere>[NSEW])",
            re.MULTILINE | re.IGNORECASE,
        )
        match = pattern.search(text)
        if not match:
            return None
        value = float(match.group("degrees")) + float(match.group("minutes")) / 60.0
        if match.group("hemisphere").upper() in {"S", "W"}:
            value = -value
        return value

    @staticmethod
    def _decode_source(
        raw: bytes,
        requested_encoding: str,
    ) -> tuple[str, str, str]:
        """Decode CNV bytes deterministically and explain an automatic fallback."""
        if requested_encoding != "auto":
            return raw.decode(requested_encoding), requested_encoding, "reader argument"

        try:
            return raw.decode("utf-8"), "utf-8", "valid UTF-8"
        except UnicodeDecodeError:
            cp1252_text = raw.decode("cp1252")
            # Sea-Bird's legacy DOS exports use CP437 for the theta glyph.  Select
            # it only when the same declaration explicitly describes sigma-theta;
            # otherwise prefer the normal Windows code page.
            theta_evidence = re.search(
                r"^#\s*name\s+\d+\s*=\s*sigma-[^\x00-\x7f]+.*"
                r"Density\s*\[sigma-theta",
                cp1252_text,
                re.MULTILINE | re.IGNORECASE,
            )
            if theta_evidence:
                return (
                    raw.decode("cp437"),
                    "cp437",
                    "legacy sigma-theta declaration",
                )
            return cp1252_text, "cp1252", "invalid UTF-8; Windows text fallback"

    @staticmethod
    def _parse_column(line: str, used_keys: set[str]) -> _CnvColumn | None:
        """Parse one ``# name`` declaration and create a collision-free key."""
        match = _NAME_RE.match(line)
        if not match:
            return None

        name = match.group("name").strip()
        details = match.group("details").strip()
        left_bracket = details.find("[")
        if left_bracket >= 0:
            description = details[:left_bracket].rstrip(" ,")
            right_bracket = details.find("]", left_bracket + 1)
            unit_end = right_bracket if right_bracket >= 0 else len(details)
            unit = details[left_bracket + 1 : unit_end].strip()
        else:
            description = details
            unit = ""

        key = name
        suffix = 1
        while key in used_keys:
            key = f"{name}_{suffix}"
            suffix += 1
        used_keys.add(key)
        return _CnvColumn(
            index=int(match.group("index")),
            name=name,
            key=key,
            original_label=details,
            description=description,
            unit=unit,
        )

    @staticmethod
    def _validate_data_line(
        line: str,
        expected_columns: int,
        line_number: int,
    ) -> None:
        """Require numeric content and the exact header-declared row width."""
        matches = list(_NUMBER_RE.finditer(line))
        remainder = _NUMBER_RE.sub("", line).strip()
        if remainder:
            raise ValueError(
                f"CNV data line {line_number} contains non-numeric content: "
                f"{remainder!r}"
            )
        if len(matches) != expected_columns:
            raise ValueError(
                f"CNV data line {line_number} has {len(matches)} values; "
                f"the header declares {expected_columns} columns"
            )

    @staticmethod
    def _replace_nvalues(
        lines: list[str],
        actual_count: int,
    ) -> tuple[list[str], bool]:
        """Return header lines with ``nvalues`` set to the validated count."""
        replaced = False
        result: list[str] = []
        for line in lines:
            if _NVALUE_RE.match(line):
                newline = "\n" if line.endswith(("\n", "\r")) else ""
                result.append(f"# nvalues = {actual_count}{newline}")
                replaced = True
            else:
                result.append(line)

        if not replaced:
            insertion = next(
                (index for index, line in enumerate(result) if _NAME_RE.match(line)),
                next(
                    (
                        index
                        for index, line in enumerate(result)
                        if line.strip().startswith("*END*")
                    ),
                    len(result),
                ),
            )
            result.insert(insertion, f"# nvalues = {actual_count}\n")
        return result, replaced

    @staticmethod
    def inspect(
        path: Path,
        *,
        requested_encoding: str,
        sanitize_input: bool,
    ) -> tuple[_CnvHeader, str | None]:
        """Inspect, validate, and optionally produce normalized upstream input."""
        source_bytes = path.read_bytes()
        text, encoding, encoding_basis = _CnvSourceParser._decode_source(
            source_bytes,
            requested_encoding,
        )
        lines = text.splitlines(keepends=True)
        end_index = next(
            (
                index
                for index, line in enumerate(lines)
                if line.strip().startswith("*END*")
            ),
            None,
        )
        if end_index is None:
            raise ValueError("CNV header does not contain an *END* marker")

        header_lines = lines[: end_index + 1]
        raw_header = "\n".join(line.rstrip("\r\n") for line in header_lines)
        used_keys: set[str] = set()
        columns = tuple(
            column
            for line in header_lines
            if (
                column := _CnvSourceParser._parse_column(
                    line.rstrip("\r\n"),
                    used_keys,
                )
            )
            is not None
        )
        columns = tuple(sorted(columns, key=lambda column: column.index))
        if not columns:
            raise ValueError("CNV header does not declare any '# name' columns")
        expected_indices = list(range(len(columns)))
        actual_indices = [column.index for column in columns]
        if actual_indices != expected_indices:
            raise ValueError(
                "CNV column indices must be contiguous and start at zero; "
                f"found {actual_indices}"
            )

        data_lines: list[tuple[int, str]] = []
        blank_data_lines = 0
        for line_number, line in enumerate(lines[end_index + 1 :], end_index + 2):
            if not line.strip():
                blank_data_lines += 1
                continue
            if line.startswith(("*", "#")):
                continue
            _CnvSourceParser._validate_data_line(
                line,
                len(columns),
                line_number,
            )
            data_lines.append((line_number, line))
        if not data_lines:
            raise ValueError("CNV file contains no numeric data rows after *END*")

        declared_count = None
        for line in header_lines:
            match = _NVALUE_RE.match(line)
            if match:
                declared_count = int(match.group("count"))
                break

        interval_kind = None
        interval_value = None
        start_time = None
        bad_flag = None
        for line in header_lines:
            stripped = line.rstrip("\r\n")
            if interval_match := _INTERVAL_RE.match(stripped):
                interval_kind = interval_match.group("kind").strip().lower()
                interval_value = float(interval_match.group("value"))
            if start_match := _START_TIME_RE.match(stripped):
                start_time = _CnvSourceParser._parse_datetime(
                    start_match.group("value")
                )
            if flag_match := _BAD_FLAG_RE.match(stripped):
                try:
                    bad_flag = float(flag_match.group("value"))
                except ValueError:
                    logger.warning(
                        "Ignoring non-numeric CNV bad_flag: %s",
                        flag_match.group("value"),
                    )

        actions: list[dict[str, Any]] = []
        needs_rewrite = encoding.lower().replace("_", "-") not in {"utf-8", "utf8"}
        if needs_rewrite:
            actions.append(
                {
                    "action": "transcode_to_utf8_for_upstream_parser",
                    "source_encoding": encoding,
                    "basis": encoding_basis,
                }
            )

        actual_count = len(data_lines)
        if declared_count != actual_count:
            if not sanitize_input:
                raise ValueError(
                    f"CNV declares nvalues={declared_count}, but contains "
                    f"{actual_count} data rows; enable sanitize_input to parse "
                    "actual rows"
                )
            needs_rewrite = True
            actions.append(
                {
                    "action": "replace_nvalues",
                    "declared": declared_count,
                    "actual": actual_count,
                }
            )

        sanitized_header = list(header_lines)
        for index, line in enumerate(sanitized_header):
            if _START_TIME_RE.match(line.rstrip("\r\n")):
                content = line.rstrip("\r\n")
                line_ending = line[len(content) :]
                trimmed = content.rstrip(" \t")
                if trimmed != content:
                    if sanitize_input:
                        sanitized_header[index] = trimmed + line_ending
                        needs_rewrite = True
                        actions.append(
                            {
                                "action": "trim_start_time_trailing_whitespace",
                                "line": index + 1,
                            }
                        )

        if blank_data_lines and sanitize_input:
            needs_rewrite = True
            actions.append(
                {
                    "action": "remove_blank_data_lines",
                    "count": blank_data_lines,
                }
            )

        if (
            declared_count != actual_count or declared_count is None
        ) and sanitize_input:
            sanitized_header, _ = _CnvSourceParser._replace_nvalues(
                sanitized_header,
                actual_count,
            )

        normalized_text = None
        if needs_rewrite and sanitize_input:
            normalized_text = "".join(sanitized_header)
            if normalized_text and not normalized_text.endswith(("\n", "\r")):
                normalized_text += "\n"
            normalized_text += "".join(line for _, line in data_lines)

        header = _CnvHeader(
            raw=raw_header,
            columns=columns,
            declared_sample_count=declared_count,
            actual_sample_count=actual_count,
            interval_kind=interval_kind,
            interval_value=interval_value,
            start_time=start_time,
            nmea_time=_CnvSourceParser._find_header_datetime(
                raw_header,
                ("NMEA UTC (Time)", "NMEA UTC"),
            ),
            system_time=_CnvSourceParser._find_header_datetime(
                raw_header,
                ("System UTC",),
            ),
            upload_time=_CnvSourceParser._find_header_datetime(
                raw_header,
                ("System UpLoad Time",),
            ),
            bad_flag=bad_flag,
            latitude=_CnvSourceParser._parse_position(raw_header, "Latitude"),
            longitude=_CnvSourceParser._parse_position(raw_header, "Longitude"),
            encoding=encoding,
            encoding_basis=encoding_basis,
            source_sha256=hashlib.sha256(source_bytes).hexdigest(),
            sanitizations=tuple(actions),
        )
        return header, normalized_text

    @staticmethod
    @contextmanager
    def upstream_input(
        path: Path,
        normalized_text: str | None,
    ) -> Iterator[Path]:
        """Yield original or temporary normalized input and always clean up."""
        if normalized_text is None:
            yield path
            return

        temporary_path: Path | None = None
        try:
            with tempfile.NamedTemporaryFile(
                mode="w",
                encoding="utf-8",
                suffix=".cnv",
                prefix="seasenselib-cnv-",
                delete=False,
            ) as handle:
                handle.write(normalized_text)
                temporary_path = Path(handle.name)
            yield temporary_path
        finally:
            if temporary_path is not None:
                temporary_path.unlink(missing_ok=True)


class _SeabirdScientificAdapter:
    """Isolate the supported ``seabirdscientific`` CNV API variants."""

    @staticmethod
    def decode(
        path: Path,
        columns: tuple[_CnvColumn, ...],
    ) -> tuple[list[np.ndarray], str, str]:
        """Return source-column arrays across the 2.x and 3.x APIs.

        The adapter deliberately returns arrays in the source header's column
        order. This prevents changes in an upstream container type from leaking
        into SeaSenseLib's dataset construction and provenance handling.
        """
        try:
            from seabirdscientific import instrument_data
        except ImportError as exc:
            raise ImportError(
                "SbeCnvSeabirdScientificReader requires seabirdscientific>=2.7.0"
            ) from exc

        try:
            version = importlib_metadata.version("seabirdscientific")
        except importlib_metadata.PackageNotFoundError:
            version = "unknown"

        try:
            if hasattr(instrument_data, "read_cnv_file"):
                parsed = instrument_data.read_cnv_file(path)
                backend = "read_cnv_file"
            else:
                parsed = instrument_data.cnv_to_instrument_data(path)
                backend = "cnv_to_instrument_data"
        except Exception as exc:
            raise ValueError(
                f"seabirdscientific {version} failed to parse CNV file: {exc}"
            ) from exc

        if isinstance(parsed, xr.Dataset):
            missing = [column.key for column in columns if column.key not in parsed]
            if missing:
                raise ValueError(
                    "seabirdscientific did not return declared CNV columns: "
                    f"{missing}"
                )
            arrays = [np.asarray(parsed[column.key].values) for column in columns]
        else:
            measurements = getattr(parsed, "measurements", None)
            if not isinstance(measurements, dict):
                raise TypeError(
                    "Unsupported seabirdscientific CNV result; expected "
                    "xarray.Dataset or InstrumentData.measurements"
                )
            missing = [
                column.key for column in columns if column.key not in measurements
            ]
            if missing:
                raise ValueError(
                    "seabirdscientific did not return declared CNV columns: "
                    f"{missing}"
                )
            arrays = [
                np.asarray(measurements[column.key].values) for column in columns
            ]
        return arrays, backend, version


class _CnvTimeCodec:
    """Pure helpers for locating, masking, and decoding CNV time channels."""

    @staticmethod
    def seconds_to_datetimes(values: np.ndarray, epoch: datetime) -> np.ndarray:
        """Interpret numeric values as elapsed seconds from ``epoch``."""
        numeric = np.asarray(values, dtype=float)
        result = pd.Timestamp(epoch) + pd.to_timedelta(numeric, unit="s")
        return np.asarray(result, dtype="datetime64[ns]")

    @staticmethod
    def julian_days_to_datetimes(
        values: np.ndarray,
        reference_year: int,
    ) -> np.ndarray:
        """Interpret one-based Julian-day values within ``reference_year``."""
        numeric = np.asarray(values, dtype=float)
        result = pd.Timestamp(
            year=reference_year,
            month=1,
            day=1,
        ) + pd.to_timedelta(numeric - 1.0, unit="D")
        return np.asarray(result, dtype="datetime64[ns]")

    @staticmethod
    def max_difference_seconds(left: np.ndarray, right: np.ndarray) -> float:
        """Return the maximum absolute difference of valid timestamp pairs."""
        valid = ~(np.isnat(left) | np.isnat(right))
        if not np.any(valid):
            return float("inf")
        differences = np.abs(
            (left[valid] - right[valid])
            .astype("timedelta64[ns]")
            .astype(np.int64)
        )
        return float(np.max(differences) / 1_000_000_000)

    @classmethod
    def progression_difference_seconds(
        cls,
        left: np.ndarray,
        right: np.ndarray,
    ) -> float:
        """Compare timestamp progressions while ignoring their absolute offset."""
        valid = ~(np.isnat(left) | np.isnat(right))
        if np.count_nonzero(valid) < 2:
            return float("inf")
        left_valid = left[valid]
        right_valid = right[valid]
        return cls.max_difference_seconds(
            left_valid - left_valid[0] + np.datetime64("2000-01-01", "ns"),
            right_valid - right_valid[0] + np.datetime64("2000-01-01", "ns"),
        )

    @staticmethod
    def find_array(
        arrays: dict[str, np.ndarray],
        names: tuple[str, ...],
    ) -> tuple[str, np.ndarray] | None:
        """Find the first named array using case-insensitive matching."""
        lowered = {name.lower(): name for name in arrays}
        for candidate in names:
            actual = lowered.get(candidate.lower())
            if actual is not None:
                return actual, arrays[actual]
        return None

    @staticmethod
    def mask_bad_value(
        values: np.ndarray,
        bad_flag: float | None,
    ) -> np.ndarray:
        """Replace the header-declared missing-value sentinel with NaN."""
        numeric = np.asarray(values, dtype=float)
        if bad_flag is None:
            return numeric
        result = numeric.copy()
        result[result == bad_flag] = np.nan
        return result


class _CnvHeaderMetadataParser:
    """Extract structured sensor facts while retaining the complete header."""

    @staticmethod
    def parse_sensors(raw_header: str) -> dict[str, dict[str, Any]]:
        """Return selected sensor XML fields keyed by channel number."""
        starts = list(
            re.finditer(
                r"#\s*<sensor\s+Channel=[\"'](?P<channel>\d+)[\"'][^>]*>",
                raw_header,
                re.IGNORECASE,
            )
        )
        result: dict[str, dict[str, Any]] = {}
        for index, match in enumerate(starts):
            end = (
                starts[index + 1].start()
                if index + 1 < len(starts)
                else len(raw_header)
            )
            block = raw_header[match.end() : end]
            metadata: dict[str, Any] = {"channel": int(match.group("channel"))}
            sensor = re.search(
                r"#\s*<(?P<type>\w+Sensor)\s+SensorID=[\"']"
                r"(?P<id>[^\"']*)[\"']",
                block,
                re.IGNORECASE,
            )
            if sensor:
                metadata["sensor_type"] = sensor.group("type")
                metadata["sensor_id"] = sensor.group("id")
            for source, target in (
                ("SerialNumber", "serial_number"),
                ("CalibrationDate", "calibration_date"),
            ):
                value = re.search(
                    rf"#\s*<{source}>(?P<value>[^<]+)</{source}>",
                    block,
                    re.IGNORECASE,
                )
                if value:
                    metadata[target] = value.group("value").strip()
            comment = re.search(r"#\s*<!--\s*(?P<value>.*?)\s*-->", block)
            if comment:
                metadata["sensor_comment"] = comment.group("value").strip()
            result[match.group("channel")] = metadata
        return result


class SbeCnvSeabirdScientificReader(AbstractReader):
    """Read a Sea-Bird CNV file with ``seabirdscientific`` (beta).

    The beta reader is selected explicitly with format key
    ``sbe-cnv-seabirdscientific``.  Automatic ``.cnv`` detection intentionally
    remains assigned to :class:`SbeCnvReader` while the new path is evaluated.
    """

    def __init__(
        self,
        input_file: str,
        sanitize_input: bool = True,
        encoding: str = "auto",
        time_source: str = "auto",
        time_q_epoch: int = 2000,
        time_n_epoch: int = 1970,
        use_default_latitude: bool | None = None,
        default_latitude: float | None = None,
        mapping: dict | None = None,
        **kwargs,
    ):
        """Configure source normalization, time interpretation, and mapping.

        ``time_q_epoch`` and ``time_n_epoch`` remain explicit because the CNV
        labels alone do not establish one consistent epoch across observed
        files and upstream documentation. No epoch is inferred from values.
        """
        default_latitude_configured = (
            use_default_latitude is not None
            or default_latitude is not None
            or "fix_missing_coords" in kwargs
        )
        if "fix_missing_coords" in kwargs:
            use_default_latitude = bool(kwargs.pop("fix_missing_coords"))
        if time_q_epoch not in {1970, 2000} or time_n_epoch not in {1970, 2000}:
            raise ValueError("time_q_epoch and time_n_epoch must be 1970 or 2000")

        super().__init__(input_file, mapping, **kwargs)
        self._sanitize_input = bool(sanitize_input)
        self._encoding = encoding
        self._time_source = time_source
        self._time_q_epoch = time_q_epoch
        self._time_n_epoch = time_n_epoch
        self._default_latitude_configured = default_latitude_configured
        self._use_default_latitude = (
            bool(use_default_latitude)
            if use_default_latitude is not None
            else default_latitude is not None
        )
        self._default_latitude = 45.0 if default_latitude is None else default_latitude
        self._raw_header = None
        self._raw_metadata_blocks: dict[str, Any] = {}
        self._raw_metadata_variables: dict[str, Any] = {}
        self._validate_file()

    @classmethod
    def _get_valid_extensions(cls) -> tuple[str, ...]:
        """Return the source extensions accepted during file validation."""
        return (".cnv",)

    @classmethod
    def reader_args(cls) -> list[dict[str, Any]]:
        """Describe reader-specific arguments for the public reader API."""
        return [
            cls._reader_arg(
                "sanitize_input",
                "bool",
                True,
                "Normalize known technical CNV parser incompatibilities.",
            ),
            cls._reader_arg(
                "encoding",
                "str",
                "auto",
                "Source encoding; auto uses UTF-8 and evidence-based legacy fallbacks.",
            ),
            cls._reader_arg(
                "time_source",
                "str",
                "auto",
                "Time channel to use, or auto for evidence-based selection.",
                choices=(
                    "auto",
                    "timeS",
                    "timeJ",
                    "timeJV2",
                    "timeSCP",
                    "timeQ",
                    "timeK",
                    "timeN",
                    "interval",
                ),
            ),
            cls._reader_arg(
                "time_q_epoch",
                "int",
                2000,
                "Epoch year for timeQ; configurable because upstream metadata "
                "conflicts.",
                choices=("1970", "2000"),
            ),
            cls._reader_arg(
                "time_n_epoch",
                "int",
                1970,
                "Epoch year for timeN; configurable because upstream metadata "
                "conflicts.",
                choices=("1970", "2000"),
            ),
        ]

    def _absolute_time_candidates(
        self,
        arrays: dict[str, np.ndarray],
        header: _CnvHeader,
    ) -> list[tuple[str, str, np.ndarray]]:
        """Decode supported absolute channels without choosing among them."""
        candidates: list[tuple[str, str, np.ndarray]] = []
        if header.start_time is not None:
            for julian_name in ("timeJ", "timeJV2", "timeSCP"):
                julian = _CnvTimeCodec.find_array(arrays, (julian_name,))
                if julian is not None:
                    name, values = julian
                    candidates.append(
                        (
                            name,
                            "julian_days",
                            _CnvTimeCodec.julian_days_to_datetimes(
                                _CnvTimeCodec.mask_bad_value(values, header.bad_flag),
                                header.start_time.year,
                            ),
                        )
                    )
        for names, source_type, epoch_year in (
            (("timeQ",), "seconds_since_configured_timeQ_epoch", self._time_q_epoch),
            (("timeK",), "seconds_since_2000", 2000),
            (("timeN",), "seconds_since_configured_timeN_epoch", self._time_n_epoch),
        ):
            found = _CnvTimeCodec.find_array(arrays, names)
            if found is not None:
                name, values = found
                candidates.append(
                    (
                        name,
                        source_type,
                        _CnvTimeCodec.seconds_to_datetimes(
                            _CnvTimeCodec.mask_bad_value(values, header.bad_flag),
                            datetime(epoch_year, 1, 1),
                        ),
                    )
                )
        return candidates

    def _elapsed_time_coordinate(
        self,
        source_name: str,
        elapsed_values: np.ndarray,
        absolute_candidates: list[tuple[str, str, np.ndarray]],
        header: _CnvHeader,
    ) -> _TimeCoordinate | None:
        """Anchor ``timeS`` to header or corroborating absolute evidence."""
        elapsed_values = _CnvTimeCodec.mask_bad_value(
            elapsed_values,
            header.bad_flag,
        )
        references = (
            ("header:nmea_utc", header.nmea_time),
            ("header:system_utc", header.system_time),
            ("header:start_time", header.start_time),
        )
        best: tuple[float, str, datetime, np.ndarray, str] | None = None
        for reference_name, reference_time in references:
            if reference_time is None:
                continue
            candidate = _CnvTimeCodec.seconds_to_datetimes(
                elapsed_values,
                reference_time,
            )
            if absolute_candidates:
                for absolute_name, _, absolute in absolute_candidates:
                    difference = _CnvTimeCodec.max_difference_seconds(
                        candidate,
                        absolute,
                    )
                    if best is None or difference < best[0]:
                        best = (
                            difference,
                            reference_name,
                            reference_time,
                            candidate,
                            absolute_name,
                        )
            elif reference_name == "header:start_time":
                return _TimeCoordinate(
                    candidate,
                    source_name,
                    "seconds_since_start_time",
                    reference_name,
                    reference_time,
                )

        if best is not None and best[0] <= _TIME_CONSISTENCY_TOLERANCE_SECONDS:
            return _TimeCoordinate(
                best[3],
                source_name,
                "seconds_since_header_reference",
                best[1],
                best[2],
                best[4],
                best[0],
            )

        # An absolute channel may be rounded to whole seconds while elapsed
        # time retains fractions.  Fuse only when both progressions agree.
        for absolute_name, _, absolute in absolute_candidates:
            elapsed_from_zero = _CnvTimeCodec.seconds_to_datetimes(
                elapsed_values - elapsed_values[0],
                datetime(2000, 1, 1),
            )
            progression_difference = _CnvTimeCodec.progression_difference_seconds(
                elapsed_from_zero,
                absolute,
            )
            if progression_difference <= _TIME_CONSISTENCY_TOLERANCE_SECONDS:
                valid = ~np.isnat(absolute)
                if not np.any(valid):
                    continue
                first = int(np.flatnonzero(valid)[0])
                fused = absolute[first] + pd.to_timedelta(
                    elapsed_values - elapsed_values[first],
                    unit="s",
                )
                return _TimeCoordinate(
                    np.asarray(fused, dtype="datetime64[ns]"),
                    source_name,
                    "elapsed_time_anchored_to_absolute_channel",
                    absolute_name,
                    None,
                    absolute_name,
                    progression_difference,
                )

        if header.start_time is not None:
            if absolute_candidates:
                logger.warning(
                    "CNV time sources disagree by more than %.3f seconds; "
                    "using timeS with the declared start_time",
                    _TIME_CONSISTENCY_TOLERANCE_SECONDS,
                )
            return _TimeCoordinate(
                _CnvTimeCodec.seconds_to_datetimes(
                    elapsed_values,
                    header.start_time,
                ),
                source_name,
                "seconds_since_start_time",
                "header:start_time",
                header.start_time,
                best[4] if best is not None else None,
                best[0] if best is not None else None,
            )
        return None

    def _calculate_time(
        self,
        arrays: dict[str, np.ndarray],
        header: _CnvHeader,
    ) -> _TimeCoordinate:
        """Select a deterministic time coordinate and record its evidence."""
        absolute = self._absolute_time_candidates(arrays, header)
        requested = self._time_source.lower()

        if requested == "auto":
            elapsed = _CnvTimeCodec.find_array(arrays, ("timeS",))
            if elapsed is not None:
                result = self._elapsed_time_coordinate(
                    elapsed[0],
                    elapsed[1],
                    absolute,
                    header,
                )
                if result is not None:
                    return result
            if absolute:
                name, source_type, values = absolute[0]
                return _TimeCoordinate(values, name, source_type)
        elif requested == "times":
            elapsed = _CnvTimeCodec.find_array(arrays, ("timeS",))
            if elapsed is None:
                raise ValueError("Requested CNV time source timeS is not present")
            result = self._elapsed_time_coordinate(
                elapsed[0],
                elapsed[1],
                absolute,
                header,
            )
            if result is None:
                raise ValueError(
                    "timeS is present but has no deterministic reference time"
                )
            return result
        elif requested == "interval":
            pass
        else:
            selected = next(
                (item for item in absolute if item[0].lower() == requested),
                None,
            )
            if selected is None:
                raise ValueError(
                    f"Requested CNV time source {self._time_source!r} is unavailable"
                )
            return _TimeCoordinate(selected[2], selected[0], selected[1])

        if (
            header.interval_kind == "seconds"
            and header.interval_value is not None
            and header.start_time is not None
        ):
            values = (
                np.arange(header.actual_sample_count, dtype=float)
                * header.interval_value
            )
            return _TimeCoordinate(
                _CnvTimeCodec.seconds_to_datetimes(values, header.start_time),
                "start_time + interval",
                "start_time_plus_seconds_interval",
                "header:start_time",
                header.start_time,
            )
        raise ValueError(
            "CNV file has no deterministic supported time reference. "
            "Provide timeS/timeJ/timeQ/timeK/timeN or a start_time plus "
            "seconds interval."
        )

    @staticmethod
    def _coordinate_array(
        arrays: dict[str, np.ndarray],
        aliases: tuple[str, ...],
        fallback: float | None,
        bad_flag: float | None,
    ) -> tuple[np.ndarray | float, str | None]:
        """Return a coordinate channel or a scalar header fallback."""
        found = _CnvTimeCodec.find_array(arrays, aliases)
        if found is not None:
            return _CnvTimeCodec.mask_bad_value(found[1], bad_flag), found[0]
        return (np.nan if fallback is None else fallback), None

    def _assign_time_metadata(self, dataset: xr.Dataset, time: _TimeCoordinate) -> None:
        """Attach the time decision and validation evidence to the dataset."""
        dataset.attrs["cnv_time_source_variable"] = time.source_name
        dataset.attrs["cnv_time_source_type"] = time.source_type
        dataset[params.TIME].attrs["source_variable"] = time.source_name
        dataset[params.TIME].attrs["source_type"] = time.source_type
        if time.reference_source:
            dataset.attrs["cnv_time_reference_source"] = time.reference_source
            dataset[params.TIME].attrs["reference_source"] = time.reference_source
        if time.reference_time:
            value = time.reference_time.isoformat()
            dataset.attrs["cnv_time_reference"] = value
            dataset[params.TIME].attrs["reference_time"] = value
        if time.validation_source:
            dataset.attrs["cnv_time_validation_source"] = time.validation_source
        if time.validation_max_difference_seconds is not None:
            dataset.attrs["cnv_time_validation_max_difference_seconds"] = (
                time.validation_max_difference_seconds
            )

    def _assign_global_metadata(
        self,
        dataset: xr.Dataset,
        header: _CnvHeader,
        backend: str,
        version: str,
    ) -> None:
        """Attach selected header facts and assemble raw metadata provenance."""
        model = re.search(
            r"^\*\s*Sea-Bird\s+SBE\s*(?P<value>.*?)\s+Data File:",
            header.raw,
            re.MULTILINE | re.IGNORECASE,
        )
        software = re.search(
            r"^\*\s*Software Version\s+(?P<value>.+?)\s*$",
            header.raw,
            re.MULTILINE | re.IGNORECASE,
        )
        if model:
            dataset.attrs["cnv_sbe_model"] = f"SBE {model.group('value').strip()}"
        if software:
            dataset.attrs["cnv_software_version"] = software.group("value").strip()
        for name, value in (
            ("cnv_start_date", header.start_time),
            ("cnv_upload_date", header.upload_time),
            ("cnv_nmea_date", header.nmea_time),
        ):
            if value is not None:
                dataset.attrs[name] = value.strftime("%Y-%m-%d %H:%M:%S")
        if header.interval_kind == "seconds" and header.interval_value is not None:
            dataset.attrs["cnv_interval_seconds"] = header.interval_value
        if header.bad_flag is not None:
            dataset.attrs["cnv_bad_flag"] = header.bad_flag

        sensors = _CnvHeaderMetadataParser.parse_sensors(header.raw)
        for channel, metadata in sensors.items():
            dataset.attrs[f"cnv_sensor_{channel}"] = metadata

        self._raw_metadata_blocks = {
            "attributes": {
                "parser": "seabirdscientific",
                "parser_api": backend,
                "parser_version": version,
                "source_encoding": header.encoding,
                "encoding_detection_basis": header.encoding_basis,
                "source_file_sha256": header.source_sha256,
                "declared_sample_count": header.declared_sample_count,
                "actual_sample_count": header.actual_sample_count,
                "interval_kind": header.interval_kind,
                "interval_value": header.interval_value,
                "bad_flag": header.bad_flag,
            },
            "calibration": {"sensors": sensors} if sensors else None,
            "other": {"sanitization": list(header.sanitizations)},
        }

    @staticmethod
    def _verify_conductivity_units(dataset: xr.Dataset) -> None:
        """Warn if values conflict with a declared conductivity unit.

        This check never changes data or units. It only surfaces a possible
        source inconsistency for review by the normal SeaSenseLib pipeline.
        """
        conductivity_names = {
            "c0mS/cm",
            "c0S/m",
            "c1mS/cm",
            "c1S/m",
            "cond0S/m",
            "cond1S/m",
            "cond0mS/cm",
            "cond1mS/cm",
        }
        for name in dataset.data_vars:
            original = dataset[name].attrs.get("cnv_original_name", name)
            declared = dataset[name].attrs.get("units", "")
            if original not in conductivity_names or not declared:
                continue
            try:
                infer_conductivity_unit(dataset[name].values, declared=declared)
            except ValueError as exc:
                logger.warning(
                    "Could not verify conductivity units for '%s': %s",
                    name,
                    exc,
                )

    def _load_data(self) -> xr.Dataset:
        """Decode the CNV source and construct the unmapped xarray dataset."""
        source_path = Path(self.input_file)
        header, normalized_text = _CnvSourceParser.inspect(
            source_path,
            requested_encoding=self._encoding,
            sanitize_input=self._sanitize_input,
        )
        self._raw_header = header.raw

        with _CnvSourceParser.upstream_input(
            source_path,
            normalized_text,
        ) as upstream_path:
            decoded_arrays, backend, version = _SeabirdScientificAdapter.decode(
                upstream_path,
                header.columns,
            )

        if len(decoded_arrays) != len(header.columns):
            raise ValueError(
                "seabirdscientific returned "
                f"{len(decoded_arrays)} channels for {len(header.columns)} CNV columns"
            )
        if any(len(values) != header.actual_sample_count for values in decoded_arrays):
            lengths = [len(values) for values in decoded_arrays]
            raise ValueError(
                "seabirdscientific returned channel lengths that do not match actual "
                f"CNV rows ({header.actual_sample_count}): {lengths}"
            )

        arrays = {
            column.key: np.asarray(values, dtype=float)
            for column, values in zip(header.columns, decoded_arrays)
        }
        time = self._calculate_time(arrays, header)
        latitude, latitude_key = self._coordinate_array(
            arrays,
            ("latitude", "lat"),
            header.latitude,
            header.bad_flag,
        )
        longitude, longitude_key = self._coordinate_array(
            arrays,
            ("longitude", "lon"),
            header.longitude,
            header.bad_flag,
        )
        latitude_is_array = np.ndim(latitude) > 0
        longitude_is_array = np.ndim(longitude) > 0
        dataset = self._get_xarray_dataset_template(
            time.values,
            None,
            np.nan if latitude_is_array else latitude,
            np.nan if longitude_is_array else longitude,
        )
        if latitude_is_array:
            dataset = dataset.assign_coords(latitude=(params.TIME, latitude))
        if longitude_is_array:
            dataset = dataset.assign_coords(longitude=(params.TIME, longitude))
        self._assign_time_metadata(dataset, time)

        self._raw_metadata_variables = {}
        coordinate_targets = {
            key: target
            for key, target in (
                (latitude_key, params.LATITUDE),
                (longitude_key, params.LONGITUDE),
            )
            if key is not None
        }
        coordinate_keys = set(coordinate_targets)
        for column in header.columns:
            values = arrays[column.key]
            if column.key not in coordinate_keys:
                # Keep the Sea-Bird scan flag as raw QC evidence.  For all
                # other variables, the declared missing-value sentinel is NaN.
                if column.name.lower() != "flag":
                    values = _CnvTimeCodec.mask_bad_value(values, header.bad_flag)
                dataset[column.key] = ([params.TIME], values)
            target_name = coordinate_targets.get(column.key, column.key)
            target = dataset[target_name]
            target.attrs["cnv_original_name"] = column.name
            target.attrs["cnv_original_label"] = column.original_label
            target.attrs["cnv_original_unit"] = column.unit
            if column.description:
                target.attrs["long_name"] = column.description
            if column.unit:
                target.attrs["units"] = column.unit
            if column.name.lower() == "flag" and header.bad_flag is not None:
                target.attrs["cnv_bad_flag"] = header.bad_flag
            self._raw_metadata_variables[column.key] = {
                "column_index": column.index,
                "original_name": column.name,
                "original_label": column.original_label,
                "description": column.description,
                "units": column.unit,
                "role": (
                    "coordinate" if column.key in coordinate_keys else "measurement"
                ),
            }

        self._assign_global_metadata(dataset, header, backend, version)
        self._verify_conductivity_units(dataset)
        return dataset

    @classmethod
    def format_key(cls) -> str:
        """Return the explicit format-selector key for this beta reader."""
        return "sbe-cnv-seabirdscientific"

    @classmethod
    def format_name(cls) -> str:
        """Return the human-readable reader name."""
        return "SeaBird CNV (seabirdscientific, beta)"

    @classmethod
    def file_extension(cls) -> None:
        """Disable automatic extension detection during the beta period."""
        # Explicit beta opt-in; avoid taking automatic .cnv detection from pycnv.
        return None

    @classmethod
    def format_mappings(cls) -> dict[str, list[str]]:
        """Map known CNV channel labels to SeaSenseLib parameters."""
        return {
            params.TEMPERATURE: ["t090C", "t068", "t190C", "t168", "tv290C"],
            params.SALINITY: ["sal00", "sal11"],
            params.CONDUCTIVITY: [
                "c0mS/cm",
                "c0S/m",
                "c1mS/cm",
                "c1S/m",
                "cond0S/m",
                "cond1S/m",
                "cond0mS/cm",
                "cond1mS/cm",
            ],
            params.PRESSURE: [
                "prdM",
                "prDM",
                "prSM",
                "prM",
                "pr50M",
                "pr200M",
                "pr350M",
            ],
            params.DEPTH: ["depSM"],
            params.OXYGEN: [
                "sbeox0V",
                "sbeox0",
                "sbeox0ML/L",
                "sbeox0Mm/Kg",
                "sbeox1V",
                "sbeox1ML/L",
                "sbeopoxML/L",
                "sbeopoxMm/Kg",
                "sbeopoxPS",
            ],
            params.TURBIDITY: ["turbWETntu0"],
            params.FLUORESCENCE: ["flECO-AFL"],
            params.DENSITY: [
                "sigma-t00",
                "sigma-t11",
                "sigma-theta00",
                "sigma-theta11",
                "sigma-Θ00",
            ],
            params.POTENTIAL_TEMPERATURE: ["potemp090C", "potemp190C"],
        }
