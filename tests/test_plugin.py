# Copyright (C) 2026 Yves Sorge
# SPDX-License-Identifier: AGPL-3.0-only

"""Installation-level tests for the SeaSenseLib reader entry point."""

from __future__ import annotations

from importlib.metadata import metadata
from pathlib import Path

import numpy as np
import seasenselib as ssl


FORMAT_KEY = "sbe-cnv-seabirdscientific"
PACKAGE_NAME = "seasenselib-seabirdscientific"


def test_installed_distribution_declares_agpl() -> None:
    """Built package metadata exposes the intended strong-copyleft license."""
    assert metadata(PACKAGE_NAME)["License-Expression"] == "AGPL-3.0-only"


def _write_fixture(path: Path) -> Path:
    """Write a minimal, deterministic CNV file."""
    path.write_text(
        "\n".join(
            (
                "* Sea-Bird SBE 37 Data File:",
                "# nvalues = 2",
                "# name 0 = timeS: Time [seconds]",
                "# name 1 = t090C: Temperature [deg C]",
                "# interval = seconds: 1",
                "# start_time = Jan 01 2024 00:00:00",
                "# bad_flag = -9.990e-29",
                "*END*",
                "0 10.25",
                "1 10.5",
                "",
            )
        ),
        encoding="utf-8",
    )
    return path


def test_reader_is_discovered_as_plugin() -> None:
    """The installed entry point appears in SeaSenseLib's public registry."""
    info = next(item for item in ssl.list_readers() if item["key"] == FORMAT_KEY)

    assert info["is_plugin"] is True
    assert info["extension"] is None

    from seasenselib.readers import get_reader_by_format_key

    reader_class = get_reader_by_format_key(FORMAT_KEY)
    assert reader_class is not None
    assert reader_class.__module__ == "seasenselib_seabirdscientific.reader"


def test_public_api_reads_cnv_with_plugin(tmp_path: Path) -> None:
    """SeaSenseLib can instantiate and run the plugin through its public API."""
    source = _write_fixture(tmp_path / "sample.cnv")

    dataset = ssl.read(
        source,
        file_format=FORMAT_KEY,
        perform_default_postprocessing=False,
    )

    np.testing.assert_array_equal(dataset["t090C"].values, [10.25, 10.5])
    np.testing.assert_array_equal(
        dataset.time.values,
        np.asarray(
            ["2024-01-01T00:00:00", "2024-01-01T00:00:01"],
            dtype="datetime64[ns]",
        ),
    )
