# SeaSenseLib CNV reader based on seabirdscientific

[![CI](https://github.com/ysorge/seasenselib-seabirdscientific/actions/workflows/ci.yml/badge.svg)](https://github.com/ysorge/seasenselib-seabirdscientific/actions/workflows/ci.yml)
[![Python 3.10–3.13](https://img.shields.io/badge/python-3.10%E2%80%933.13-blue.svg)](https://www.python.org/downloads/)
[![License: AGPL-3.0](https://img.shields.io/badge/license-AGPL--3.0-blue.svg)](LICENSE)
[![SeaSenseLib](https://img.shields.io/badge/SeaSenseLib-GitHub-181717?logo=github)](https://github.com/ocean-uhh/seasenselib)
[![SeaSenseLib on PyPI](https://img.shields.io/pypi/v/seasenselib?label=SeaSenseLib%20on%20PyPI)](https://pypi.org/project/seasenselib/)
[![SeaSenseLib DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20044197.svg)](https://doi.org/10.5281/zenodo.20044197)

Minimal beta plugin providing the reader `sbe-cnv-seabirdscientific`
for [SeaSenseLib](https://github.com/ocean-uhh/seasenselib).
Numeric CNV decoding is delegated to the `seabirdscientific`
package; SeaSenseLib integration adds strict source validation,
deterministic time handling, metadata and provenance.

## Installation

The currently supported Python versions are 3.10 through 3.13. See
[Known limitations](#known-limitations) for the constraints that determine
this range.

```bash
python -m pip install .
```

For development:

```bash
python -m pip install -e '.[test]'
pytest
```

The package registers itself through the `seasenselib.readers` entry-point
group.

Verify installation with:

```bash
seasenselib list readers
```


## Usage

```python
import seasenselib as ssl

dataset = ssl.read(
    "cast.cnv",
    file_format="sbe-cnv-seabirdscientific",
)
```

## Caveats

The reader remains an explicit beta and does not claim automatic
`.cnv` detection.

The development version of SeaSenseLib may already bundle the same beta
reader. When both copies are present, plugin discovery intentionally selects
the installed plugin class and may report that it overrides the built-in
class.

## Known limitations

- Full fixture regression currently covers `seabirdscientific` 2.7.8. The
  documented 3.x xarray return type is supported by an adapter but remains
  outside the dependency range until it has passed the complete fixture
  corpus.
- A real `timeN` fixture is still needed to independently verify its epoch.
- Julian-day channels require a reference year from `start_time`; automatic
  rollover across New Year is not inferred.
- The strict validation pass makes the reader slower than the
  pycnv-backed reader on large files.
- Python 3.9 is not supported because SeaSenseLib 0.6.0 requires
  `mhkit>=1.0.0` (for Sea-Bird HEX reader), whose available releases require Python 3.10 or newer.
- The plugin currently requires exactly `seabirdscientific==2.7.8`.
  Versions 2.7.9 and 2.7.10 constrain SciPy to `~=1.13.1`, which conflicts
  with the `scipy>=1.14.0` requirement of MHKiT 1.x. SciPy 1.13.1 also has no
  Python 3.13 wheel. The pin must remain until a newer compatible manufacturer
  release has passed the complete regression corpus.
- Python 3.13 is covered by this plugin's CI and its dependency graph resolves
  to binary wheels, but `seabirdscientific` 2.7.8 does not explicitly list
  Python 3.13 in its PyPI classifiers.

## Related project and citation

This reader is an independently maintained plugin for
[SeaSenseLib](https://github.com/ocean-uhh/seasenselib). SeaSenseLib is also
available [on PyPI](https://pypi.org/project/seasenselib/) and has a
[citable Zenodo record](https://doi.org/10.5281/zenodo.20044197).

The DOI `10.5281/zenodo.20044197` identifies **SeaSenseLib**, not this plugin.
This repository does not currently have its own DOI. Citation metadata for
the plugin and the referenced SeaSenseLib software is provided in
[`CITATION.cff`](CITATION.cff).
