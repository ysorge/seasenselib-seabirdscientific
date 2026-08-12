# SeaSenseLib CNV reader based on seabirdscientific

Minimal beta plugin providing the reader `sbe-cnv-seabirdscientific`
for [SeaSenseLib](https://github.com/ocean-uhh/seasenselib).
Numeric CNV decoding is delegated to the `seabirdscientific`
package; SeaSenseLib integration adds strict source validation,
deterministic time handling, metadata and provenance.

## Installation

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

- Full fixture regression currently covers `seabirdscientific` 2.7.x. The
  documented 3.x xarray return type is supported by an adapter but has not yet
  passed the complete fixture corpus.
- A real `timeN` fixture is still needed to independently verify its epoch.
- Julian-day channels require a reference year from `start_time`; automatic
  rollover across New Year is not inferred.
- The strict validation pass makes the reader slower than the
  pycnv-backed reader on large files.
