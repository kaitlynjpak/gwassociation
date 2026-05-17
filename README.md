# gwassociation: GW–EM Association Framework

`gwassociation` is a lightweight Python package for evaluating candidate
associations between gravitational-wave (GW) sky maps and electromagnetic (EM)
transients. It computes sky-position, distance, and timing overlap terms and
combines them into posterior odds.

## Install

```console
pip install -e .
```

Optional integrations are available as extras:

```console
pip install -e ".[healpy,ligo,hdf5,dev]"
```

## Minimal Python API

```python
from gwassociation import Association

assoc = Association(
    "fits_files/S190425z_bayestar.fits.gz,0",
    {
        "name": "AT2024abc",
        "ra": 120.5,
        "dec": -30.2,
        "z": 0.05,
        "time": 1234567890.0,
        "gw_time": 1234567880.0,
    },
)
results = assoc.compute_odds(em_model="kilonova")
print(results["confidence"], results["posterior_odds"])
```

## Command line

```console
gwassociation \
  --gw-file fits_files/S190425z_bayestar.fits.gz,0 \
  --ra 120.5 \
  --dec -30.2 \
  --z 0.05 \
  --time 1234567890 \
  --gw-time 1234567880 \
  --out results
```

## Examples

The maintained example is intentionally small and uses the public API:

```console
python examples/minimal_script.py fits_files/S190425z_bayestar.fits.gz,0
```

## Development

```console
pip install -e ".[dev]"
pytest -q tests
```

The repository intentionally excludes generated plots, ad-hoc analysis outputs,
and notebooks so the package stays compact. The compact FITS fixture under `fits_files/` is retained because it is required for regression coverage.

## Core modules

- `gwassociation.association.Association`: high-level API.
- `gwassociation.io`: sky-map and transient containers/loaders.
- `gwassociation.analysis`: spatial, line-of-sight distance, temporal, and odds calculations.
- `gwassociation.stats`: lower-level overlap and prior utilities.
- `gwassociation.plots` / `gwassociation.plotting`: optional plotting helpers.

## License

See the project license file if one is added by the maintainers.
