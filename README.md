# gwassociation: GW–EM Association Framework

`gwassociation` is a Python package for evaluating associations between
**gravitational-wave (GW)** events and **electromagnetic (EM)** transients using
a Bayesian framework.

Given a GW skymap and one or more EM candidates, `gwassociation` computes:

- **Spatial overlap** between the EM position and the GW localization
- **Distance overlap** between GW distance posteriors and EM redshift
- **Temporal overlap** between GW time and EM detection time
- **Posterior odds** and **association probability** for each candidate

---

## Quick Start

### Install

From PyPI:

```console
pip install gwassociation
```

From a local checkout:

```console
git clone https://github.com/<your-username>/gwassociation.git
cd gwassociation
pip install -e .
```

### Minimal Example (Python API)

```python
from gwassociation import Association

# GW skymap FITS file (2D or 3D)
gw_file = "S250818k_bayestar.fits.gz"

# EM transient information
transient_info = {
    "name": "AT2024abc",
    "ra": 192.42625,         # deg
    "dec": 34.82472,         # deg
    "z": 0.438,              # redshift
    "z_err": 0.005,          # redshift uncertainty
    "time": 1242442967.447,  # EM detection time (GPS)
    "gw_time": 1242442965.0  # GW trigger time (GPS)
}

assoc = Association(gw_file, transient_info)
results = assoc.compute_odds(
    em_model="kilonova",         # "kilonova", "grb", or "afterglow"
    prior_odds=1.0,              # prior odds for association
    chance_coincidence_rate=1e-4 # expected chance coincidence rate
)

print(f"P(Associated)   = {results['confidence']:.1%}")
print(f"Posterior odds  = {results['posterior_odds']:.3e}")
print(f"Bayes factor    = {results['bayes_factor']:.3e}")
print(f"Log10(odds)     = {results['log_posterior_odds']:.2f}")
print(f"Decision        = {'ASSOCIATED' if results['associated'] else 'NOT ASSOCIATED'}")
```

### Rank Multiple Candidates

```python
candidates = [
    {
        "name": "AT2024abc",
        "ra": 192.42625,
        "dec": 34.82472,
        "z": 0.438,
        "z_err": 0.005,
        "time": 1242442967.447,
        "gw_time": 1242442965.0,
    },
    {
        "name": "AT2024def",
        "ra": 192.43000,
        "dec": 34.82000,
        "z": 0.440,
        "z_err": 0.006,
        "time": 1242442970.0,
        "gw_time": 1242442965.0,
    },
]

assoc = Association(gw_file, {"gw_time": 1242442965.0})
rankings = assoc.rank_candidates(candidates)

for i, r in enumerate(rankings, 1):
    cand = r["candidate"]
    print(
        f"{i}. {cand.name}: "
        f"P(assoc) = {r['probability']:.1%}, "
        f"log10(odds) = {r['log_odds']:.2f}"
    )
```

---

## Command-Line Interface

A simple CLI is provided via `console_scripts`:

```console
gwassociation \
  --gw-file S250818k_bayestar.fits.gz \
  --ra 192.42625 \
  --dec 34.82472 \
  --z 0.438 \
  --z-err 0.005 \
  --time 1242442967.447 \
  --gw-time 1242442965.0 \
  --model kilonova \
  --out results/ \
  --verbose
```

This will:

- Load the GW skymap
- Compute spatial, distance, and temporal overlaps
- Print association statistics
- Save results and plots in `results/`

---

## What the Package Does

Given a GW event and EM candidate(s), `gwassociation` computes:

- **Spatial overlap** \(I_\Omega\): how well the EM position lies within the GW localization.
- **Distance overlap** \(I_{DL}\): agreement between GW distance posterior and the EM distance inferred from redshift.
- **Temporal overlap** \(I_t\): likelihood of the observed time delay between GW and EM signals under a chosen EM model.

These are combined into a Bayes factor and posterior odds:

```text
Posterior Odds = Prior Odds × (I_Ω × I_DL × I_t) / P_chance
P(associated | data) = Posterior Odds / (1 + Posterior Odds)
```

---

## API Overview

The core high-level API is the `Association` class.

```python
from gwassociation import Association

assoc = Association(
    gw_file="skymap.fits",
    transient_info={"ra": 120.5, "dec": -30.2, "time": 1234567890.0},
)
```

Key methods:

- `compute_odds(**kwargs) -> dict`: compute posterior odds and return result statistics.
- `rank_candidates(candidates_list) -> list`: rank multiple candidate transients by association probability.
- `plot_skymap(out_file: str = "skymap.png")`: plot the GW skymap and highlight the transient position.

## Plotting

The `gwassociation.plots` module provides helper functions:

- `plot_association_summary(results_dict, output_file="association_summary.png")`
- `plot_candidate_ranking(candidates_results, output_file="candidate_ranking.png")`

## Dependencies

Core dependencies are declared in `setup.cfg` and include `numpy`, `scipy`,
`matplotlib`, `astropy`, `healpy`, `ligo.skymap`, and `h5py`.

## Development

Clone and install in editable mode:

```console
git clone https://github.com/<your-username>/gwassociation.git
cd gwassociation
pip install -e .
```

Run tests:

```console
pytest
```

Run the working example:

```console
python working_example.py
```

## Citation

If you use `gwassociation` in your work, please cite:

```tex
@software{gwassociation,
  author  = {Ignacio Ma{\~n}ana, Kaitlyn Pak},
  title   = {gwassociation: GW–EM Association Framework},
  year    = {2025},
  url     = {https://github.com/<your-username>/gwassociation},
  version = {0.1.0}
}
```

And cite the relevant methodology papers on GW–EM associations, GW localization,
and EM counterpart modeling.

## License

This project is licensed under the MIT License. See `LICENSE` for details.
