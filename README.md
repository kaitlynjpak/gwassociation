# GW-Assoc: GW–EM Association Framework

A Python package for evaluating associations between **gravitational-wave (GW)** events and **electromagnetic (EM)** transients using a Bayesian framework.

Given a GW skymap and one or more EM candidates, `gw-assoc` computes:
- **Spatial overlap** between the EM position and the GW localization
- **Distance overlap** between GW distance posteriors and EM redshift
- **Temporal overlap** between GW time and EM detection time
- **Posterior odds** and **association probability** for each candidate

---

## Quick Start

### Install

From PyPI:

pip install gw-assoc

git clone https://github.com/<your-username>/gw-assoc.git
cd gw-assoc/gwPackage
pip install -e .### Minimal Example (Python API)

from gw_assoc import Association

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

# Create association object
assoc = Association(gw_file, transient_info)

# Compute association odds
results = assoc.compute_odds(
    em_model="kilonova",         # "kilonova", "grb", or "afterglow"
    prior_odds=1.0,              # prior odds for association
    chance_coincidence_rate=1e-4 # expected chance coincidence rate
)

print(f"P(Associated)   = {results['confidence']:.1%}")
print(f"Posterior odds  = {results['posterior_odds']:.3e}")
print(f"Bayes factor    = {results['bayes_factor']:.3e}")
print(f"Log10(odds)     = {results['log_posterior_odds']:.2f}")
print(f"Decision        = {'ASSOCIATED' if results['associated'] else 'NOT ASSOCIATED'}")### Rank Multiple Candidates

candidates = [
    {"name": "AT2024abc", "ra": 192.42625, "dec": 34.82472, "z": 0.438, "z_err": 0.005,
     "time": 1242442967.447, "gw_time": 1242442965.0},
    {"name": "AT2024def", "ra": 192.43000, "dec": 34.82000, "z": 0.440, "z_err": 0.006,
     "time": 1242442970.0, "gw_time": 1242442965.0},
]

assoc = Association(gw_file, {"gw_time": 1242442965.0})
rankings = assoc.rank_candidates(candidates)

for i, r in enumerate(rankings, 1):
    cand = r["candidate"]
    print(
        f"{i}. {cand.name}: "
        f"P(assoc) = {r['probability']:.1%}, "
        f"log10(odds) = {r['log_odds']:.2f}"
    )---

## Command-Line Interface

A simple CLI is provided via `console_scripts`:

gw-assoc \
  --gw-file S250818k_bayestar.fits.gz \
  --ra 192.42625 \
  --dec 34.82472 \
  --z 0.438 \
  --z-err 0.005 \
  --time 1242442967.447 \
  --gw-time 1242442965.0 \
  --model kilonova \
  --out results/ \
  --verboseThis will:
- Load the GW skymap
- Compute spatial, distance, and temporal overlaps
- Print association statistics
- Save results and plots in `results/`

---

## What the Package Does

Given a GW event and EM candidate(s), `gw-assoc` computes:

- **Spatial overlap** \(I_\Omega\):  
  How well the EM position lies within the GW localization.

- **Distance overlap** \(I_{DL}\):  
  Agreement between GW distance posterior (from 3D skymaps) and the EM distance inferred from redshift.

- **Temporal overlap** \(I_t\):  
  Likelihood of the observed time delay between GW and EM signals under a chosen EM model:
  - `kilonova`: peaks ~1 day after merger
  - `grb`: prompt emission ~seconds after merger
  - `afterglow`: peaks around days and decays as a power law

These are combined into a Bayes factor and posterior odds:

\[
\text{Posterior Odds} = \text{Prior Odds} \times \frac{I_\Omega \, I_{DL} \, I_t}{P_{\text{chance}}}
\]

\[
P(\text{associated} \mid \text{data}) = \frac{\text{Posterior Odds}}{1 + \text{Posterior Odds}}
\]

where \(P_{\text{chance}}\) encodes the chance-coincidence probability.

---

## API Overview

The core high-level API is the `Association` class.

### `Association`

from gw_assoc import Association**Constructor**

Association(
    gw_file: str,
    transient_info: dict | None = None,
    secondary_skymap: str | None = None,
    secondary_event_time: float | None = None,
)- `gw_file`: path to primary GW skymap (FITS; 2D or 3D)
- `transient_info`: dict with EM transient data:
  - `ra`, `dec` (deg; required for point sources)
  - `z`, `z_err` (optional redshift info)
  - `time` (EM detection time, GPS)
  - `gw_time` (GW event time, GPS; optional)
- `secondary_skymap`: optional second skymap for skymap–skymap coincidence
- `secondary_event_time`: event time for the secondary skymap

**Methods**

- `compute_odds(**kwargs) -> dict`  
  Compute posterior odds and return a results dictionary.

  Key kwargs:
  - `em_model`: `"kilonova"`, `"grb"`, `"afterglow"`
  - `prior_odds`: prior association odds (float)
  - `chance_coincidence_rate`: chance-coincidence rate
  - `H0_uncertainty`: uncertainty on H0 (km/s/Mpc) used in distance overlap

- `rank_candidates(candidates_list) -> list`  
  Rank multiple candidate transients by association probability.

- `plot_skymap(out_file: str = "skymap.png")`  
  Plot the GW skymap and highlight the EM transient position.

---

## Plotting

The `gw_assoc.plots` module provides helper functions:

- `plot_association_summary(results_dict, output_file="association_summary.png")`  
  Creates a 2x2 panel summary:
  - Overlap integrals bar chart (I_Ω, I_DL, I_t)
  - Bayes factor on Jeffreys scale
  - Posterior odds and association probability
  - Summary table of key statistics

- `plot_candidate_ranking(candidates_results, output_file="candidate_ranking.png")`  
  Visualize multiple candidates’ odds and probabilities.

---

## Dependencies

Core dependencies (installed automatically via `pip`):

- `numpy`
- `scipy`
- `matplotlib`
- `astropy`
- `healpy`
- `ligo.skymap`
- `h5py`

Optional / recommended:

- `pandas` (for convenient table handling)
- `seaborn` (for prettier plots)

See `setup.cfg` / `pyproject.toml` for exact version pins.

---

## Development

Clone and install in editable mode:

git clone https://github.com/<your-username>/gw-assoc.git
cd gw-assoc/gwPackage
pip install -e ".[dev]"   # if you define dev extrasRun tests:

python test_gw_assoc.py
# or, if using pytest:
pytestRun the working example:

python working_example.pyThis will:
- Create test data (if needed)
- Run several association scenarios
- Generate demonstration plots

---

## Real-World Workflow (Typical Use)

1. **Get GW skymap** from GraceDB (`bayestar.fits.gz` or `bilby.fits.gz`).
2. **Prepare EM candidates** with RA/Dec, redshift, and times.
3. **Create `Association`** with the GW file and transient info.
4. **Compute odds and rank candidates.**
5. **Generate plots** (skymaps, summary, candidate rankings).
6. **Use results** to prioritize follow-up or interpret associations.

---

## Citation

If you use `gw-assoc` in your work, please cite:
tex
@software{gw_assoc,
  author  = {Ignacio Ma{\~n}ana, Kaitlyn Pak},
  title   = {GW-Assoc: GW–EM Association Framework},
  year    = {2025},
  url     = {https://github.com/<your-username>/gw-assoc},
  version = {0.1.0}
}And cite the relevant methodology papers on GW–EM associations, GW localization, and EM counterpart modeling (see your paper’s bibliography section).

---

## License

This project is licensed under the MIT License. See `LICENSE` for details.
