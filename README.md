# Post-liquefaction recovery screen

Public reproducibility package for the manuscript:

> Post-liquefaction recovery is too slow to lower re-liquefaction hazard within an earthquake sequence: a demand-referenced screen and field test from Canterbury

Author: Lijian REN, Inner Mongolia University of Technology / Hohai University.

This repository contains the public code, derived tables, generated figures, and effective-stress simulation support files for a demand-referenced post-liquefaction recovery screen and a Canterbury sequence field test. The manuscript is under submission preparation; this repository does not contain active manuscript files, cover letters, review correspondence, or raw third-party data.

## Contents

- `analysis/`: Python analysis and figure-generation scripts.
- `opensees_campaign/`: OpenSees/PM4Sand effective-stress simulation scripts and small derived simulation outputs used for the simulation figure.
- `outputs/`: derived non-sensitive tables and module summaries used to document the analysis chain.
- `figures/`: generated figure exports in PDF/PNG.
- `data/parquet_coord_bridge.csv`: derived coordinate bridge table used by the open pipeline.
- `DATASETS_AND_LINKS.csv`: source registry for external data that must be downloaded separately.
- `REPRODUCIBLE_RUNBOOK.md`: environment, data placement, and run-order notes.

## Data boundary

Large or third-party inputs are not redistributed here. Download or query them from their original public sources, then place them under `data/` following `REPRODUCIBLE_RUNBOOK.md`.

The analysis-ready Canterbury CPTu Parquet input is archived at:

- https://doi.org/10.5281/zenodo.20839217

Mapped liquefaction polygons and USGS ShakeMap products are public sources but are not copied into this repository.

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python analysis/fig1_mechanism.py
python analysis/fig2_evidence.py
python analysis/graphical_abstract.py
python opensees_campaign/fig_simulation.py
```

The figure scripts above regenerate the included submitted figure exports from fixed values and derived outputs. Full analytical reruns require the external data listed in `DATASETS_AND_LINKS.csv`; rerunning the nonlinear OpenSees sweeps also requires `openseespy`.

## Citation

Use `CITATION.cff` for software citation metadata. Cite the input CPTu dataset separately using DOI `10.5281/zenodo.20839217`. Versioned code releases are archived from GitHub through Zenodo when the release workflow is enabled.

## License

Code in this repository is released under the MIT License. Third-party input data are governed by their own source licenses and access terms. Derived tables and figures are provided for research reproducibility; do not treat them as a redistribution of the external input datasets.
