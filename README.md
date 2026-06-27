# Post-liquefaction recovery screen

Public reproducibility package for the manuscript:

> Post-liquefaction recovery is too slow to lower re-liquefaction hazard within an earthquake sequence: a demand-referenced screen and field test from Canterbury

Author: Lijian REN, Inner Mongolia University of Technology / Hohai University.

This repository contains the public code, derived tables, and generated figures for a demand-referenced post-liquefaction recovery screen and a Canterbury sequence field test. The manuscript is under submission preparation; this repository does not contain active manuscript files, cover letters, reviewer material, or raw third-party data.

## Contents

- `analysis/`: Python analysis and figure-generation scripts.
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
python analysis/fig_converged.py
python analysis/graphical_abstract.py
```

The figure scripts above regenerate the included figure exports from fixed values and derived outputs. Full analytical reruns require the external data listed in `DATASETS_AND_LINKS.csv`.

## Citation

Use `CITATION.cff` for software citation metadata. Cite the input CPTu dataset separately using DOI `10.5281/zenodo.20839217`. A code archive DOI should be minted from a GitHub release through Zenodo once repository archiving is enabled.

## License

Code in this repository is released under the MIT License. Third-party input data are governed by their own source licenses and access terms. Derived tables and figures are provided for research reproducibility; do not treat them as a redistribution of the external input datasets.
