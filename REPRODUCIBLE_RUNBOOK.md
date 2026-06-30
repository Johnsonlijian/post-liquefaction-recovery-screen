# Reproducible Runbook

## Environment

Tested with Python 3.11 on Windows. Install the Python dependencies:

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## External data placement

This repository includes derived outputs and a small coordinate bridge table. It does not redistribute large or third-party inputs.

Place external inputs as follows when rerunning the full pipeline:

```text
data/
  canterbury_cptu.parquet
  ecan_mapped_liquefaction_sept2010.geojson
  ecan_mapped_liquefaction_feb2011.geojson
  parquet_coord_bridge.csv
  raw/
    public_demand_manifestation_cache_2026-06-08/
  processed/
    nzgd_profile_parse_2026-06-08/
```

Required public sources and access notes are listed in `DATASETS_AND_LINKS.csv`.

## Run order

The final open-data field test is represented by:

```bash
python analysis/mod07_field_scale_density_paradox.py
python analysis/mod08_LSN_controlled_test.py
python analysis/mod09_buffered_validated_test.py
python analysis/fig1_mechanism.py
python analysis/fig2_evidence.py
python analysis/graphical_abstract.py
python opensees_campaign/fig_simulation.py
```

Earlier diagnostic modules (`mod01` to `mod06` and `eda*`) document the development chain and require additional processed repeat-CPT inputs under `data/processed/nzgd_profile_parse_2026-06-08/`. Their derived output tables are included in `outputs/` for auditability, but their raw processed inputs are not redistributed.

The effective-stress simulation support files are in `opensees_campaign/`. The included small CSV, JSON, and NumPy files are derived outputs used by `fig_simulation.py` to regenerate the manuscript simulation figure. Re-running the nonlinear PM4Sand/PDMY02 sweeps requires `openseespy` and is separated from ordinary figure regeneration.

## Expected outputs

- `outputs/mod07_site_table.csv`
- `outputs/mod08_site_table.csv`
- `outputs/mod09_site_table.csv`
- `outputs/mod09_SUMMARY.md`
- `figures/Fig1_mechanism.pdf`
- `figures/Fig1_mechanism.png`
- `figures/Fig2_evidence.pdf`
- `figures/Fig2_evidence.png`
- `figures/GraphicalAbstract.pdf`
- `figures/GraphicalAbstract.png`
- `opensees_campaign/Fig_simulation.pdf`
- `opensees_campaign/Fig_simulation.png`

## Reproducibility notes

- Random sampling is not used in the final deterministic modules.
- Coordinate and manifestation joins use public identifiers and fixed tolerance choices described in the scripts.
- USGS and public ArcGIS services can change URL structure or service metadata; if a live source has moved, update `DATASETS_AND_LINKS.csv` and rerun the relevant data-fetch step without changing the derived-output filenames.
- Versioned repository archives should be minted from tagged GitHub releases through Zenodo or another archive.
