# 02 — Data cleaning and QC

This folder contains **notebooks** that turn the combined barcode table from data collection into cleaned, analysis-ready tables, **QC plots**, and **bootstrap input** files.

**Workflow (high level):**

| Step | Notebook | Role |
|------|----------|------|
| 1 | [`01-data-preprocessing.ipynb`](01-data-preprocessing.ipynb) | Load `gRNA_clonalbarcode_combined.csv`, metadata, and gRNA info; filter and aggregate; write per-vector and summary outputs under **`data/`** (e.g. `Immunoediting_*_final_df.csv`, `Immunoediting_sample_summary_df.csv`). |
| 2 | [`02-QC-plotting.ipynb`](02-QC-plotting.ipynb) | QC figures and exploratory analysis on the preprocessed data (uses [`UltraSeq_QC_functions.py`](UltraSeq_QC_functions.py)). |
| 3 | [`03-input_data_for_BT.ipynb`](03-input_data_for_BT.ipynb) | From the NI / MI / HI final tables, build **parquet** inputs for bootstrapping (`Immunoediting_*_raw_final_df.parquet` under **`data/`**), used by [`03_bootstrapping/`](../03_bootstrapping/). |

**Shared code:** [`UltraSeq_QC_functions.py`](UltraSeq_QC_functions.py) — QC helpers used by the QC notebook.

**Environment:** use the project conda env from [`UltraSeqImmune.yml`](../UltraSeqImmune.yml) (see root [`README.md`](../README.md)).
