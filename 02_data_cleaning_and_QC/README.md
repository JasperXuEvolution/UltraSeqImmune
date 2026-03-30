# 02 — Data cleaning and QC

This folder contains **notebooks** to turn raw combined tables from data collection into cleaned, analysis-ready objects, then **QC plots**.

**Workflow (high level):**

| Step | Notebook | Role |
|------|----------|------|
| 1 | [`01-data-preprocessing.ipynb`](01-data-preprocessing.ipynb) | Load `gRNA_clonalbarcode_combined.csv`, metadata, and gRNA info; filter and aggregate; write per-vector and summary outputs under **`data/`**. |
| 2 | [`02-QC-plotting.ipynb`](02-QC-plotting.ipynb) | QC figures and exploratory analysis on the preprocessed data. |
