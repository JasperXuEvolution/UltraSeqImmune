# Immunoediting — Ultra-Seq analysis

This repository contains data analysis scripts for the study **“Tumor suppressor genotype influences the extent and mode of immunosurveillance in lung cancer”**.

The code implements an **Ultra-Seq** workflow: paired-end reads are reduced to a **gRNA × clonal barcode × sample** table, then cleaned and summarized, then **bootstrapped** to quantify uncertainty in tumor metrics, with **contrasts between lentiviral vectors** (different immune microenvironments, e.g. NI / MI / HI).

---

## 1. Barcode extraction

* See **[`01_data_collection/`](01_data_collection/)** — [`README`](01_data_collection/README.md).
* Main pipeline: [`01-info_extraction.bash`](01_data_collection/01-info_extraction.bash).
* **Vector ID**, **sgRNA** and **clonal barcode** are extracted; Python helpers are **`UltraSeq_Step2.py`**, **`UltraSeq_Step3.py`**, **`UltraSeq_Step5.py`** under [`01_data_collection/Python_scripts/`](01_data_collection/Python_scripts/).
* **Output:** a combined table (e.g. `gRNA_clonalbarcode_combined.csv`) with sgRNA, clonal barcode, sample ID, and read counts.

---

## 2. QC and preprocessing

* See **[`02_data_cleaning_and_QC/`](02_data_cleaning_and_QC/)** — [`README`](02_data_cleaning_and_QC/README.md).
* **Step 1 —** [`01-data-preprocessing.ipynb`](02_data_cleaning_and_QC/01-data-preprocessing.ipynb): load the combined barcode table and metadata, filter and aggregate, and write per-vector tables and summaries under **`02_data_cleaning_and_QC/data/`**.
* **Step 2 —** [`02-QC-plotting.ipynb`](02_data_cleaning_and_QC/02-QC-plotting.ipynb): QC figures and exploratory analysis (uses [`UltraSeq_QC_functions.py`](02_data_cleaning_and_QC/UltraSeq_QC_functions.py)).
* **Step 3 —** [`03-input_data_for_BT.ipynb`](02_data_cleaning_and_QC/03-input_data_for_BT.ipynb): build **parquet** inputs for bootstrapping (`Immunoediting_*_raw_final_df.parquet` in **`data/`**), passed to [`03_bootstrapping/`](03_bootstrapping/).

---

## 3. Bootstrapping and vector contrasts

* See **[`03_bootstrapping/`](03_bootstrapping/)** — [`README`](03_bootstrapping/README.md).
* **Step A:** [`UltraSeq_Boostrapping_Immunoediting.py`](03_bootstrapping/Python_scripts/UltraSeq_Boostrapping_Immunoediting.py) — bootstrap resampling for cohort-normalized sgRNA/gene metrics.
* **Step B:** [`Immunoediting_BT_TreatmentEffect.py`](03_bootstrapping/Python_scripts/Immunoediting_BT_TreatmentEffect.py) — compare **two** bootstrap outputs (reference vs contrast vector).
* Example Slurm jobs: [`03_bootstrapping/Job_submission/`](03_bootstrapping/Job_submission/).

---

## 4. Figures

* See **[`04_plotting/`](04_plotting/)** — notebooks `Immunoediting-figures_final_p1.ipynb` … `p3.ipynb` and example tables under **`04_plotting/data/`**.

---

## Config and environment

* **[`config.sh`](config.sh)**.
* **Conda:** [`UltraSeqImmune.yml`](UltraSeqImmune.yml).
