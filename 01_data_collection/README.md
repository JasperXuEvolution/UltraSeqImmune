# 01 — Data collection

This folder contains the steps needed for **barcode extraction** from sequencing reads: merge paired FASTQs (AdapterRemoval), build Bartender inputs (`UltraSeq_Step2.py`), cluster sgRNA and clonal barcodes (Bartender + `UltraSeq_Step3.py`), then combine samples (`UltraSeq_Step5.py`).

**Output (master table):** Step 5 writes `gRNA_clonalbarcode_combined.csv` with one row per `(gRNA, clonal barcode, sample)` and these columns:

| Column | Meaning |
|--------|--------|
| `gRNA` | sgRNA sequence |
| `Clonal_barcode` | Clonal barcode sequence (10 bp) |
| `Sample_ID` | Sample identifier |
| `Frequency` | Read count for a clonal tumor |
| `Vector_ID` | 4 bp vector tag |

- **Config:** [`../config.sh`](../config.sh) — set **`DATA_COLLECTION`** to the **absolute path of this `01_data_collection` directory** (the folder that contains this README). **`Merging/`**, **`Bartender/`**, and **`Processed_data/`** are **subfolders inside** `01_data_collection` (created by the pipeline; merged FASTQs, intermediate clustering, and final combined tables).
- **Pipeline:** `01-info_extraction.bash` sources that config and runs the loop (Slurm `#SBATCH` lines at top if you use `sbatch`).
- **Conda env:** [`../UltraSeqImmune.yml`](../UltraSeqImmune.yml).

**`NGS_address`**: each line lists **paths to your NGS FASTQs** for one sample—comma-separated **`R1_path,R2_path,Sample_ID`**. The pipeline reads it from **`${DATA_COLLECTION}/NGS_address`**. Example layout: [`NGS_address`](NGS_address) in this folder shows the same format (absolute paths to `.fastq.gz` files).

In that same directory: **`gRNA_information.csv`**, **`NGS_address`**, and **`Python_scripts/`** (with `UltraSeq_Step*.py`) live **next to** `01-info_extraction.bash`, as in this repo layout.
