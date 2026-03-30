# 03 — Bootstrapping

This folder runs **bootstrap resampling** on preprocessed tumor tables (from [**02**](../02_data_cleaning_and_QC/README.md)) to get **uncertainty and p-values** for cohort-normalized sgRNA/gene metrics. A second step **contrasts bootstrap outputs between lentiviral vectors** (different **immune microenvironments**, e.g. NI / MI / HI).

**Workflow (high level):**

| Step | What | Role |
|------|------|------|
| 1 | [`Python_scripts/UltraSeq_Boostrapping_Immunoediting.py`](Python_scripts/UltraSeq_Boostrapping_Immunoediting.py) | Read a parquet of gRNA-level tumors; compare focal vs reference genotypes (or plasmid); resample mice; write **intermediate** and **summary** tables under **`Output_data/`**. |
| 2 | [`Python_scripts/Immunoediting_BT_TreatmentEffect.py`](Python_scripts/Immunoediting_BT_TreatmentEffect.py) | Read **two** bootstrap **intermediate** table (one vector as reference, one as contrast); compute **fold effects** and **empirical p-values / FDR** at gRNA or gene level. |

- **Config:** [`../config.sh`](../config.sh) — set **`BT_ANALYSIS`** to the **absolute path of this `03_bootstrapping` directory**. Expect **`Input_data/`** and **`Output_data/`** under that root (parquet inputs from step 02; outputs from the Python scripts).
- **Conda env:** [`../UltraSeqImmune.yml`](../UltraSeqImmune.yml).
