#!/bin/bash
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=Immunoediting_MIvsNI_Normal_R10000_gene
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=50g 
#SBATCH --time=4:00:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=batch

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${SCRIPT_DIR}/../../../config.sh"
[[ -f "$CONFIG" ]] || { echo "Missing $CONFIG" >&2; exit 1; }
# shellcheck source=/dev/null
source "$CONFIG"

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate UltraSeqImmune



# command
mkdir -p "$BT_ANALYSIS/Output_data/GenotypeEffect/"

python3 "$BT_ANALYSIS/Python_scripts/Immunoediting_BT_TreatmentEffect.py" \
--a0 "$BT_ANALYSIS/Output_data/Normal_Method/Immunoediting_NI_BT_NormalMethod_KTC_N100_R10000_gene_intermediate" \
--a1 "$BT_ANALYSIS/Output_data/Normal_Method/Immunoediting_MI_BT_NormalMethod_KTC_N100_R10000_gene_intermediate" \
--a2 Targeted_gene_name \
--o1 "$BT_ANALYSIS/Output_data/GenotypeEffect/Immunoediting_MIvsNI_Normal_R10000" 


sacct --format=JobID,JobName%20,Submit,Start,End,State,Partition,CPUTime,MaxRSS --units=G -j $SLURM_JOBID