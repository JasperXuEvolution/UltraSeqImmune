#!/bin/bash
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=Immunoediting_HI_BT_NMethod_CE_N100_R10000_K
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=5g 
#SBATCH --time=12:00:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=batch

source ../../../config.sh

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate UltraSeqImmune



# command
mkdir -p "$BT_ANALYSIS/Output_data/K_mice/"

python3 "$BT_ANALYSIS/Python_scripts/UltraSeq_Boostrapping_Immunoediting.py" \
--a0 "$BT_ANALYSIS/Input_data/Immunoediting_HI_raw_final_df.parquet" \
--a2 100 --a3 100 --a4 10000 --a5 KT --a6 KT --a7 10 \
--o1 "$BT_ANALYSIS/Output_data/K_mice/Immunoediting_HI_BT" \
--o2 "$BT_ANALYSIS/Output_data/K_mice/Immunoediting_HI_BT" \
--l1 30 40 50 60 70 80 90 95 96 97 98 99 \
--m 'N' --c 'Yes'

sacct --format=JobID,JobName,Submit,Start,End,State,Partition,ReqTRES%30,CPUTime,MaxRSS,NodeList%30 --units=M -j $SLURM_JOBID