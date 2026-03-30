#!/bin/bash
#SBATCH --job-name=UltraSeqImmune_IE
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=50g
#SBATCH --time=24:00:00
#SBATCH --partition=YOUR_PARTITION
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=YOUR_EMAIL@example.com

set -euo pipefail

source ../config.sh
: "${DATA_COLLECTION:?}"

# Run root (from config) and paths under it
RUN="$DATA_COLLECTION"
PY="${RUN}/Python_scripts"
MANIFEST="${RUN}/NGS_address"
MERGE="${RUN}/Merging"
BT="${RUN}/Bartender"
PROC="${RUN}/Processed_data"
GRNA="${RUN}/gRNA_information.csv"

module load adapterremoval/2.3.1
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate UltraSeqImmune

mkdir -p "$MERGE" "$BT" "$PROC"

while read -r line; do
  [[ -z "${line// }" || "$line" =~ ^# ]] && continue
  r1=$(echo "$line" | cut -d',' -f1)
  r2=$(echo "$line" | cut -d',' -f2)
  sampleID=$(echo "$line" | cut -d',' -f3)

  # Step 1: AdapterRemoval merge
  d1="${MERGE}/${sampleID}"
  mkdir -p "$d1"
  AdapterRemoval --file1 "$r1" --file2 "$r2" \
    --basename "${d1}/Merged" --collapse --gzip
  echo "Merged: ${d1}"

  # Step 2: Bartender input
  d2="${BT}/${sampleID}"
  mkdir -p "$d2"
  python3 "${PY}/UltraSeq_Step2.py" --a "${MERGE}/${sampleID}/Merged.collapsed.gz" --o "$d2"

  # Step 3: sgRNA clustering (bartender + Step3)
  bartender_single_com -z -1 -d 2 -l 5 -f "${d2}/gRNA.bartender" -o "${d2}/gRNA"
  mkdir -p "${BT}/${sampleID}/Clonal_barcode"
  python3 "${PY}/UltraSeq_Step3.py" \
    --a1 "${d2}/gRNA_barcode.csv" --a2 "${d2}/gRNA_cluster.csv" --a3 "$GRNA" \
    --a5 "${d2}/gRNA.bartender" --a6 "${d2}/clonalbarcode.bartender" --o "${d2}/"

  # Step 4: clonal barcode clustering
  while read -r bpath; do
    out="${bpath/.bartender/}"
    bartender_single_com -z 5 -d 1 -l 5 -f "$bpath" -o "$out"
  done < "${BT}/${sampleID}/Bartender_input_address"

  mkdir -p "${PROC}/${sampleID}"
done < "$MANIFEST"

# Step 5: combine samples
python3 "${PY}/UltraSeq_Step5.py" --a "$BT" --o "${PROC}/"
