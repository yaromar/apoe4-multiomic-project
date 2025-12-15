#!/usr/bin/env bash
# ==============================================================================
# 03_quip_to_gz.sh
#
# Purpose
# -------
# Convert Quip-compressed FASTQ files (*.qp) to gzip-compressed FASTQs (*.gz)
# using Quip + pigz, executed as a Slurm array via dSQ.
#
# This script:
#   1) Scans ~/scratch60 for *.qp files
#   2) Writes one command per file into a DSQ joblist
#   3) Uses pre-determined, validated Slurm parameters
#
# Confirmed parameters (from dry run on heaviest samples)
# -------------------------------------------------------
#   CPUs per task:        10
#   Walltime:             10:00:00
#   Max concurrent jobs:  20
#   Array size:           auto (from joblist; e.g. 0–1817)
#
# ==============================================================================

set -euo pipefail

# -----------------------------
# Configuration
# -----------------------------
INPUT_ROOT="$HOME/scratch60"
JOBLIST="joblist_quip_to_gz.txt"
THREADS=10          # pigz threads; must match Slurm -c
REMOVE_QP=1         # delete .qp after successful conversion

echo "[INFO] Input root:  $INPUT_ROOT"
echo "[INFO] Joblist:     $JOBLIST"
echo "[INFO] Threads:     $THREADS"
echo "[INFO] Remove .qp:  $REMOVE_QP"
echo

# Start fresh
: > "$JOBLIST"

# Find all .qp files
mapfile -t QP_FILES < <(find "$INPUT_ROOT" -type f -name "*.qp" | sort)

if [[ "${#QP_FILES[@]}" -eq 0 ]]; then
  echo "[WARN] No .qp files found. Exiting."
  exit 0
fi

echo "[INFO] Found ${#QP_FILES[@]} .qp files"

# Build DSQ joblist
for qp in "${QP_FILES[@]}"; do
  gz="${qp%.qp}.gz"

  if [[ "$REMOVE_QP" -eq 1 ]]; then
    echo "bash -lc 'module load Quip; module load pigz; \
quip -d --stdout \"${qp}\" | pigz -p ${THREADS} > \"${gz}\" && \
test -s \"${gz}\" && rm -f \"${qp}\"'" >> "$JOBLIST"
  else
    echo "bash -lc 'module load Quip; module load pigz; \
quip -d --stdout \"${qp}\" | pigz -p ${THREADS} > \"${gz}\" && \
test -s \"${gz}\"'" >> "$JOBLIST"
  fi
done

echo "[INFO] Joblist written: $JOBLIST"
head -n 3 "$JOBLIST"
echo

# -----------------------------
# Generate DSQ Slurm script
# -----------------------------
module load dSQ

dsq \
  --job-file "$JOBLIST" \
  -c 10 \
  -t 10:00:00 \
  --max-jobs 20 \
  --output DSQ_results

echo
echo "[INFO] DSQ Slurm script generated."
echo "[INFO] Submit with: sbatch dsq-joblist-*.slurm"
