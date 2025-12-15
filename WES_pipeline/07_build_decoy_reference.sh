#!/bin/bash
# 07_build_decoy_reference.sh
#
# Run INTERACTIVELY (not DSQ/array).
#
# Purpose:
#   Build a GRCh38 + decoy fasta for BWA indexing by concatenating:
#     - common abundant/contaminant sequences (adapters, rDNA, phix, polyA/C, etc.)
#     - the primary GRCh38 genome.fa
#
# Notes:
#   - BWA index uses `-b 300000000` to handle large references (~300M chars).
#   - request a bit more for safety.

set -euo pipefail

# --- (1) Start an interactive session (example) ---
# srun --pty -c 1 -p interactive --mem=8G -t 2:00:00 bash

# --- (2) Choose a working directory for the decoy reference ---
# (Change this to wherever you keep references for the project.)
WORKDIR="$PWD/decoy_reference"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

# --- (3) Build decoy.fa (abundant sequences + genome.fa) ---
ABUNDANT_DIR="/gpfs/ysm/datasets/genomes/Homo_sapiens/NCBI/GRCh38/Sequence/AbundantSequences"
GENOME_FA="/gpfs/ysm/datasets/genomes/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/genome.fa"

# Abundant/contaminant fasta files to include
# (These are expected to exist in $ABUNDANT_DIR)
cat \
  "$ABUNDANT_DIR/adapter_contam1.fa" \
  "$ABUNDANT_DIR/hum5SrDNA.fa" \
  "$ABUNDANT_DIR/humRibosomal.fa" \
  "$ABUNDANT_DIR/phix.fa" \
  "$ABUNDANT_DIR/polyA.fa" \
  "$ABUNDANT_DIR/polyC.fa" \
  > decoy.fa

# Append GRCh38 genome
cat "$GENOME_FA" >> decoy.fa

# Optional sanity check
echo "decoy.fa created at: $WORKDIR/decoy.fa"
echo "Number of sequences:"
grep -c "^>" decoy.fa || true

# --- (4) Build BWA index ---
module load BWA

# You noted: ~4.6GB for ~300,000,000 characters; use -b 300000000 and bwtsw
bwa index -b 300000000 -a bwtsw decoy.fa

echo "BWA index complete:"
ls -lh decoy.fa* 2>/dev/null || true
