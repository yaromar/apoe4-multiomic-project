#!/bin/bash
# 14_apply_bqsr.sh
#
# Purpose:
#   Apply GATK BQSR to each sample BAM and remove the pre-BQSR BAM + index.
#
# DSQ usage:
#   1) Generate joblist:
#        bash 14_apply_bqsr.sh > joblist.txt
#   2) Run DSQ (parameters already validated):
#        module load dSQ
#        dsq --job-file joblist.txt \
#            -c 1 -t 7:00:00 \
#            --max-jobs 60 \
#            --output DSQ_results
#
# Notes:
#   - Assumes recal table is named *_recal_data.table next to the BAM.

REF=~/project/AD_WES_yaro/reference_genome/decoy.fa

for folder in ~/scratch60/*; do
  # Expect exactly one BAM per folder at this stage
  for bam in "$folder"/*.bam; do
    # Skip if glob doesn't match
    [[ -e "$bam" ]] || continue

    recal="${bam%.bam}_recal_data.table"
    out="${bam%.bam}_bqsr.bam"
    bai="${bam%.bam}.bai"

    # Only emit a command if recal table exists (prevents DSQ jobs that immediately fail)
    # If you'd rather fail loudly, remove this check.
    echo "if [[ -f \"$recal\" ]]; then \
module load GATK; \
gatk ApplyBQSR -R \"$REF\" -I \"$bam\" --bqsr-recal-file \"$recal\" -O \"$out\" && \
rm -f \"$bam\" \"$bai\"; \
else echo \"Missing recal table: $recal\" 1>&2; exit 1; fi"
  done
done
