#!/bin/bash
# 23_merge_sample_vcfs.sh
#
# Purpose:
#   Merge per-sample VCFs into one cohort-level VCF per brain region (example: CBM).
#
# Expected inputs:
#   BGZF-compressed + tabix-indexed VCFs (*.vcf.gz + *.tbi) produced in step #22.
#
#
# Notes:
#   - If you merged earlier due to header issues, point INPUT_GLOB to that stage
#     (e.g., *filtered_snpeff.vcf.gz). Otherwise, prefer the final annotated stage
#     (e.g., *filtered_snpeff_dbsnp_clinvar_gwas.vcf.gz).

set -euo pipefail

module load BCFtools
module load HTSlib

THREADS=8
OUTDIR="$(pwd)/merged_vcfs"
mkdir -p "$OUTDIR"

# Choose which stage you are merging:
#INPUT_GLOB="*CBM*/*filtered_snpeff.vcf.gz"
INPUT_GLOB="*CBM*/*filtered_snpeff_dbsnp_clinvar_gwas.vcf.gz"

OUT_VCF="$OUTDIR/merged_CBM.vcf.gz"

# Build an explicit file list 
LIST="$OUTDIR/cbm_vcfs.list"
ls -1 $INPUT_GLOB > "$LIST"

echo "Merging $(wc -l < "$LIST") VCFs into: $OUT_VCF"
echo "Input list: $LIST"

# Merge
# -m none: do not merge different alleles into multi-allelic records 
# -Oz: bgzip-compressed output
bcftools merge \
  --threads "$THREADS" \
  -m none \
  -l "$LIST" \
  -Oz -o "$OUT_VCF"

# Index merged VCF
tabix -f -p vcf "$OUT_VCF"

echo "Done."
echo "Merged VCF: $OUT_VCF"
echo "Index: ${OUT_VCF}.tbi"
