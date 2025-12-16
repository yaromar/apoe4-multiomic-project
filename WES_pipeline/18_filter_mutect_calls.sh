#!/bin/bash
# 18_filter_mutect_calls.sh
#
# Purpose:
#   Create a “cleaned” VCF by removing variants that failed Mutect2 filtering
#   (based on FILTER tags / reasons embedded in the VCF records).
#
# What it does:
#   - Writes all header lines (starting with '#') to *_filtered.vcf
#   - Writes only variant records that do NOT match any of the listed filter keywords
#
# DSQ usage:
#   1) Generate joblist:
#        bash 18_filter_mutect_calls.sh > joblist.txt
#   2) Run DSQ (parameters already validated; very light jobs):
#        module load dSQ
#        dsq --job-file joblist.txt \
#            -c 1 --mem-per-cpu=1G -t 0:10:00 \
#            --max-jobs 200 \
#            --output DSQ_results
#
# Resources (validated):
#   - 1 CPU
#   - 1 GB
#   - essentially instant

for folder in ~/scratch60/*; do
  file=("$folder"/*contFiltered.vcf.gz)
  [[ -e "${file[0]}" ]] || continue

  # If multiple matches exist, loop through them
  for vcf_gz in "${file[@]}"; do
    output="${vcf_gz/_mt2_contFiltered.vcf.gz/_filtered.vcf}"

    echo "module load GATK; \
zcat \"$vcf_gz\" | grep '^#' > \"$output\"; \
zcat \"$vcf_gz\" | grep -v '^#' | grep -v -E 'base_qual|clustered_events|contamination|fragment|germline|haplotype|map_qual|multiallelic|orientation|panel_of_normals|position|slippage|strand_bias' >> \"$output\""
  done
done
