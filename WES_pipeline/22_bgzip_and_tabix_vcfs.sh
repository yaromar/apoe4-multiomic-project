#!/bin/bash
# 22_bgzip_and_tabix_vcfs.sh
#
# Purpose:
#   BGZF-compress and tabix-index the VCFs for merging downstream.
#
# Rationale:
#   - bcftools merge/concat and region queries typically expect .vcf.gz (BGZF) + .tbi.
#
# Chosen input pattern:
#   *_filtered_snpeff_dbsnp_clinvar_gwas.vcf  (output of step #20)
#
# Usage pattern (DSQ):
#   bash 22_bgzip_and_tabix_vcfs.sh > joblist_22.txt
#   module load dSQ
#   dsq --job-file joblist_22.txt -c 2 --mem-per-cpu=20G -t 10:00:00 --max-jobs <SET_ME> --output DSQ_results_22
#   sbatch dsq-joblist_22.slurm

for folder in ~/scratch60/*; do
  file=($folder/*filtered_snpeff_dbsnp_clinvar_gwas.vcf)

  # Skip folders without a match
  if [[ ! -e "${file[0]}" ]]; then
    continue
  fi

  vcf="${file[0]}"
  vcfgz="${vcf}.gz"

  # Emit one command per sample for DSQ
  # -f: overwrite if exists (useful if re-running)
  # -p vcf: tabix preset
  echo "module load HTSlib; \
bgzip -f \"$vcf\"; \
tabix -f -p vcf \"$vcfgz\""
done
