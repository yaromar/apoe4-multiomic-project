#!/bin/bash
# 05_multiqc_raw_QC.sh
#
# Purpose:
#   Run MultiQC interactively to aggregate *raw* FastQC outputs (and any other QC logs)
#   into a single report, then clean up per-file FastQC reports to save space.
#
# This step is NOT a DSQ/array job — run interactively.
#
# Interactive command (validated):
#   srun --pty -c 1 -p interactive --mem-per-cpu=15G bash
#
# Then:
#   module load miniconda
#   conda activate py37_dev
#
# Recommended usage from the directory that contains your FastQC outputs, e.g.:
#   cd ~/scratch60/fastqc_raw
#   bash /path/to/WES_pipeline/05_multiqc_raw_QC.sh
#
# Resources:
#   1 CPU, 15G RAM, ~5 minutes walltime (typical)

set -euo pipefail

OUTDIR="multiqc_report"
mkdir -p "${OUTDIR}"

# Run MultiQC (recursively scans current directory)
multiqc * -o "${OUTDIR}"

# Cleanup: remove individual FastQC reports after aggregation
# This removes both *_fastqc.html and *_fastqc.zip files inside subdirectories
rm -f */*fastqc*

echo "MultiQC report written to: ${OUTDIR}/"
echo "Per-file FastQC reports removed."
