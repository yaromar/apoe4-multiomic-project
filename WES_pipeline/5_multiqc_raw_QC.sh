#!/bin/bash
# 05_multiqc_raw_QC.sh
#
# Purpose:
#   Run MultiQC interactively to aggregate *raw* FastQC outputs (and any other QC logs)
#   into a single report.
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
#   1 CPU, 15G RAM, ~5 minutes walltime 

set -euo pipefail

OUTDIR="multiqc_report"
mkdir -p "${OUTDIR}"

# MultiQC scans the current directory recursively by default when given "*".
# If you want to be explicit, you can replace "*" with ".".
multiqc * -o "${OUTDIR}"

echo "MultiQC report written to: ${OUTDIR}/"
