#!/usr/bin/env bash
# ==============================================================================
# 00_download_fastq.sh
#
# Purpose
# -------
# Download WES FASTQ files from YCGA/FCB HTTP directories using wget.
# Designed for "transfer node" usage (per YCGA note: archived data is only
# accessible from transfer node).
#
# Notes from YCGA (Chris)
# -----------------------
# 1) Lane is denoted in the FASTQ filename/path as L001, L002, L003, L004.
# 2) If multiple lanes exist for a sample and index combo is the same,
#    it is the same library (re-capture + re-sequencing of same prepared library).
# 3) The last 6 digits correspond to I7 and I5 index identifiers used for the sample
#    (e.g., Sample_10260309_ST_BA22_I7_I5). There are 384 possible 10nt sequences
#    numbered 1–384.
#    - Same index combination => re-capture + re-seq of same library.
#    - Different index => library re-prepped from DNA.
#
#
# Usage
# -----
#   # On transfer node:
#   bash 00_download_fastq.sh
#
# Optional:
#   OUTDIR=/path/to/store/data bash 00_download_fastq.sh
#
# Requirements
# ------------
# - bash
# - wget
#
# ==============================================================================

set -euo pipefail

# Where to store downloads (use scratch)
OUTDIR="${OUTDIR:-$PWD/fastq_downloads}"
mkdir -p "$OUTDIR"


WGET_OPTS=(
  "-r"                 # recursive
  "-np"                # no parent
  "-nH"                # no host directories
  "--cut-dirs=2"       # remove leading path components 
  "-R" "index.html*"   # don't download directory index pages
  "--timestamping"     # only download if remote is newer
  "--progress=dot:giga"
)

# YCGA/FCB directories (try only a few at a time -- very heavy)
URLS=(
  "http://fcb.ycga.yale.edu:3010/YXPRL7H9jUWRHX84XOTbx06hlRHcf/sample_dir_000003076/"
  "http://fcb.ycga.yale.edu:3010/5zWewUw0reZ1rA0HCFz9IcO8fCjNV/sample_dir_000003091/"
  "http://fcb.ycga.yale.edu:3010/1QdfSJQMaPf4mbF73rnFrJaasp0ZO/sample_dir_000003260/"
  "http://fcb.ycga.yale.edu:3010/MqwVYi4pvDLiJwKOclruBhVIcORD8/sample_dir_000003343/"
  "http://fcb.ycga.yale.edu:3010/32pBoFUKdqWHYkNVjdYSWP6eZVJZj/sample_dir_000003145/"
  "http://fcb.ycga.yale.edu:3010/f33inzmB21X3Adzb7SFQ8JUW68vhU/sample_dir_000003038/"
)

log() { printf "[%s] %s\n" "$(date '+%F %T')" "$*" >&2; }

log "Starting FASTQ download"
log "Output directory: $OUTDIR"
log "Reminder: run this on the TRANSFER NODE (archived data access requirement)."

cd "$OUTDIR"

for url in "${URLS[@]}"; do
  # Extract the sample_dir name for clearer logs
  sample_dir="$(basename "${url%/}")"
  log "------------------------------------------------------------"
  log "Downloading: $sample_dir"
  log "From URL:    $url"

  # Keep each sample_dir in its own folder to avoid collisions
  mkdir -p "$sample_dir"
  (
    cd "$sample_dir"
    wget "${WGET_OPTS[@]}" "$url" 2>&1 | tee "../${sample_dir}.wget.log"
  )

  log "Finished: $sample_dir"
done

log "All downloads finished."

# ------------------------------------------------------------------------------
# Sanity checks / quick inspection helpers
# ------------------------------------------------------------------------------

log "Sanity check: list FASTQ files sorted by size (largest last)"
# This roughly matches:
#   ls -lR | grep 'fastq.gz' | tr -s ' ' | sort -t' ' -k4 -V
# but is more robust and produces size in bytes.
find . -type f -name "*.fastq.gz" -printf "%s\t%p\n" | sort -n | tail -n 50

log "Sanity check: summarize lanes seen (L001–L004) across downloaded FASTQs"
find . -type f -name "*.fastq.gz" \
  | sed -n 's/.*_L\\([0-9][0-9][0-9]\\)_.*/L\\1/p' \
  | sort | uniq -c | sort -nr

log "Done."
