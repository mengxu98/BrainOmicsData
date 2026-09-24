#!/bin/bash

# Download script for GSE186538 public GEO supplementary files.
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186538
# ftp: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186538/suppl/

set -e

source "$(dirname "$0")/../functions/utils.sh"

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_ROOT="${BRAINOMICS_DATA_ROOT:-$(cd "$REPO_DIR/../.." && pwd)/data/BrainOmicsData}"
DATA_DIR="$DATA_ROOT/raw/GSE186538"

log_message "Starting GSE186538 data download..."

DOWNLOAD_LIST="
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186538/suppl/GSE186538_Human_cell_meta.txt.gz|GSE186538_Human_cell_meta.txt.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186538/suppl/GSE186538_Human_counts.mtx.gz|GSE186538_Human_counts.mtx.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186538/suppl/GSE186538_Human_genes.txt.gz|GSE186538_Human_genes.txt.gz|0
"

batch_download "$DOWNLOAD_LIST" "$DATA_DIR" 5
cleanup_temp_files "$DATA_DIR"

log_message "GSE186538 data download completed!" --message-type success
