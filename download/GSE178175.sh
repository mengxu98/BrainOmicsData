#!/bin/bash

# Download script for GSE178175 public GEO supplementary files.
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178175
# ftp: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178175/suppl/

set -e

source "$(dirname "$0")/../functions/utils.sh"

DATA_DIR="../../data/BrainOmicsData/raw/GSE178175"

log_message "Starting GSE178175 data download..."

DOWNLOAD_LIST="
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178175/suppl/GSE178175_RAW.tar|GSE178175_RAW.tar|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178175/suppl/filelist.txt|filelist.txt|0
"

batch_download "$DOWNLOAD_LIST" "$DATA_DIR" 5

for archive in "$DATA_DIR"/*.tar "$DATA_DIR"/*.tar.gz "$DATA_DIR"/*.tgz; do
    [ -f "$archive" ] || continue
    marker="$DATA_DIR/.${archive##*/}_extracted"
    if [ ! -f "$marker" ]; then
        log_message "Extracting {.file ${archive##*/}}..."
        tar -xf "$archive" -C "$DATA_DIR"
        touch "$marker"
    fi
done

cleanup_temp_files "$DATA_DIR"

log_message "GSE178175 data download completed!" --message-type success
