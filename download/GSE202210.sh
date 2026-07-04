#!/bin/bash

# Download script for GSE202210 public GEO supplementary files.
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202210
# ftp: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202210/suppl/

set -e

source "$(dirname "$0")/../functions/utils.sh"

DATA_DIR="../../data/BrainOmicsData/raw/GSE202210"

log_message "Starting GSE202210 data download..."

DOWNLOAD_LIST="
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202210/suppl/GSE202210_RAW.tar|GSE202210_RAW.tar|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202210/suppl/filelist.txt|filelist.txt|0
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

log_message "GSE202210 data download completed!" --message-type success
