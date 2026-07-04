#!/bin/bash

# Download script for GSE212606 public GEO supplementary files.
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212606
# ftp: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE212nnn/GSE212606/suppl/

set -e

source "$(dirname "$0")/../functions/utils.sh"

DATA_DIR="../../data/BrainOmicsData/raw/GSE212606"

log_message "Starting GSE212606 data download..."

DOWNLOAD_LIST="
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE212nnn/GSE212606/suppl/GSE212606_RAW.tar|GSE212606_RAW.tar|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE212nnn/GSE212606/suppl/GSE212606_processed_data_file_descriptions.txt.gz|GSE212606_processed_data_file_descriptions.txt.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE212nnn/GSE212606/suppl/filelist.txt|filelist.txt|0
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

log_message "GSE212606 data download completed!" --message-type success
