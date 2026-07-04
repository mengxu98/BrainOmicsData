#!/bin/bash

# Download script for GSE103723 public GEO supplementary files.
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103723
# ftp: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103723/suppl/

set -e

source "$(dirname "$0")/../functions/utils.sh"

DATA_DIR="../../data/BrainOmicsData/raw/GSE103723"

log_message "Starting GSE103723 data download..."

DOWNLOAD_LIST="
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103723/suppl/GSE103723_RAW.tar|GSE103723_RAW.tar|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103723/suppl/GSE103723_Region_Sample_Barcode_Information.xlsx|GSE103723_Region_Sample_Barcode_Information.xlsx|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103723/suppl/filelist.txt|filelist.txt|0
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

log_message "GSE103723 data download completed!" --message-type success
