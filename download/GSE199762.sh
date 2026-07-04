#!/bin/bash

# Download script for GSE199762 public GEO supplementary files.
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199762
# restricted companion: dbGaP phs003509.v1.p1
# ftp: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE199nnn/GSE199762/suppl/

set -e

source "$(dirname "$0")/../functions/utils.sh"

DATA_DIR="../../data/BrainOmicsData/raw/GSE199762"

log_message "Starting GSE199762 data download..."

DOWNLOAD_LIST="
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE199nnn/GSE199762/suppl/GSE199762_RAW.tar|GSE199762_RAW.tar|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE199nnn/GSE199762/suppl/GSE199762_samples_from_GSE186538.tar.gz|GSE199762_samples_from_GSE186538.tar.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE199nnn/GSE199762/suppl/filelist.txt|filelist.txt|0
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

log_message "GSE199762 data download completed!" --message-type success
