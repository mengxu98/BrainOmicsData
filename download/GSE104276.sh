#!/bin/bash

# Download script for GSE104276 public GEO supplementary files.
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104276
# ftp: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104276/suppl/

set -e

source "$(dirname "$0")/../functions/utils.sh"

DATA_DIR="../../data/BrainOmicsData/raw/GSE104276"

log_message "Starting GSE104276 data download..."

DOWNLOAD_LIST="
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104276/suppl/GSE104276_RAW.tar|GSE104276_RAW.tar|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104276/suppl/GSE104276_all_pfc_2394_UMI_TPM_NOERCC.xls.gz|GSE104276_all_pfc_2394_UMI_TPM_NOERCC.xls.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104276/suppl/GSE104276_all_pfc_2394_UMI_count_NOERCC.xls.gz|GSE104276_all_pfc_2394_UMI_count_NOERCC.xls.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104276/suppl/GSE104276_readme_sample_barcode.xlsx|GSE104276_readme_sample_barcode.xlsx|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104276/suppl/filelist.txt|filelist.txt|0
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

log_message "GSE104276 data download completed!" --message-type success
