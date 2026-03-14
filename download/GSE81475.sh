#!/bin/bash

# Download script for GSE81475 data
# paper: https://doi.org/10.1038/s41467-022-34975-2
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/36509746
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81475


set -e

source "$(dirname "$0")/../functions/utils.sh"

DATA_DIR="../../data/BrainOmicsData/raw/GSE81475"

log_message "Starting GSE81475 data download..."

DOWNLOAD_LIST="
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE81475&format=file&file=GSE81475%5FZika%2EGEO%2EhumanBrain%2EsingleCell%2Egene%2Ecount%2Etxt%2Egz|GSE81475_Zika.GEO.humanBrain.singleCell.gene.count.txt.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE81475&format=file&file=GSE81475%5FZika%2EGEO%2EhumanBrain%2EsingleCell%2Egene%2ERPKM%2Etxt%2Egz|GSE81475_Zika.GEO.humanBrain.singleCell.gene.RPKM.txt.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81475/matrix/GSE81475_series_matrix.txt.gz|GSE81475_series_matrix.txt.gz|0
"

# Perform batch download
batch_download "$DOWNLOAD_LIST" "$DATA_DIR" 5

# Extract files directly to GSE81475 directory and rename them
if [ ! -f "$DATA_DIR/.GSE81475_extracted" ]; then
    log_message "Extracting $DATA_DIR/GSE81475_Zika.GEO.humanBrain.singleCell.gene.count.txt.gz to $DATA_DIR/counts.txt ..."
    mkdir -p "$DATA_DIR"
    gunzip -c "$DATA_DIR/GSE81475_Zika.GEO.humanBrain.singleCell.gene.count.txt.gz" > "$DATA_DIR/counts.txt"
    log_message "Extracting $DATA_DIR/GSE81475_Zika.GEO.humanBrain.singleCell.gene.RPKM.txt.gz to $DATA_DIR/RPKM.txt ..."
    gunzip -c "$DATA_DIR/GSE81475_Zika.GEO.humanBrain.singleCell.gene.RPKM.txt.gz" > "$DATA_DIR/RPKM.txt"
    touch "$DATA_DIR/.GSE81475_extracted"
    log_message "GSE81475 data extraction completed!" --message-type success
else
    log_message "Archive already extracted (marker exists), skipping extraction."
fi


cleanup_temp_files "$DATA_DIR"

log_message "GSE81475 data download completed!" --message-type success
