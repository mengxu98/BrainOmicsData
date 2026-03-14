#!/bin/bash

# Download script for GSE204682 data
# GSE204683 (Multiome: snRNA-seq + snATAC-seq (GSE204684))
# paper: https://doi.org/10.1126/sciadv.adg3754
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/37824614
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE204683
# ATAC-seq: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE204684
# CELLxGENE (RRID: SCR_021059) data (h5ad):
# https://cellxgene.cziscience.com/collections/ceb895f4-ff9f-403a-b7c3-187a9657ac2c
# code: https://doi.org/10.5281/zenodo.7703253


set -e

source "$(dirname "$0")/../functions/utils.sh"

DATA_DIR="../../data/BrainOmicsData/raw/GSE204682"

log_message "Starting GSE204682 data download..."

DOWNLOAD_LIST="
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE204682&format=file&file=GSE204682%5Fbarcodes%2Etsv%2Egz|GSE204682_barcodes.tsv.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE204682&format=file&file=GSE204682%5Fcount%5Fmatrix%2ERDS%2Egz|GSE204682_count_matrix.RDS.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE204682&format=file&file=GSE204682%5Fpeaks%2Ebed%2Egz|GSE204682_peaks.bed.gz|0
"
# https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE204682&format=file&file=GSE204682%5Fbigwigs%5Fper%5Fcell%5Ftype%2Etar%2Egz|GSE204682_bigwigs_per_cell_type.tar.gz|0
# bigwigs_per_cell_type.tar.gz is a large file (~42.4GB)

# Perform batch download
batch_download "$DOWNLOAD_LIST" "$DATA_DIR" 5

# Extract all .gz files
log_message "Extracting .gz files..."
for gzfile in "$DATA_DIR"/*.gz; do
    if [ -f "$gzfile" ]; then
        outfile="${gzfile%.gz}"
        if [ ! -f "$outfile" ]; then
            log_message "Decompressing $(basename "$gzfile") -> $(basename "$outfile")"
            gunzip -c "$gzfile" > "$outfile"
        else
            log_message "File $(basename "$outfile") already exists, skipping..."
        fi
    fi
done


cleanup_temp_files "$DATA_DIR"

log_message "GSE204682 data download and organization completed!" --message-type success
