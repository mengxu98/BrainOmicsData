#!/bin/bash

# Download script for GSE207334 public GEO supplementary files.
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207334
# ftp: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207334/suppl/

set -e

source "$(dirname "$0")/../functions/utils.sh"

DATA_DIR="../../data/BrainOmicsData/raw/GSE207334"

log_message "Starting GSE207334 data download..."

DOWNLOAD_LIST="
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207334/suppl/GSE207334_Multiome_atac_counts.mtx.gz|GSE207334_Multiome_atac_counts.mtx.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207334/suppl/GSE207334_Multiome_atac_peaks.txt.gz|GSE207334_Multiome_atac_peaks.txt.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207334/suppl/GSE207334_Multiome_cell_meta.txt.gz|GSE207334_Multiome_cell_meta.txt.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207334/suppl/GSE207334_Multiome_rna_counts.mtx.gz|GSE207334_Multiome_rna_counts.mtx.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207334/suppl/GSE207334_Multiome_rna_genes.txt.gz|GSE207334_Multiome_rna_genes.txt.gz|0
"

batch_download "$DOWNLOAD_LIST" "$DATA_DIR" 5
cleanup_temp_files "$DATA_DIR"

log_message "GSE207334 data download completed!" --message-type success
