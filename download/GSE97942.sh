#!/bin/bash

# Download script for GSE97942 data
# paper: https://doi.org/10.1038/nbt.4038
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/29227469
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97942


set -e

source "$(dirname "$0")/../functions/utils.sh"

DATA_DIR="../../data/BrainOmicsData/raw/GSE97942"

log_message "Starting GSE97942 data download..."

DOWNLOAD_LIST="
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE97930&format=file&file=GSE97930%5FCerebellarHem%5FsnDrop%2Dseq%5FUMI%5FCount%5FMatrix%5F08%2D01%2D2017%2Etxt%2Egz|GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE97930&format=file&file=GSE97930%5FFrontalCortex%5FsnDrop%2Dseq%5FUMI%5FCount%5FMatrix%5F08%2D01%2D2017%2Etxt%2Egz|GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE97930&format=file&file=GSE97930%5FVisualCortex%5FsnDrop%2Dseq%5FUMI%5FCount%5FMatrix%5F08%2D01%2D2017%2Etxt%2Egz|GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE97930&format=file&file=GSE97930%5FREADME%2Etxt|README.txt|0
"

# Perform batch download
batch_download "$DOWNLOAD_LIST" "$DATA_DIR" 5

# Extract and rename files
if [ ! -f "$DATA_DIR/.GSE97942_extracted" ]; then
    log_message "Extracting and renaming count matrix files..."
    mkdir -p "$DATA_DIR"
    
    # Extract CerebellarHem
    if [ -f "$DATA_DIR/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz" ]; then
        log_message "Extracting CerebellarHem counts..."
        gunzip -c "$DATA_DIR/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz" > "$DATA_DIR/CerebellarHem_counts.txt"
    fi
    
    # Extract FrontalCortex
    if [ -f "$DATA_DIR/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz" ]; then
        log_message "Extracting FrontalCortex counts..."
        gunzip -c "$DATA_DIR/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz" > "$DATA_DIR/FrontalCortex_counts.txt"
    fi
    
    # Extract VisualCortex
    if [ -f "$DATA_DIR/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz" ]; then
        log_message "Extracting VisualCortex counts..."
        gunzip -c "$DATA_DIR/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz" > "$DATA_DIR/VisualCortex_counts.txt"
    fi
    
    touch "$DATA_DIR/.GSE97942_extracted"
    log_message "GSE97942 data extraction completed!" --message-type success
else
    log_message "Files already extracted (marker exists), skipping extraction."
fi


cleanup_temp_files "$DATA_DIR"

log_message "GSE97942 data download and organization completed!" --message-type success
