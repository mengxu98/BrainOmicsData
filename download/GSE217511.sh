#!/bin/bash

# Download script for GSE217511 data
# paper: https://doi.org/10.1038/s41467-022-34975-2
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/36509746
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217511


set -e

source "$(dirname "$0")/../functions/utils.sh"

DATA_DIR="../../data/BrainOmicsData/raw/GSE217511"

log_message "Starting GSE217511 data download..."

DOWNLOAD_LIST="
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE217511&format=file&file=GSE217511%5FCombined%5Flog%2Dnormalized%5FAveExp%5FsnRNAseq%5FHuman%2Exlsx|GSE217511_Combined_log-normalized_AveExp_snRNAseq_Human.xlsx|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE217511&format=file&file=GSE217511%5FCortex%5FSeuratmetadata%2Ecsv%2Egz|GSE217511_Cortex_Seuratmetadata.csv.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE217511&format=file&file=GSE217511%5FCortex%5FUMAPembedding%2Ecsv%2Egz|GSE217511_Cortex_UMAPembedding.csv.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE217511&format=file&file=GSE217511%5FCorticalPlate%5FSeuratmetadata%2Ecsv%2Egz|GSE217511_CorticalPlate_Seuratmetadata.csv.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE217511&format=file&file=GSE217511%5FCorticalPlate%5FUMAPembedding%2Ecsv%2Egz|GSE217511_CorticalPlate_UMAPembedding.csv.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE217511&format=file&file=GSE217511%5FGerminalMatrix%5FSeuratmetadata%2Ecsv%2Egz|GSE217511_GerminalMatrix_Seuratmetadata.csv.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE217511&format=file&file=GSE217511%5FGerminalMatrix%5FUMAPembedding%2Ecsv%2Egz|GSE217511_GerminalMatrix_UMAPembedding.csv.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE217511&format=file|GSE217511_RAW.tar|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE217511&format=file&file=GSE217511%5FSVZ%5FCaudate%5FSeuratmetadata%2Ecsv%2Egz|GSE217511_SVZ_Caudate_Seuratmetadata.csv.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE217511&format=file&file=GSE217511%5FSVZ%5FCaudate%5FUMAPembedding%2Ecsv%2Egz|GSE217511_SVZ_Caudate_UMAPembedding.csv.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE217511&format=file&file=GSE217511%5FsubCorticalPlate%5FSeuratmetadata%2Ecsv%2Egz|GSE217511_subCorticalPlate_Seuratmetadata.csv.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE217511&format=file&file=GSE217511%5FsubCorticalPlate%5FUMAPembedding%2Ecsv%2Egz|GSE217511_subCorticalPlate_UMAPembedding.csv.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE217511&format=file&file=GSE217511%5FsubGerminalMatrix%5FSeuratmetadata%2Ecsv%2Egz|GSE217511_subGerminalMatrix_Seuratmetadata.csv.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE217511&format=file&file=GSE217511%5FsubGerminalMatrix%5FUMAPembedding%2Ecsv%2Egz|GSE217511_subGerminalMatrix_UMAPembedding.csv.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217511/matrix/GSE217511_series_matrix.txt.gz|GSE217511_series_matrix.txt.gz|0
"

# Perform batch download
batch_download "$DOWNLOAD_LIST" "$DATA_DIR" 5

# Extract the tar file to GSE217511/GSE217511 directory
TARGET_DIR="$DATA_DIR/GSE217511"
if [ ! -f "$DATA_DIR/.GSE217511_extracted" ]; then
    log_message "Extracting $DATA_DIR/GSE217511_RAW.tar to $TARGET_DIR ..."
    mkdir -p "$TARGET_DIR"
    tar -xvf "$DATA_DIR/GSE217511_RAW.tar" -C "$TARGET_DIR"
    
    # Organize files by sample (GSM ID)
    log_message "Organizing files by sample (GSM ID)..."
    # Save current directory
    ORIGINAL_DIR=$(pwd)
    cd "$TARGET_DIR"
    
    # Find all unique GSM IDs from file names
    for file in *.gz; do
        if [ -f "$file" ]; then
            # Extract GSM ID (e.g., GSM6720852 from GSM6720852_9C_barcodes.tsv.gz)
            gsm_id=$(echo "$file" | sed -n 's/^\(GSM[0-9]*\)_.*/\1/p')
            if [ -n "$gsm_id" ]; then
                # Create directory for this GSM ID if it doesn't exist
                mkdir -p "$gsm_id"
                # Move file to GSM directory
                mv "$file" "$gsm_id/" 2>/dev/null || true
            fi
        fi
    done
    
    # Return to original directory
    cd "$ORIGINAL_DIR"
    
    # Create a marker so we don't re-extract next runs
    touch "$DATA_DIR/.GSE217511_extracted"
    log_message "Files organized by GSM ID" --message-type success
else
    log_message "Archive already extracted (marker exists), skipping tar extraction."
fi

# Extract all .csv.gz files
log_message "Extracting .csv.gz files..."
for gzfile in "$DATA_DIR"/*.csv.gz; do
    if [ -f "$gzfile" ]; then
        csvfile="${gzfile%.gz}"
        if [ ! -f "$csvfile" ]; then
            log_message "Decompressing $(basename "$gzfile") -> $(basename "$csvfile")"
            gunzip -c "$gzfile" > "$csvfile"
        else
            log_message "File $(basename "$csvfile") already exists, skipping..."
        fi
    fi
done


cleanup_temp_files "$DATA_DIR"

log_message "GSE217511 data download and organization completed!" --message-type success
