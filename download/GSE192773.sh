#!/bin/bash

# Download script for GSE192773 data
# paper:
#   https://doi.org/10.1038/s41586-023-06338-4
#   https://doi.org/10.1038/s41467-025-60665-w
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE192773
# code:

set -e

source "$(dirname "$0")/../functions/utils.sh"

DATA_DIR="../../data/BrainOmicsData/raw/GSE192773"

log_message "Starting GSE192773 data download..."

DOWNLOAD_LIST="
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE192773&format=file&file=GSE192773%5FAll%5FSpecies%5FUMI%5FCount%5FMatrix%5FBarcodes%2Etxt%2Egz|GSE192773_All_Species_UMI_Count_Matrix_Barcodes.txt.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE192773&format=file&file=GSE192773%5FAll%5FSpecies%5FUMI%5FCount%5FMatrix%5FGenes%2Etxt%2Egz|GSE192773_All_Species_UMI_Count_Matrix_Genes.txt.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE192773&format=file&file=GSE192773%5FAll%5FSpecies%5FUMI%5FCount%5FMatrix%5FMetaData%2Etxt%2Egz|GSE192773_All_Species_UMI_Count_Matrix_MetaData.txt.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE192773&format=file&file=GSE192773%5FAll%5FSpecies%5FUMI%5FCount%5FMatrix%5FsnRNAseq%2Emtx%2Egz|GSE192773_All_Species_UMI_Count_Matrix_snRNAseq.mtx.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE192773&format=file&file=GSE192773%5FSeq%5FTemplate%5FEvolutionPaper%5FRNA%5FUpdated%2Exlsx|GSE192773_Seq_Template_EvolutionPaper_RNA_Updated.xlsx|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE192773&format=file&file=GSE192773%5FSeurat%5FAstrocyte%5FRNA%2ERDS%2Egz|GSE192773_Seurat_Astrocyte_RNA.RDS.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE192773&format=file&file=GSE192773%5FSeurat%5FExcitatory%5FRNA%2ERDS%2Egz|GSE192773_Seurat_Excitatory_RNA.RDS.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE192773&format=file&file=GSE192773%5FSeurat%5FInhibitory%5FRNA%2ERDS%2Egz|GSE192773_Seurat_Inhibitory_RNA.RDS.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE192773&format=file&file=GSE192773%5FSeurat%5FMicroglia%5FRNA%2ERDS%2Egz|GSE192773_Seurat_Microglia_RNA.RDS.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE192773&format=file&file=GSE192773%5FSeurat%5FOPCOli%5FRNA%2ERDS%2Egz|GSE192773_Seurat_OPCOli_RNA.RDS.gz|0
"

# Function to decompress gzip files
decompress_gz_file() {
    local gz_file="$1"
    local decompressed_file="${gz_file%.gz}"
    
    # Check if compressed file exists
    if [ ! -f "$gz_file" ]; then
        log_message "Compressed file {.file $(basename "$gz_file")} not found, skipping decompression" --message-type warning
        return 1
    fi
    
    # Check if decompressed file already exists
    if [ -f "$decompressed_file" ]; then
        log_message "Decompressed file {.file $(basename "$decompressed_file")} already exists, skipping decompression"
        return 0
    fi
    
    # Decompress the file
    log_message "Decompressing {.file $(basename "$gz_file")}..."
    if gunzip -c "$gz_file" > "$decompressed_file"; then
        local size=$(stat -f%z "$decompressed_file" 2>/dev/null || echo "0")
        local size_mb=$((size / 1024 / 1024))
        log_message "Successfully decompressed {.file $(basename "$decompressed_file")} ({.val ${size_mb}}MB)" --message-type success
        return 0
    else
        log_message "Failed to decompress {.file $(basename "$gz_file")}" --message-type error || true
        return 1
    fi
}

# Function to extract compressed archives (keeps original files)
extract_compressed_file() {
    local file="$1"
    local dir
    dir="$(dirname "$file")"
    
    if [ ! -f "$file" ]; then
        log_message "Compressed file {.file $(basename "$file")} not found, skipping extraction" --message-type warning
        return 1
    fi
    
    case "$file" in
        *.tar.gz|*.tgz)
            log_message "Extracting {.file $(basename "$file")}..."
            if tar -xzf "$file" -C "$dir" -k; then
                log_message "Successfully extracted {.file $(basename "$file")}" --message-type success
                return 0
            else
                log_message "Failed to extract {.file $(basename "$file")}" --message-type error || true
                return 1
            fi
            ;;
        *.tar)
            log_message "Extracting {.file $(basename "$file")}..."
            if tar -xf "$file" -C "$dir" -k; then
                log_message "Successfully extracted {.file $(basename "$file")}" --message-type success
                return 0
            else
                log_message "Failed to extract {.file $(basename "$file")}" --message-type error || true
                return 1
            fi
            ;;
        *.zip)
            log_message "Extracting {.file $(basename "$file")}..."
            if unzip -n -q "$file" -d "$dir"; then
                log_message "Successfully extracted {.file $(basename "$file")}" --message-type success
                return 0
            else
                log_message "Failed to extract {.file $(basename "$file")}" --message-type error || true
                return 1
            fi
            ;;
        *.gz)
            decompress_gz_file "$file"
            ;;
        *)
            return 0
            ;;
    esac
}

batch_download_parallel "$DOWNLOAD_LIST" "$DATA_DIR" 5 3

log_message "Starting extraction of compressed files..."
compressed_count=0
extracted_count=0

while IFS='|' read -r url filename expected_size; do
    [[ -z "$url" || "$url" =~ ^[[:space:]]*# ]] && continue
    
    if [[ "$filename" == *.tar.gz || "$filename" == *.tgz || "$filename" == *.tar || "$filename" == *.zip || "$filename" == *.gz ]]; then
        compressed_count=$((compressed_count + 1))
        compressed_file="$DATA_DIR/$filename"
        if extract_compressed_file "$compressed_file"; then
            extracted_count=$((extracted_count + 1))
        fi
    fi
done <<< "$DOWNLOAD_LIST"

log_message "Extraction completed: {.val ${extracted_count}}/{.val ${compressed_count}} files extracted"


cleanup_temp_files "$DATA_DIR"

log_message "GSE192773 data download and extraction completed!" --message-type success
