#!/bin/bash

# Download script for GSE235493 data
# paper: https://doi.org/10.1126/sciadv.adu7944
# data: 
#	fetal whole-brain: https://doi.org/10.5281/zenodo.13879662
#	midbrain atlases: https://doi.org/10.5281/zenodo.13879918
# code: https://doi.org/10.5281/zenodo.15243470

set -e

source "$(dirname "$0")/../functions/utils.sh"

DATA_DIR="../../data/BrainOmicsData/raw/BrainSTEM"

log_message "Starting BrainSTEM data download..."

DOWNLOAD_LIST="
https://zenodo.org/records/17291276/files/fetalMidbrainAtlas.noRNA.rds?download=1|fetalMidbrainAtlas.noRNA.rds|0
https://zenodo.org/records/17291276/files/fetalWholeBrainAtlas.noRNA.rds?download=1|fetalWholeBrainAtlas.noRNA.rds|0
https://zenodo.org/api/records/13879918/files-archive|midbrain_atlases.zip|0
"

# Perform batch download
batch_download_parallel "$DOWNLOAD_LIST" "$DATA_DIR" 10 10

# Verify downloaded files
log_message "Verifying downloaded files..."
verification_failed=0

while IFS='|' read -r url filename expected_size; do
    [[ -z "$url" || "$url" =~ ^[[:space:]]*# ]] && continue
    
    file_path="$DATA_DIR/$filename"
    if [ ! -f "$file_path" ]; then
        log_message "File {.file ${filename}} not found after download" --message-type error || true
        verification_failed=$((verification_failed + 1))
        continue
    fi
    
    if ! verify_file_integrity "$file_path" "$expected_size"; then
        verification_failed=$((verification_failed + 1))
    fi
done <<< "$DOWNLOAD_LIST"

if [ $verification_failed -gt 0 ]; then
    log_message "File verification failed for {.val ${verification_failed}} file(s)" --message-type error || true
    exit 1
fi

log_message "All files verified successfully!" --message-type success

# Extract zip files (idempotent)
extract_zip_file() {
    local zip_file="$1"
    local extract_dir="$2"
    local marker_file="$3"
    
    if [ ! -f "$zip_file" ]; then
        log_message "Zip file {.file $(basename "$zip_file")} not found, skipping extraction" --message-type warning
        return 1
    fi
    
    if [ -f "$marker_file" ]; then
        log_message "Archive already extracted (marker exists), skipping extraction: {.file $(basename "$zip_file")}" 
        return 0
    fi
    
    log_message "Extracting {.file $(basename "$zip_file")}..."
    if unzip -q -o "$zip_file" -d "$extract_dir" 2>/dev/null || unzip -o "$zip_file" -d "$extract_dir"; then
        touch "$marker_file"
        log_message "Successfully extracted {.file $(basename "$zip_file")}" --message-type success
        return 0
    else
        log_message "Failed to extract {.file $(basename "$zip_file")}" --message-type error || true
        return 1
    fi
}

# Extract zip archives
while IFS='|' read -r url filename expected_size; do
    [[ -z "$url" || "$url" =~ ^[[:space:]]*# ]] && continue
    
    if [[ "$filename" == *.zip ]]; then
        zip_file="$DATA_DIR/$filename"
        marker_file="$DATA_DIR/.${filename%.zip}_extracted"
        extract_zip_file "$zip_file" "$DATA_DIR" "$marker_file"
    fi
done <<< "$DOWNLOAD_LIST"

# Check extracted data integrity
log_message "Checking extracted data integrity..."
total_files=0
found_files=0

while IFS='|' read -r url filename expected_size; do
    [[ -z "$url" || "$url" =~ ^[[:space:]]*# ]] && continue
    
    if [[ "$filename" == *.zip ]]; then
        zip_file="$DATA_DIR/$filename"
        marker_file="$DATA_DIR/.${filename%.zip}_extracted"
        
        if [ -f "$marker_file" ]; then
            total_files=$((total_files + 1))
            # Check if extracted directory exists and has content
            extracted_dir="${filename%.zip}"
            if [ -d "$DATA_DIR/$extracted_dir" ] && [ "$(ls -A "$DATA_DIR/$extracted_dir" 2>/dev/null)" ]; then
                found_files=$((found_files + 1))
                file_count=$(find "$DATA_DIR/$extracted_dir" -type f | wc -l | tr -d ' ')
                log_message "Extracted {.path $extracted_dir} contains {.val ${file_count}} files"
            else
                log_message "Extracted directory {.path $extracted_dir} is empty or missing" --message-type warning
            fi
        fi
    fi
done <<< "$DOWNLOAD_LIST"

if [ $total_files -gt 0 ]; then
    log_message "Data extraction check: {.val ${found_files}}/{.val ${total_files}} archives extracted successfully"
fi


cleanup_temp_files "$DATA_DIR"

log_message "BrainSTEM data download and verification completed!" --message-type success
