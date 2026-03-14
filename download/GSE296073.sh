#!/bin/bash

# Download script for GSE296073 data
# journal: Nature
# date: 2025
# paper: https://doi.org/10.1038/s41586-025-09362-8
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/40770097
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE296073
# GSE274829

set -e

source "$(dirname "$0")/../functions/utils.sh"

DATA_DIR="../../data/BrainOmicsData/raw/GSE296073"

log_message "Starting GSE296073 data download..."

DOWNLOAD_LIST="
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE296073&format=file|GSE296073_RAW.tar|0
https://zenodo.org/records/15299853/files/h_pre_peri_DY.rds?download=1|h_pre_peri_DY.rds|0
https://zenodo.org/records/15299853/files/INS.rds?download=1|INS.rds|0
"

# https://zenodo.org/records/15299853/files/organoid6w_DY.rds?download=1
# https://zenodo.org/records/15299853/files/img.rds?download=1

# Perform batch download
batch_download "$DOWNLOAD_LIST" "$DATA_DIR" 5


if [ -z "$(find "$DATA_DIR" -type f -name "*.gz" 2>/dev/null | head -1)" ]; then
	log_message "Extracting $DATA_DIR/GSE296073_RAW.tar ..."
	tar -xvf "$DATA_DIR/GSE296073_RAW.tar" -C "$DATA_DIR"
	
	log_message "Renaming files: replacing '-' with '_'..."
	find "$DATA_DIR" -maxdepth 1 -type f -name "*.gz" | while read -r file; do
		dir=$(dirname "$file")
		base=$(basename "$file")
		if [[ "$base" == *"-"* ]]; then
			new_base=$(echo "$base" | sed 's/-/_/g')
			mv "$file" "$dir/$new_base"
		fi
	done
	log_message "Archive extracted and files renamed" --message-type success
else
	log_message "Archive already extracted (extracted files exist), skipping tar extraction."
fi

ORIGINAL_DIR=$(pwd)
cd "$DATA_DIR"
GZ_FILES_IN_ROOT=$(find . -maxdepth 1 -type f -name "*.gz" 2>/dev/null | wc -l | tr -d ' ')
GZ_FILES_IN_SUBDIRS=$(find . -mindepth 2 -type f -name "*.gz" 2>/dev/null | wc -l | tr -d ' ')

if [ "$GZ_FILES_IN_ROOT" -gt 0 ] && [ "$GZ_FILES_IN_SUBDIRS" -eq 0 ]; then
	log_message "Organizing files by sample ID..."

	find . -maxdepth 1 -type f -name "*.gz" | while read -r file; do
		base=$(basename "$file")
		# Extract folder name: keep full prefix including GSM (e.g., GSM8964678_2014028cx_all)
		# Pattern: GSM8964678_2014028cx_all_barcodes.tsv.gz -> GSM8964678_2014028cx_all
		folder_name=$(echo "$base" | sed -E 's/_(barcodes|features|matrix)\.(tsv|mtx)\.gz$//')
		# Extract 10X standard filename (e.g., barcodes.tsv.gz)
		# Pattern: GSM8964678_2014028cx_all_barcodes.tsv.gz -> barcodes.tsv.gz
		standard_name=$(echo "$base" | sed -E 's/^.*_(barcodes|features|matrix)\.(tsv|mtx)\.gz$/\1.\2.gz/')
		
		if [ -n "$folder_name" ] && [ -n "$standard_name" ]; then
			mkdir -p "$folder_name"
			mv "$file" "$folder_name/$standard_name" 2>/dev/null || true
		fi
	done
	
	log_message "Files organized by sample ID" --message-type success
elif [ "$GZ_FILES_IN_SUBDIRS" -gt 0 ]; then
	log_message "Files already organized (found .gz files in subdirectories), skipping organization step."
else
	log_message "No .gz files found to organize."
fi

cd "$ORIGINAL_DIR"

cleanup_temp_files "$DATA_DIR"

log_message "GSE296073 data download and organization completed!" --message-type success
