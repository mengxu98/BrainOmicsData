#!/bin/bash

# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261983

set -e

source "$(dirname "$0")/../functions/utils.sh"

DATA_DIR="../../data/BrainOmicsData/raw/GSE261983"

log_message "Starting GSE261983 data download..."

DOWNLOAD_LIST="
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE261983&format=file|GSE261983_RAW.tar|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE261983&format=file&file=GSE261983%5FATAC%5Fpeakset%2Ebed%2Egz|GSE261983_ATAC_peakset.bed.gz|0
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE261983&format=file&file=GSE261983%5FRNA%5Fgene%5Flist%2Etxt%2Egz|GSE261983_RNA_gene_list.txt.gz|0
"

# Perform batch download
batch_download "$DOWNLOAD_LIST" "$DATA_DIR" 5


TAR_FILE="$DATA_DIR/GSE261983_RAW.tar"
# 检查是否有样本文件夹（解压后的标志）
SAMPLE_DIRS_EXIST=$(find "$DATA_DIR" -maxdepth 1 -type d -name "GSM*_*_RNA" 2>/dev/null | wc -l | tr -d ' ')

if [ -f "$TAR_FILE" ] && [ "$SAMPLE_DIRS_EXIST" -eq 0 ]; then
	log_message "Extracting $TAR_FILE ..."
	tar -xvf "$TAR_FILE" -C "$DATA_DIR"
	
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
elif [ "$SAMPLE_DIRS_EXIST" -gt 0 ]; then
	log_message "Archive already extracted (sample directories exist), skipping tar extraction."
else
	log_message "Warning: Tar file not found or extraction check failed."
fi

ORIGINAL_DIR=$(pwd)
cd "$DATA_DIR"

# 排除已知的辅助文件，只处理样本数据文件
EXCLUDE_PATTERNS="GSE261983_ATAC_peakset.bed.gz|GSE261983_RNA_gene_list.txt.gz|GSE261983_RAW.tar"

GZ_FILES_IN_ROOT=$(find . -maxdepth 1 -type f -name "*.gz" 2>/dev/null | grep -vE "$EXCLUDE_PATTERNS" | wc -l | tr -d ' ')
GZ_FILES_IN_SUBDIRS=$(find . -mindepth 2 -type f -name "*.gz" 2>/dev/null | wc -l | tr -d ' ')

if [ "$GZ_FILES_IN_ROOT" -gt 0 ] && [ "$GZ_FILES_IN_SUBDIRS" -eq 0 ]; then
	log_message "Organizing files by sample ID and renaming to 10X standard names..."
	
	find . -maxdepth 1 -type f -name "*.gz" | grep -vE "$EXCLUDE_PATTERNS" | while read -r file; do
		base=$(basename "$file")
		# Extract folder name: keep full prefix including GSM (e.g., GSM8155483_RT00372N_RNA)
		# Pattern: GSM8155483_RT00372N_RNA_barcodes.txt.gz -> GSM8155483_RT00372N_RNA
		folder_name=$(echo "$base" | sed -E 's/_(barcodes|features|matrix)\.(txt|tsv|mtx)\.gz$//')
		
		# Extract file type and determine 10X standard filename
		# Pattern: GSM8155483_RT00372N_RNA_barcodes.txt.gz -> barcodes.tsv.gz
		# Pattern: GSM8155483_RT00372N_RNA_matrix.mtx.gz -> matrix.mtx.gz
		if [[ "$base" =~ _(barcodes|features|matrix)\.(txt|tsv|mtx)\.gz$ ]]; then
			file_type="${BASH_REMATCH[1]}"
			file_ext="${BASH_REMATCH[2]}"
			
			# Convert to 10X standard: barcodes.txt.gz -> barcodes.tsv.gz, matrix.mtx.gz -> matrix.mtx.gz
			if [ "$file_type" == "barcodes" ] && [ "$file_ext" == "txt" ]; then
				standard_name="barcodes.tsv.gz"
			elif [ "$file_type" == "features" ] && [ "$file_ext" == "txt" ]; then
				standard_name="features.tsv.gz"
			elif [ "$file_type" == "matrix" ]; then
				standard_name="matrix.mtx.gz"
			else
				standard_name="${file_type}.${file_ext}.gz"
			fi
			
			if [ -n "$folder_name" ] && [ -n "$standard_name" ]; then
				mkdir -p "$folder_name"
				mv "$file" "$folder_name/$standard_name" 2>/dev/null || true
			fi
		fi
	done
	
	# Count processed directories
	PROCESSED_DIRS=$(find . -maxdepth 1 -type d -name "GSM*_*_RNA" 2>/dev/null | wc -l | tr -d ' ')
	log_message "Files organized by sample ID and renamed to 10X standard names (organized $PROCESSED_DIRS sample directories)" --message-type success
elif [ "$GZ_FILES_IN_SUBDIRS" -gt 0 ]; then
	log_message "Files already organized (found .gz files in subdirectories), skipping organization step."
else
	log_message "No sample data files found to organize."
fi

cd "$ORIGINAL_DIR"

cleanup_temp_files "$DATA_DIR"

log_message "GSE261983 data download and organization completed!" --message-type success
