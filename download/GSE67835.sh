#!/bin/bash

# Download script for GSE67835 data
# paper: https://doi.org/10.1073/pnas.1507125112
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/26060301
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67835

set -e

source "$(dirname "$0")/../functions/utils.sh"

DATA_DIR="../../data/BrainOmicsData/raw/GSE67835"

log_message "Starting GSE67835 data download..."

DOWNLOAD_LIST="
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE67835&format=file|GSE67835_RAW.tar|0
"

# Perform batch download
batch_download "$DOWNLOAD_LIST" "$DATA_DIR" 5


# Extract the tar file (idempotent)
if [ ! -f "$DATA_DIR/.GSE67835_extracted" ]; then
	log_message "Extracting $DATA_DIR/GSE67835_RAW.tar ..."
	# Create GSE67835 subdirectory
	mkdir -p "$DATA_DIR/GSE67835"
	tar -xvf "$DATA_DIR/GSE67835_RAW.tar" -C "$DATA_DIR/GSE67835"
	# create a marker so we don't re-extract next runs
	touch "$DATA_DIR/.GSE67835_extracted"
else
	log_message "Archive already extracted (marker exists), skipping tar extraction."
fi

# Extract all .gz files (idempotent)
if [ ! -f "$DATA_DIR/.GSE67835_gz_extracted" ]; then
	log_message "Extracting .gz files from $DATA_DIR/GSE67835/ ..."
	find "$DATA_DIR/GSE67835" -name "*.csv.gz" -type f | while read gz_file; do
		# Extract to same directory, remove .gz extension
		gunzip -f "$gz_file"
	done
	# create a marker so we don't re-extract next runs
	touch "$DATA_DIR/.GSE67835_gz_extracted"
else
	log_message "GZ files already extracted (marker exists), skipping gz extraction."
fi


cleanup_temp_files "$DATA_DIR"

log_message "GSE67835 data download and organization completed!" --message-type success
