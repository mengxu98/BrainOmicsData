#!/bin/bash

# Enhanced download script for ROSMAP UCSC snATAC-seq Epigenomic data
# Usage: ./download_rosmap_ucsc_snATAC-seq_Epigenomic.sh

set -e

source "$(dirname "$0")/../functions/utils.sh"

# https://cells.ucsc.edu/?ds=rosmap-ad-aging-brain+ad-atac+erosion

DATA_DIR="../../data/BrainOmicsData/raw/ROSMAP/ATAC_Epigenomic"

log_message "Starting ROSMAP UCSC snATAC-seq Epigenomic data download..."

DOWNLOAD_LIST="
https://cells.ucsc.edu/rosmap-ad-aging-brain/ad-atac/erosion/matrix.mtx.gz|matrix.mtx.gz|0
https://cells.ucsc.edu/rosmap-ad-aging-brain/ad-atac/erosion/features.tsv.gz|features.tsv.gz|0
https://cells.ucsc.edu/rosmap-ad-aging-brain/ad-atac/erosion/barcodes.tsv.gz|barcodes.tsv.gz|0
https://cells.ucsc.edu/rosmap-ad-aging-brain/ad-atac/erosion/meta.tsv|meta.tsv|0
https://cells.ucsc.edu/rosmap-ad-aging-brain/ad-atac/erosion/UMAP_coordinates.coords.tsv.gz|UMAP_coordinates.coords.tsv.gz|0
"

# Perform batch download
batch_download "$DOWNLOAD_LIST" "$DATA_DIR" 5


cleanup_temp_files "$DATA_DIR"

log_message "ROSMAP UCSC snATAC-seq Epigenomic data download completed!" --message-type success
