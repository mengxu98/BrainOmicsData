#!/bin/bash

# Download script for the public Brain Transcriptome Single-cell Atlas files.
# Source: https://zenodo.org/records/10939707
# Current processing/BTSatlas_1.py requires BTS_atlas.h5ad and splits several
# integration datasets from it, including AllenM1, EGAD00001006049,
# EGAS00001006537, GSE144136, GSE168408, GSE178175, and GSE202210.

set -e

source "$(dirname "$0")/../functions/utils.sh"

DATA_DIR="../../data/BrainOmicsData/raw/BTSatlas"

log_message "Starting BTSatlas data download..."

DOWNLOAD_LIST="
https://zenodo.org/api/records/10939707/files/BTS_atlas.h5ad/content|BTS_atlas.h5ad|10254118224
https://zenodo.org/api/records/10939707/files/BTS_atlas_CellTypist_model.pkl/content|BTS_atlas_CellTypist_model.pkl|328657
"

# RiskGene_multipanel.zip is public but not required by the current processing
# workflow and is therefore omitted by default.

batch_download "$DOWNLOAD_LIST" "$DATA_DIR" 5

cleanup_temp_files "$DATA_DIR"

log_message "BTSatlas data download completed!" --message-type success
