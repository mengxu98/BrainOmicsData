#!/bin/bash

# Download script for public HYPOMAP CELLxGENE data.
# Collection: https://cellxgene.cziscience.com/collections/d0941303-7ce3-4422-9249-cf31eb98c480
# h5ad: https://datasets.cellxgene.cziscience.com/af06d70d-a2bc-43c1-8378-bc545bc8ee2c.h5ad

set -e

source "$(dirname "$0")/../functions/utils.sh"

DATA_DIR="../../data/BrainOmicsData/raw/HYPOMAP"

log_message "Starting HYPOMAP data download..."

DOWNLOAD_LIST="
https://datasets.cellxgene.cziscience.com/af06d70d-a2bc-43c1-8378-bc545bc8ee2c.h5ad|human_HYPOMAP_snRNASeq.h5ad|5231283640
"

batch_download "$DOWNLOAD_LIST" "$DATA_DIR" 5

cleanup_temp_files "$DATA_DIR"

log_message "HYPOMAP data download completed!" --message-type success
