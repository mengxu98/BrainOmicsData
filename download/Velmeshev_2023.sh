#!/bin/bash

# Download the complete public CELLxGENE object for Velmeshev et al. 2023.
# Do not subset Dataset or donor IDs at download time. Preprocessing retains the
# complete object, then creates the analysis object from dataset == "Velmeshev"
# while excluding donor 5936 as a duplicate of GSE204683.
#
# The paper's NIHMS1972877-supplement-Data_S1.xlsx is not placed in
# DOWNLOAD_LIST because PMC returns a browser challenge to unattended curl.
# Supporting metadata record (not a processing input):
# URL: https://pmc.ncbi.nlm.nih.gov/articles/instance/11005279/bin/NIHMS1972877-supplement-Data_S1.xlsx
# bytes: 19895
# SHA-256: 72edad3894eac42fe64da8f746eb181549d144341d0ed7a7ead06953a1139c39
# The complete, directly downloadable UCSC cell metadata below is the reviewed
# processing input.

set -e

source "$(dirname "$0")/../functions/utils.sh"

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_ROOT="$(cd "$REPO_DIR/../.." && pwd)/data/BrainOmicsData"
DATA_DIR="$DATA_ROOT/raw/Velmeshev_2023"

DOWNLOAD_LIST="
https://datasets.cellxgene.cziscience.com/2001fbc7-de84-4e84-9a76-4290bcd65415.h5ad|Velmeshev_2023_full_cellxgene.h5ad|5976097241|31709ba1b4273f40da49571f86aca50012a2424ad85292bf1fdbe952398c2221
https://cells.ucsc.edu/pre-postnatal-cortex/all/rna/meta.tsv|Velmeshev_2023_ucsc_metadata.tsv|83775768|dae25e8fd27c327707eae88a7d4c4aae54d704cd3a8ee13d07c356e89978c18b
"

if [ "$#" -gt 0 ]; then
    download_query "$DOWNLOAD_LIST" "$@"
    exit $?
fi

log_message "Starting Velmeshev et al. 2023 full-object download..."

batch_download "$DOWNLOAD_LIST" "$DATA_DIR" 5

cleanup_temp_files "$DATA_DIR"

log_message "Velmeshev et al. 2023 full-object download completed!" --message-type success
