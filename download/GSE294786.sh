#!/bin/bash

# Download the complete processed GSE294786 source files.
# GEO packages Homo sapiens, Pan troglodytes, and Macaca mulatta together in
# the same counts and metadata files, so species filtering must happen after
# download.

set -e

source "$(dirname "$0")/../functions/utils.sh"

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_ROOT="${BRAINOMICS_DATA_ROOT:-$(cd "$REPO_DIR/../.." && pwd)/data/BrainOmicsData}"
DATA_DIR="$DATA_ROOT/raw/GSE294786"

DOWNLOAD_LIST="
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE294nnn/GSE294786/suppl/GSE294786_all_counts_transposed.tsv.gz|GSE294786_all_counts_transposed.tsv.gz|143357306|3838612a176223f87a800424d4e8d21830b7b92594bd6d132c2d7a2a4074ffd9
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE294nnn/GSE294786/suppl/GSE294786_all_meta_progenitor.tsv.gz|GSE294786_all_meta_progenitor.tsv.gz|2549099|533aff101a5c03acb42ea31fafe32ebb3df64beb8601b17f9b75fd65dccba907
"

if [ "$#" -gt 0 ]; then
    download_query "$DOWNLOAD_LIST" "$@"
    exit $?
fi

log_message "Starting GSE294786 data download..."

batch_download "$DOWNLOAD_LIST" "$DATA_DIR" 5

cleanup_temp_files "$DATA_DIR"

log_message "GSE294786 data download completed!" --message-type success
