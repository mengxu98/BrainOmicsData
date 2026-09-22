#!/bin/bash

# Download all public processed human 10x Multiome assets for Wang et al. 2025.
# Source:
# https://cellxgene.cziscience.com/collections/ad2149fc-19c5-41de-8cfe-44710fbada73
# Raw sequencing data are recorded at https://assets.nemoarchive.org/dat-oiif74w.

set -e

source "$(dirname "$0")/../functions/utils.sh"

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_ROOT="$(cd "$REPO_DIR/../.." && pwd)/data/BrainOmicsData"
DATA_DIR="$DATA_ROOT/raw/Wang_2025"

DOWNLOAD_LIST="
https://datasets.cellxgene.cziscience.com/a4310202-4dc8-4e1b-a96d-d9675f5b14d1.h5ad|Wang_2025_RNA.h5ad|2782116565|6d7b039bf95853abf022b2672d118dd39f111825ac35228acf0e1c2bdd36f970
https://datasets.cellxgene.cziscience.com/c9ba15d8-45a1-4f68-ac91-a59ad2d44a2d-fragment.tsv.bgz|Wang_2025_ATAC_fragments.tsv.bgz|42321523373|5b459e69f822b91c14a5800bff7fb6a2f8986ee488d55f5e7094de519c8de37f
https://datasets.cellxgene.cziscience.com/c9ba15d8-45a1-4f68-ac91-a59ad2d44a2d-fragment.tsv.bgz.tbi|Wang_2025_ATAC_fragments.tsv.bgz.tbi|6057453|6693a36b37d6c787fcc37c235633796b0962ae4d1fd364dcf542e3fbf056fdee
"

if [ "$#" -gt 0 ]; then
    download_query "$DOWNLOAD_LIST" "$@"
    exit $?
fi

log_message "Starting Wang et al. 2025 data download..."

batch_download "$DOWNLOAD_LIST" "$DATA_DIR" 5

cleanup_temp_files "$DATA_DIR"

log_message "Wang et al. 2025 data download completed!" --message-type success
