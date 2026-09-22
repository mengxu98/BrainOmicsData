#!/bin/bash

# Download script for AllenM1 public Allen Institute human M1 10x files.
# Source page:
# https://brain-map.org/our-research/cell-types-taxonomies/cell-types-database-rna-seq-data/human-m1-10x
# The processing workflow consumes these original Allen files directly.  The
# Formal input is the complete author-released Allen human M1 matrix.

set -e

source "$(dirname "$0")/../functions/utils.sh"

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_ROOT="$(cd "$REPO_DIR/../.." && pwd)/data/BrainOmicsData"
DATA_DIR="$DATA_ROOT/raw/AllenM1"

DOWNLOAD_LIST="
https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/matrix.csv|matrix.csv|7726025779|5a300c987fc7a29feb671a587a14f203efe3a4d9cbac4da5fbffc7e66e4cc263
https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/metadata.csv|metadata.csv|23963087|8fc132a419683e812d099e3235624c7b8c04a4658edbd6d636087655624af2de
https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/tsne.csv|tsne.csv|4728663|b6e8c05d634431cddee7b541b7364e3e7551f64d1edb76f9fb4a93366790539e
https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/dend.json|dend.json|417271|449df48e7cf89b8326f48776bcade313042c804143fff58aeba35d0ca8290ee4
https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/medians.csv|medians.csv|26579242|c5d4e7d12c4a2b971652ae76ee5a4e458d64466fc0518d7774b3f42383049480
https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/trimmed_means.csv|trimmed_means.csv|31180368|0d87f976a722f4874f57241a39cb70a3e5b18e23603d382bba7068ecfc51545d
https://allen-brain-map-cms-802451596237-us-west-2.s3.amazonaws.com/legacy%20files/human_dendrogram.rds|human_dendrogram.rds|7072|ce9ea2a1e058fa64e8f486d231712762594615fe70af9695da9472ee1cd0d351
https://allen-brain-map-cms-802451596237-us-west-2.s3.amazonaws.com/legacy%20files/sample-exp_component_mapping_human_10x_apr2020.zip|sample-exp_component_mapping_human_10x_apr2020.zip|604388|88a4f3776fccf8dae4ea1a7e3073a45496bfe3c7a4ddd8f7f80f4129790d1197
https://cdn.prod.website-files.com/689cfbd308fa7373b604d290/68f12d80c5338c3dc5f012b7_taxonomy.txt|taxonomy.txt|410|3989bc517e96cb572218dd97d98af4a6b782c18dbf98f1ddb0e32048d7795a5c
"

if [ "$#" -gt 0 ]; then
    download_query "$DOWNLOAD_LIST" "$@"
    exit $?
fi

log_message "Starting AllenM1 source data download..."

batch_download "$DOWNLOAD_LIST" "$DATA_DIR" 5

for archive in "$DATA_DIR"/*.zip; do
    [ -f "$archive" ] || continue
    marker="$DATA_DIR/.${archive##*/}_extracted"
    if [ ! -f "$marker" ]; then
        log_message "Extracting {.file ${archive##*/}}..."
        unzip -q -o "$archive" -d "$DATA_DIR"
        touch "$marker"
    fi
done

cleanup_temp_files "$DATA_DIR"

log_message "AllenM1 source data download completed!" --message-type success
