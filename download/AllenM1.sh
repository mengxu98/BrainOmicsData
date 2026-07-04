#!/bin/bash

# Download script for AllenM1 public Allen Institute human M1 10x files.
# Source page:
# https://brain-map.org/our-research/cell-types-taxonomies/cell-types-database-rna-seq-data/human-m1-10x
# Current preprocessing consumes the BTSatlas-derived AllenM1_raw.rds split;
# these downloads record the original public Allen source files.

set -e

source "$(dirname "$0")/../functions/utils.sh"

DATA_DIR="../../data/BrainOmicsData/raw/AllenM1"

log_message "Starting AllenM1 source data download..."

DOWNLOAD_LIST="
https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/matrix.csv|matrix.csv|7726025779
https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/metadata.csv|metadata.csv|23963087
https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/tsne.csv|tsne.csv|4728663
https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/dend.json|dend.json|417271
https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/medians.csv|medians.csv|26579242
https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/trimmed_means.csv|trimmed_means.csv|31180368
https://allen-brain-map-cms-802451596237-us-west-2.s3.amazonaws.com/legacy%20files/human_dendrogram.rds|human_dendrogram.rds|0
https://allen-brain-map-cms-802451596237-us-west-2.s3.amazonaws.com/legacy%20files/sample-exp_component_mapping_human_10x_apr2020.zip|sample-exp_component_mapping_human_10x_apr2020.zip|0
https://cdn.prod.website-files.com/689cfbd308fa7373b604d290/68f12d80c5338c3dc5f012b7_taxonomy.txt|taxonomy.txt|0
"

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
