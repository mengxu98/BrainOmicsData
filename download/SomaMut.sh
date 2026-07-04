#!/bin/bash

# Download script for public SomaMut processed files.
# Main page: https://publications.wenglab.org/SomaMut/
# Human Brain Aging: https://publications.wenglab.org/SomaMut/Jeffries_Yu_BrainAging_2025/
# Parkinson's Disease: https://publications.wenglab.org/SomaMut/Ziegenfuss_Yu_PD_2025/
# Dog Brain Aging: https://publications.wenglab.org/SomaMut/Class_Yu_Dog_2026/
# Current processing/SomaMut.R uses Jeffries_Yu_BrainAging_2025/pfc.clean.rds.

set -e

source "$(dirname "$0")/../functions/utils.sh"

DATA_DIR="../../data/BrainOmicsData/raw/SomaMut"

log_message "Starting SomaMut processed data download..."

DOWNLOAD_LIST="
https://publications.wenglab.org/SomaMut/data/Jeffries_Yu_BrainAging_2025/pfc.clean.rds|pfc.clean.rds|28330539224
https://publications.wenglab.org/SomaMut/data/Jeffries_Yu_BrainAging_2025/pfc.clean.filtered_metadata.txt|pfc.clean.filtered_metadata.txt|80149123
https://publications.wenglab.org/SomaMut/data/Jeffries_Yu_BrainAging_2025/sample_information.xlsx|sample_information.xlsx|45536
https://publications.wenglab.org/SomaMut/data/Jeffries_Yu_BrainAging_2025/pfc.clean.cpm.tar.gz|pfc.clean.cpm.tar.gz|62586658
https://publications.wenglab.org/SomaMut/data/Jeffries_Yu_BrainAging_2025/pfc.clean.log_cpm.tar.gz|pfc.clean.log_cpm.tar.gz|62657513
https://publications.wenglab.org/SomaMut/data/Jeffries_Yu_BrainAging_2025/pfc.clean.pct.tar.gz|pfc.clean.pct.tar.gz|35157858
https://publications.wenglab.org/SomaMut/data/Jeffries_Yu_BrainAging_2025/pfc.clean.cv.tar.gz|pfc.clean.cv.tar.gz|56746497
https://publications.wenglab.org/SomaMut/data/Jeffries_Yu_BrainAging_2025/MERFISH_gene_panel.txt|MERFISH_gene_panel.txt|7855
https://publications.wenglab.org/SomaMut/data/Jeffries_Yu_BrainAging_2025/all.multianno.moreInfo.txt|all.multianno.moreInfo.txt|16772058
https://publications.wenglab.org/SomaMut/data/Jeffries_Yu_BrainAging_2025/signature_mutationPattern.txt|signature_mutationPattern.txt|4077
https://publications.wenglab.org/SomaMut/data/Class_Yu_DogNeuronAging_2026/dog.somatic_mutation.tar.gz|Class_Yu_DogNeuronAging_2026/dog.somatic_mutation.tar.gz|71886
https://publications.wenglab.org/SomaMut/data/Class_Yu_DogNeuronAging_2026/human.somatic_mutation.tar.gz|Class_Yu_DogNeuronAging_2026/human.somatic_mutation.tar.gz|703219
https://publications.wenglab.org/SomaMut/data/Class_Yu_DogNeuronAging_2026/signatures.indel.rds|Class_Yu_DogNeuronAging_2026/signatures.indel.rds|131809
https://publications.wenglab.org/SomaMut/data/Class_Yu_DogNeuronAging_2026/signatures.snv.rds|Class_Yu_DogNeuronAging_2026/signatures.snv.rds|163400
https://publications.wenglab.org/SomaMut/data/Class_Yu_DogNeuronAging_2026/summary.dog+human.txt|Class_Yu_DogNeuronAging_2026/summary.dog+human.txt|41922
"

batch_download "$DOWNLOAD_LIST" "$DATA_DIR" 5

for archive in "$DATA_DIR"/*.tar.gz "$DATA_DIR"/Class_Yu_DogNeuronAging_2026/*.tar.gz; do
    [ -f "$archive" ] || continue
    marker="$DATA_DIR/.${archive#"$DATA_DIR"/}_extracted"
    marker="${marker//\//_}"
    if [ ! -f "$marker" ]; then
        log_message "Extracting {.file ${archive#"$DATA_DIR"/}}..."
        tar -xf "$archive" -C "$(dirname "$archive")"
        touch "$marker"
    fi
done

cleanup_temp_files "$DATA_DIR"

log_message "SomaMut processed data download completed!" --message-type success
