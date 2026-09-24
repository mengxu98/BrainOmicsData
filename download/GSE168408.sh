#!/usr/bin/env bash

# title: Human prefrontal cortex gene regulatory dynamics from gestation to
# adulthood at single-cell resolution
# paper: https://doi.org/10.1016/j.cell.2022.09.039
# project: https://brain.listerlab.org/
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168408
# code: https://github.com/ListerLab/pfc_development
#
# The GEO archive is downloaded in full and contains all released human
# snRNA-seq, snATAC-seq and organoid source files. RNA preprocessing
# uses the author-released RNA-all H5AD, whose default X matrix is full raw
# counts, together with its independent barcode and gene metadata files.

set -e

source "$(dirname "$0")/../functions/utils.sh"

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_ROOT="${BRAINOMICS_DATA_ROOT:-$(cd "$REPO_DIR/../.." && pwd)/data/BrainOmicsData}"
DATA_DIR="$DATA_ROOT/raw/GSE168408"

DOWNLOAD_LIST="
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE168nnn/GSE168408/suppl/GSE168408_RAW.tar|GSE168408_RAW.tar|23884810240|
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE168nnn/GSE168408/suppl/filelist.txt|GSE168408_filelist.txt|10818|54200f3e0f3dea94a31ac7d39c29580f40ade5e4834eb0956db75fad908516d2
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE168nnn/GSE168408/matrix/GSE168408-GPL21697_series_matrix.txt.gz|GSE168408-GPL21697_series_matrix.txt.gz|2945|3ceb549785023fad08c407287cbf3ba64ad60310187bd541ceeeb934f555b032
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE168nnn/GSE168408/matrix/GSE168408-GPL24676_series_matrix.txt.gz|GSE168408-GPL24676_series_matrix.txt.gz|5696|f4065503cc195b56ae96a2cb434de2c8a93223f95229dee52906809b0eb8f674
https://storage.googleapis.com/neuro-dev/Processed_data/RNA-all_full-counts-and-downsampled-CPM.h5ad|RNA-all_full-counts-and-downsampled-CPM.h5ad|4963815084|5ca0b1ebe638f13e4302b721b80f41755a759c87f7fe43da6aeb7e8f9d0fd9ce
https://storage.googleapis.com/neuro-dev/Processed_data/RNA-all_BCs-meta-data.csv|RNA-all_BCs-meta-data.csv|99546082|2265c228291af73099aaa3092f5bea6ae4c713bc057676b7a0f3ed9e3c26a30a
https://storage.googleapis.com/neuro-dev/Processed_data/RNA-all_genes-meta-data.csv|RNA-all_genes-meta-data.csv|4754987|1b0652e10800cb4d867da39b12abe893da8855958cc6613617dfc02cbac683fb
"

if [ "$#" -gt 0 ]; then
    download_query "$DOWNLOAD_LIST" "$@"
    exit $?
fi

log_message "Starting complete GSE168408 data download..."

batch_download "$DOWNLOAD_LIST" "$DATA_DIR" 8

archive_file="$DATA_DIR/GSE168408_RAW.tar"
extract_marker="$DATA_DIR/.GSE168408_RAW.tar_extracted"
if [ ! -f "$extract_marker" ]; then
    log_message "Extracting complete GSE168408 GEO archive..."
    tar -xf "$archive_file" -C "$DATA_DIR"
    touch "$extract_marker"
fi

rna_matrix_count=$(
    find "$DATA_DIR" -type f \
        -name '*filtered*.h5*' \
        | wc -l | tr -d ' '
)
atac_fragment_count=$(
    find "$DATA_DIR" -type f -name '*_snATAC_fragments.tsv.gz' \
        | wc -l | tr -d ' '
)
if [ "$rna_matrix_count" -lt 32 ] || [ "$atac_fragment_count" -lt 17 ]; then
    log_message \
        "GSE168408 GEO archive lacks the complete RNA or ATAC inventory" \
        --message-type error || true
    exit 1
fi

cleanup_temp_files "$DATA_DIR"

log_message \
    "GSE168408 complete source download finished" \
    --message-type success
