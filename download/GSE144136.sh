#!/bin/bash

# Download script for GSE144136 public GEO supplementary files.
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144136
# ftp: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144136/suppl/

set -e

source "$(dirname "$0")/../functions/utils.sh"

DATA_DIR="../../data/BrainOmicsData/raw/GSE144136"

log_message "Starting GSE144136 data download..."

DOWNLOAD_LIST="
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144136/suppl/GSE144136_Astros_3_filtered_cells.csv.gz|GSE144136_Astros_3_filtered_cells.csv.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144136/suppl/GSE144136_CellNames.csv.gz|GSE144136_CellNames.csv.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144136/suppl/GSE144136_Ex_10_L2_4_filtered_cells.csv.gz|GSE144136_Ex_10_L2_4_filtered_cells.csv.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144136/suppl/GSE144136_Ex_1_L5_6_filtered_cells.csv.gz|GSE144136_Ex_1_L5_6_filtered_cells.csv.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144136/suppl/GSE144136_Ex_3_L4_5_filtered_cells.csv.gz|GSE144136_Ex_3_L4_5_filtered_cells.csv.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144136/suppl/GSE144136_Ex_4_L_6_filtered_cells.csv.gz|GSE144136_Ex_4_L_6_filtered_cells.csv.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144136/suppl/GSE144136_Ex_5_L5_filtered_cells.csv.gz|GSE144136_Ex_5_L5_filtered_cells.csv.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144136/suppl/GSE144136_Ex_6_L4_6_filtered_cells.csv.gz|GSE144136_Ex_6_L4_6_filtered_cells.csv.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144136/suppl/GSE144136_Ex_7_L4_6_filtered_cells.csv.gz|GSE144136_Ex_7_L4_6_filtered_cells.csv.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144136/suppl/GSE144136_Ex_8_L5_6_filtered_cells.csv.gz|GSE144136_Ex_8_L5_6_filtered_cells.csv.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144136/suppl/GSE144136_Ex_9_L5_6_filtered_cells.csv.gz|GSE144136_Ex_9_L5_6_filtered_cells.csv.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144136/suppl/GSE144136_GRCh38-1.2.0_premrna.tar.gz|GSE144136_GRCh38-1.2.0_premrna.tar.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144136/suppl/GSE144136_GeneBarcodeMatrix_Annotated.mtx.gz|GSE144136_GeneBarcodeMatrix_Annotated.mtx.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144136/suppl/GSE144136_GeneNames.csv.gz|GSE144136_GeneNames.csv.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144136/suppl/GSE144136_Inhib_1_filtered_cells.csv.gz|GSE144136_Inhib_1_filtered_cells.csv.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144136/suppl/GSE144136_Inhib_2_VIP_filtered_cells.csv.gz|GSE144136_Inhib_2_VIP_filtered_cells.csv.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144136/suppl/GSE144136_Inhib_3_SST_filtered_cells.csv.gz|GSE144136_Inhib_3_SST_filtered_cells.csv.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144136/suppl/GSE144136_Inhib_5_filtered_cells.csv.gz|GSE144136_Inhib_5_filtered_cells.csv.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144136/suppl/GSE144136_Inhib_6_SST_filtered_cells.csv.gz|GSE144136_Inhib_6_SST_filtered_cells.csv.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144136/suppl/GSE144136_Inhib_7_PVALB_filtered_cells.csv.gz|GSE144136_Inhib_7_PVALB_filtered_cells.csv.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144136/suppl/GSE144136_Inhib_8_PVALB_filtered_cells.csv.gz|GSE144136_Inhib_8_PVALB_filtered_cells.csv.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144136/suppl/GSE144136_Micro_Macro_filtered_cells.csv.gz|GSE144136_Micro_Macro_filtered_cells.csv.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144136/suppl/GSE144136_OPCs_1_filtered_cells.csv.gz|GSE144136_OPCs_1_filtered_cells.csv.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144136/suppl/GSE144136_OPCs_2_filtered_cells.csv.gz|GSE144136_OPCs_2_filtered_cells.csv.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144136/suppl/GSE144136_Oligos_1_filtered_cells.csv.gz|GSE144136_Oligos_1_filtered_cells.csv.gz|0
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144136/suppl/GSE144136_Oligos_3_filtered_cells.csv.gz|GSE144136_Oligos_3_filtered_cells.csv.gz|0
"

batch_download "$DOWNLOAD_LIST" "$DATA_DIR" 5

for archive in "$DATA_DIR"/*.tar "$DATA_DIR"/*.tar.gz "$DATA_DIR"/*.tgz; do
    [ -f "$archive" ] || continue
    marker="$DATA_DIR/.${archive##*/}_extracted"
    if [ ! -f "$marker" ]; then
        log_message "Extracting {.file ${archive##*/}}..."
        tar -xf "$archive" -C "$DATA_DIR"
        touch "$marker"
    fi
done

cleanup_temp_files "$DATA_DIR"

log_message "GSE144136 data download completed!" --message-type success
