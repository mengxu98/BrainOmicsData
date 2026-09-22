#!/usr/bin/env bash

# Download the complete public author-processed matrices and per-nucleus
# metadata for Cameron et al.  The EGA FASTQ files are controlled, but the five
# brain-region raw count matrices and author cell-type annotations are public on
# Figshare and are sufficient for the formal processed-data reconstruction.
#
# Paper: https://doi.org/10.1016/j.biopsych.2022.06.033
# Figshare: https://figshare.com/articles/dataset/11629311
# EGA raw data: https://ega-archive.org/studies/EGAS00001006537

set -e

source "$(dirname "$0")/../functions/utils.sh"

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_ROOT="$(cd "$REPO_DIR/../.." && pwd)/data/BrainOmicsData"
DATA_DIR="$DATA_ROOT/raw/EGAS00001006537"

DOWNLOAD_LIST="
https://ndownloader.figshare.com/files/36712263|cameron_2022_snRNAseq_Cer_metadata.txt.gz|146690|55edf8c2be22a41d13d4d21bfc804f5d9df79e3469b1bca13b953cdde4b5c10c
https://ndownloader.figshare.com/files/36712266|cameron_2022_snRNAseq_Cer_raw_count_gEX_matrix.txt.gz|68748175|319d7996f851955d06531f358ccf1bdf7a0fd4c3b18f24fbabda6308a82eea30
https://ndownloader.figshare.com/files/36712269|cameron_2022_snRNAseq_FC_metadata.txt.gz|66717|02db94cb17b95b221fde335bb83fc84dc4b68653639d07592b11b48a00b19b97
https://ndownloader.figshare.com/files/36712272|cameron_2022_snRNAseq_FC_raw_count_gEX_matrix.txt.gz|30885000|a256f55033553df8205a62a72fbf53664c7d3cb25bd45941f3f35d98f46873b2
https://ndownloader.figshare.com/files/36712275|cameron_2022_snRNAseq_GE_metadata.txt.gz|53680|e27e92598522dac567131a65ebdf0fecdd8fc8cd251cbb8c67c771422165e66b
https://ndownloader.figshare.com/files/36712278|cameron_2022_snRNAseq_GE_raw_count_gEX_matrix.txt.gz|20347328|a1edc5cccf7d960574b0dd13b61d2e544b2d127a6dd22c84b9dd609e7bb127cb
https://ndownloader.figshare.com/files/36712281|cameron_2022_snRNAseq_Hipp_metadata.txt.gz|77101|94fd1a45508070ca30e8a731caed63c0b5162d90b659e97206af17bcd54011e7
https://ndownloader.figshare.com/files/36712284|cameron_2022_snRNAseq_Hipp_raw_count_gEX_matrix.txt.gz|22902156|bdabbd4ecbd79499105965755ecc395142a79a1d293d84c068afe752b60cce21
https://ndownloader.figshare.com/files/36712287|cameron_2022_snRNAseq_Thal_metadata.txt.gz|114448|49fc799834819718a9c334c400b14b05bb81bbd844cc7c6b6db0689ea13edb89
https://ndownloader.figshare.com/files/36712290|cameron_2022_snRNAseq_Thal_raw_count_gEX_matrix.txt.gz|45408602|b16bca510e53d5fc220390c49758fcb8766bcbf58aaead26ba59e085112906af
"

if [ "$#" -gt 0 ]; then
    download_query "$DOWNLOAD_LIST" "$@"
    exit $?
fi

log_message "Starting EGAS00001006537 public Figshare source download..."

batch_download "$DOWNLOAD_LIST" "$DATA_DIR" 5

cleanup_temp_files "$DATA_DIR"

log_message "EGAS00001006537 public Figshare source download completed!" --message-type success
