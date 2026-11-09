#!/bin/bash

set -e

source "$(dirname "$0")/../functions/utils.sh"

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_ROOT="${BRAINOMICS_DATA_ROOT:-$(cd "$REPO_DIR/../.." && pwd)/data/BrainOmicsData}"
DATA_DIR="$DATA_ROOT/raw/Li_et_al_2018/rawData"

SOURCE_URL="http://development.psychencode.org/files/processed_data/scRNA-seq"
FETAL="$SOURCE_URL/Sestan.fetalHuman.Psychencode.Rdata|Sestan.fetalHuman.Psychencode.Rdata|28475135"
ADULT="$SOURCE_URL/Sestan.adultHumanNuclei.Psychencode.Rdata|Sestan.adultHumanNuclei.Psychencode.Rdata|182212983"

if [ "$#" -gt 0 ]; then
    download_query "$FETAL"$'\n'"$ADULT" "$@"
    exit $?
fi

batch_download "$FETAL" "$DATA_DIR/Prenatal_scRNA_seq" 5
batch_download "$ADULT" "$DATA_DIR/Adult_snRNA_seq" 5
