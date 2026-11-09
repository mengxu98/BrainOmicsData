#!/bin/bash

set -e

source "$(dirname "$0")/../functions/utils.sh"

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_ROOT="${BRAINOMICS_DATA_ROOT:-$(cd "$REPO_DIR/../.." && pwd)/data/BrainOmicsData}"
DATA_DIR="$DATA_ROOT/raw/GSE207334"
OBJECT="$DATA_DIR/Ma_Sestan_mat.rds"
ARCHIVE="$OBJECT.gz"
SOURCE_URL="https://brainscope.gersteinlab.org/data/annotation_files/Ma_Sestan_mat.rds.gz"

if [ "$#" -gt 0 ]; then
    download_query "$SOURCE_URL|Ma_Sestan_mat.rds.gz|2169272289" "$@"
    exit $?
fi

if [ -f "$OBJECT" ]; then
    exit 0
fi

if [ ! -f "$ARCHIVE" ] || [ -f "$ARCHIVE.resume" ]; then
    download_with_resume "$SOURCE_URL" "$ARCHIVE" 5
fi

partial="$(mktemp "$OBJECT.part.XXXXXX")"
trap 'rm -f "$partial"' EXIT
gzip -dc "$ARCHIVE" > "$partial"
mv -n "$partial" "$OBJECT"
