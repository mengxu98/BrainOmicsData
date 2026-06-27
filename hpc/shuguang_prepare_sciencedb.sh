#!/usr/bin/env bash
set -euo pipefail

unset LC_ALL LC_CTYPE LANG LANGUAGE
export LANG=C

REPO_DIR="${REPO_DIR:-/work/home/mengxu1310/repositories/BrainOmicsData}"
INTEGRATION_DIR="${INTEGRATION_DIR:-/work/home/mengxu1310/data/BrainOmicsData/integration}"
SCIENCEDB_OUT_DIR="${SCIENCEDB_OUT_DIR:-/work/home/mengxu1310/data/ScienceDB/HumanBrain_sc_snRNAseq_age_resource}"
OBJECT_FILE="${OBJECT_FILE:-$INTEGRATION_DIR/objects_celltypes.rds}"
LISI_RESULTS_FILE="${LISI_RESULTS_FILE:-$REPO_DIR/results/lisi/lisi_results.rds}"
RSCRIPT_DIR="${RSCRIPT_DIR:-/work/home/mengxu1310/miniforge3/bin}"
SKIP_HEAVY="${SKIP_HEAVY:-0}"
DRY_RUN=0
if [[ "${1:-}" == "--dry-run" ]]; then
  DRY_RUN=1
fi

export PATH="$RSCRIPT_DIR:$PATH"
export INTEGRATION_DIR
export OUT_DIR="$SCIENCEDB_OUT_DIR"
export OBJECT_FILE
export LISI_RESULTS_FILE
export SKIP_HEAVY

require_file() {
  local path="$1"
  if [[ ! -f "$path" ]]; then
    echo "Missing required file: $path" >&2
    exit 1
  fi
}

require_dir() {
  local path="$1"
  if [[ ! -d "$path" ]]; then
    echo "Missing required directory: $path" >&2
    exit 1
  fi
}

require_dir "$REPO_DIR"
require_dir "$INTEGRATION_DIR"
require_file "$INTEGRATION_DIR/metadata_filtered.rds"
require_file "$INTEGRATION_DIR/objects_celltype_plot.rds"
if [[ "$SKIP_HEAVY" != "1" ]]; then
  require_file "$OBJECT_FILE"
fi

mkdir -p "$SCIENCEDB_OUT_DIR"
cd "$REPO_DIR"

echo "REPO_DIR=$REPO_DIR"
echo "INTEGRATION_DIR=$INTEGRATION_DIR"
echo "SCIENCEDB_OUT_DIR=$SCIENCEDB_OUT_DIR"
echo "OBJECT_FILE=$OBJECT_FILE"
echo "LISI_RESULTS_FILE=$LISI_RESULTS_FILE"
echo "SKIP_HEAVY=$SKIP_HEAVY"

if [[ "$DRY_RUN" == "1" ]]; then
  echo "Dry run completed; ScienceDB export was not started."
  exit 0
fi

bash 05_sciencedb.sh
