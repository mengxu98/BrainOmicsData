#!/usr/bin/env bash
# Serialize the frozen results into the deposit package: 22 expression shards,
# 22-column metadata, eight embeddings, per-cell LISI and the provenance tables.
#
# Requirements: 100 GB memory; no GPU.
set -euo pipefail

BRAINOMICS_STAGE=11_sciencedb_export
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/functions/pipeline_lib.sh"
source "functions/utils.sh"

cd "$BRAINOMICS_REPO_ROOT"
brainomics_require_executor
check_command Rscript

frozen="${BRAINOMICS_FROZEN_DIR:-$BRAINOMICS_REPO_ROOT/results/frozen_run}"
analysis="${BRAINOMICS_ANALYSIS_DIR:-$BRAINOMICS_REPO_ROOT/results/analysis_run}"
package="${BRAINOMICS_PACKAGE_DIR:-$BRAINOMICS_DATA_ROOT/ScienceDB}"
work="${BRAINOMICS_WORK_DIR:-$BRAINOMICS_REPO_ROOT/results/work}"
lisi_dir="${BRAINOMICS_LISI_DIR:-$BRAINOMICS_RESULTS_DIR/evaluation/lisi_checkpoints_double_precision}"

brainomics_require_dir "$frozen"
brainomics_require_dir "$analysis"
brainomics_require_dir "$lisi_dir"
mkdir -p "$work"

if [ ! -f "$work/percent_mito_by_cell.tsv.gz" ]; then
  brainomics_log "computing percent_mito from the original source counts"
  bash processing/percent_mito_build.sh "$BRAINOMICS_DATA_ROOT" "$work" "$analysis"
fi

if [ ! -f "$work/final_cell_annotations.tsv.gz" ]; then
  brainomics_log "exporting the adopted per-cell annotation table"
  Rscript functions/export_cell_annotations.R "$frozen" "$analysis" "$work"
fi

for input in percent_mito_by_cell.tsv.gz final_cell_annotations.tsv.gz; do
  brainomics_require_file "$work/$input"
done

Rscript functions/export_package.R "$frozen" "$analysis" "$package" "$work"
Rscript functions/export_lisi.R "$analysis" "$package" "$work" "$lisi_dir"

brainomics_log "deposit package exported and verified at $package"
python3 functions/export_readme.py "$package"
python3 functions/export_feature_metadata.py . "$package"
python3 functions/merge_dataset_manifest.py . "$package"
python3 functions/normalize_package.py "$package" "$work/package_seal.json"
