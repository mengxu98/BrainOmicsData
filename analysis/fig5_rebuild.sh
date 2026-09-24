#!/usr/bin/env bash
# Recompute Figure 5 from the current frozen atlas and processed source objects.
set -euo pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo"
analysis_dir="${BRAINOMICS_ANALYSIS_DIR:-$repo/results/analysis_run}"
metadata="$analysis_dir/07_downstream/revision_20260918/01_metadata/metadata_working.rds"
test -f "$metadata"

Rscript --vanilla -e 'x <- readRDS(commandArgs(TRUE)[1]); writeLines(sort(unique(x$Dataset)))' "$metadata" |
while IFS= read -r dataset; do
  test -n "$dataset" || continue
  Rscript --vanilla analysis/fig5_extract_full_panel.R "$dataset"
done

Rscript --vanilla analysis/fig5_prepare_sources.R
Rscript --vanilla plotting/fig5.R
