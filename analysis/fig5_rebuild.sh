#!/usr/bin/env bash
# Recompute Figure 5 from the prepared atlas and processed source objects.
set -euo pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo"
analysis_dir="${BRAINOMICS_ANALYSIS_DIR:-$repo/results/analysis_run}"
rscript="${BRAINOMICS_RSCRIPT:-Rscript}"
metadata="${BRAINOMICS_METADATA_FILE:-}"
if [ -z "$metadata" ]; then
  metadata="$("$rscript" --vanilla -e 'root<-commandArgs(TRUE)[1];x<-list.files(root,recursive=TRUE,full.names=TRUE);x<-x[basename(x)=="metadata_working.rds"];if(length(x)!=1L)stop("Expected one metadata_working.rds under ",root,"; found ",length(x));cat(x)' "$analysis_dir")"
fi
test -f "$metadata"

"$rscript" --vanilla -e 'x <- readRDS(commandArgs(TRUE)[1]); writeLines(sort(unique(x$Dataset)))' "$metadata" |
while IFS= read -r dataset; do
  test -n "$dataset" || continue
  "$rscript" --vanilla analysis/fig5_extract_full_panel.R "$dataset"
done

"$rscript" --vanilla analysis/fig5_prepare_sources.R
"$rscript" --vanilla plotting/fig5.R
