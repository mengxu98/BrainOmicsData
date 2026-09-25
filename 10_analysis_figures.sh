#!/usr/bin/env bash
# Draw all main and supplementary figures from the analysis outputs.
# The statistical source tables must be prepared before this stage.
#
# Requirements: 50 GB memory; no GPU.
set -euo pipefail

BRAINOMICS_STAGE=10_analysis_figures
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/functions/pipeline_lib.sh"
source "functions/utils.sh"

cd "$BRAINOMICS_REPO_ROOT"
brainomics_require_executor

if [ "$BRAINOMICS_EXECUTOR" = "slurm" ]; then
  brainomics_die "figure assembly uses local analysis inputs; set BRAINOMICS_EXECUTOR=local"
fi

check_command "$BRAINOMICS_RSCRIPT"

# Figures 1-5 and S1-S4 use the same selected R environment.
"$BRAINOMICS_RSCRIPT" --vanilla plotting/fig1.R
"$BRAINOMICS_RSCRIPT" --vanilla plotting/fig2.R
"$BRAINOMICS_RSCRIPT" --vanilla plotting/fig3.R
"$BRAINOMICS_RSCRIPT" --vanilla plotting/fig4.R
"$BRAINOMICS_RSCRIPT" --vanilla plotting/fig5.R
"$BRAINOMICS_RSCRIPT" --vanilla plotting/figS1.R
"$BRAINOMICS_RSCRIPT" --vanilla plotting/figS2.R
"$BRAINOMICS_RSCRIPT" --vanilla plotting/figS3.R --redraw
"$BRAINOMICS_RSCRIPT" --vanilla plotting/figS4.R

brainomics_log "main and supplementary figures generated from the supplied analysis inputs"
