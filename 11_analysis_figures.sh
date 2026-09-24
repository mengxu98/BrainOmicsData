#!/usr/bin/env bash
# Draw manuscript figures from the current frozen analysis run.
# The former table-refresh phase used a separate integration directory with
# historical cohort results. Current source tables live in results/analysis_run.
#
# Requirements: 50 GB memory; no GPU.
set -euo pipefail

BRAINOMICS_STAGE=10_annotation_figures
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/functions/pipeline_lib.sh"
source "functions/utils.sh"

cd "$BRAINOMICS_REPO_ROOT"
brainomics_require_executor

if [ "$BRAINOMICS_EXECUTOR" = "slurm" ]; then
  brainomics_die "annotation and figure refresh runs on the storage host; use BRAINOMICS_EXECUTOR=local"
fi

check_command Rscript

# Figures 1-5 and their panels, read from the frozen analysis run.
Rscript --vanilla plotting/fig1.R
Rscript --vanilla plotting/fig2.R
Rscript --vanilla plotting/fig3.R
Rscript --vanilla plotting/fig4.R
Rscript --vanilla plotting/fig5.R
Rscript --vanilla plotting/figS1.R
Rscript --vanilla plotting/figS2.R
Rscript --vanilla plotting/figS3.R

brainomics_log "main, S1 and S2 figures regenerated; S3 PNG exported from its existing PDF; S4 is rendered locally"
