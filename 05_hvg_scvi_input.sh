#!/usr/bin/env bash
# Select the 3,000 variable features from the common reference assembled in
# stage 04 and export scVI input shards to BRAINOMICS_RESULTS_DIR.
#
# Requirements: 512 GB memory; no GPU.
set -euo pipefail

BRAINOMICS_STAGE=05_hvg_scvi_input
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/functions/pipeline_lib.sh"

cd "$BRAINOMICS_REPO_ROOT"
brainomics_require_executor

overwrite="${1:-F}"
export BRAINOMICS_SCVI_INPUT_OVERWRITE="$overwrite"

if [ "$BRAINOMICS_EXECUTOR" = "slurm" ]; then
  brainomics_run_sbatch hpc/integration.sbatch "OVERWRITE=$overwrite"
  exit 0
fi

brainomics_require_command "$BRAINOMICS_RSCRIPT"
brainomics_require_file "$BRAINOMICS_RESULTS_DIR/objects_filtered.rds"
brainomics_require_file "$BRAINOMICS_RESULTS_DIR/objects_list_processed.rds"
"$BRAINOMICS_RSCRIPT" integration/datasets_scvi_input.R

brainomics_log "variable features and scVI input exported under $BRAINOMICS_RESULTS_DIR"
