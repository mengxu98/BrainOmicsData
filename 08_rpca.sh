#!/usr/bin/env bash
# Reciprocal PCA on the shared 3,000-HVG / 50-PC input, batch model "dataset",
# seed 20260730, dimensions 1-30 for the anchors.
#
# Checkpoints are read and written under BRAINOMICS_RESULTS_DIR.
#
# Requirements: 1 TB memory; no GPU.
set -euo pipefail

BRAINOMICS_STAGE=08_rpca
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/functions/pipeline_lib.sh"

cd "$BRAINOMICS_REPO_ROOT"
brainomics_require_executor
brainomics_require_file "$BRAINOMICS_RESULTS_DIR/r_checkpoints/pca_input.rds"
brainomics_require_file "$BRAINOMICS_RESULTS_DIR/r_checkpoints/objects_pca.rds"

if [ "$BRAINOMICS_EXECUTOR" = "slurm" ]; then
  BRAINOMICS_SLURM_MEMORY="${BRAINOMICS_RPCA_MEMORY:-1T}" \
  brainomics_run_sbatch hpc/parallel_integration.sbatch \
    "STAGE=rpca"
  exit 0
fi

brainomics_log "running the RPCA stage"
brainomics_require_command "$BRAINOMICS_RSCRIPT"
"$BRAINOMICS_RSCRIPT" "$BRAINOMICS_REPO_ROOT/functions/run_parallel_integration.R" \
  "$BRAINOMICS_RUN_ROOT" rpca

brainomics_log "RPCA reduction written to $BRAINOMICS_RESULTS_DIR"
