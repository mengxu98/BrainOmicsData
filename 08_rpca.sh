#!/usr/bin/env bash
# Robust PCA on the shared 3,000-HVG / 50-PC input, batch model "dataset",
# seed 20260730, dimensions 1-30 for the anchors.
#
# Requires BRAINOMICS_RUN_ROOT: the run directory whose pipeline/ directory holds
# this repository and whose integration/ directory holds the checkpoints.
#
# Requirements: 1 TB memory; no GPU.
set -euo pipefail

BRAINOMICS_STAGE=08_rpca
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/functions/pipeline_lib.sh"

cd "$BRAINOMICS_REPO_ROOT"
brainomics_require_executor
[ -n "${BRAINOMICS_RUN_ROOT:-}" ] ||
  brainomics_die "set BRAINOMICS_RUN_ROOT to the run directory holding pipeline/ and integration/"
export BRAINOMICS_RUN_ROOT

if [ "$BRAINOMICS_EXECUTOR" = "slurm" ]; then
  brainomics_run_sbatch hpc/parallel_integration.sbatch \
    "STAGE=rpca" "RUN_ROOT=$BRAINOMICS_RUN_ROOT"
  exit 0
fi

brainomics_log "running the RPCA stage"
Rscript "$BRAINOMICS_REPO_ROOT/functions/run_parallel_integration.R" \
  "$BRAINOMICS_RUN_ROOT" rpca

brainomics_log "RPCA reduction written to $BRAINOMICS_RUN_ROOT/integration"
