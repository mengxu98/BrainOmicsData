#!/usr/bin/env bash
# Raw PCA, raw UMAP, Harmony and checkpoint collection on the shared 3,000-HVG
# input. Each stage requires the previous stage's checkpoints.
#
# Requires BRAINOMICS_RUN_ROOT: the run directory whose pipeline/ directory holds
# this repository and whose integration/ directory holds the checkpoints.
#
# Requirements: 200 GB memory; 32-64 cores; no GPU.
set -euo pipefail

BRAINOMICS_STAGE=06_pca_harmony
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/functions/pipeline_lib.sh"

cd "$BRAINOMICS_REPO_ROOT"
brainomics_require_executor
[ -n "${BRAINOMICS_RUN_ROOT:-}" ] ||
  brainomics_die "set BRAINOMICS_RUN_ROOT to the run directory holding pipeline/ and integration/"
export BRAINOMICS_RUN_ROOT

for stage in pca raw_umap harmony collect; do
  if [ "$BRAINOMICS_EXECUTOR" = "hpc" ]; then
    brainomics_run_sbatch hpc/hpc_parallel_integration.sbatch \
      "STAGE=$stage" "RUN_ROOT=$BRAINOMICS_RUN_ROOT"
  else
    brainomics_log "running stage $stage"
    Rscript "$BRAINOMICS_REPO_ROOT/functions/run_parallel_integration.R" \
      "$BRAINOMICS_RUN_ROOT" "$stage"
  fi
done

brainomics_log "raw, Harmony and collected reductions written to $BRAINOMICS_RUN_ROOT/integration"
