#!/usr/bin/env bash
# Raw PCA, raw UMAP and Harmony on the shared 3,000-HVG
# input. Each stage requires the previous stage's checkpoints.
#
# Checkpoints are read and written under BRAINOMICS_RESULTS_DIR.
#
# Requirements: 200 GB memory; 32-64 cores; no GPU.
set -euo pipefail

BRAINOMICS_STAGE=06_pca_harmony
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/functions/pipeline_lib.sh"

cd "$BRAINOMICS_REPO_ROOT"
brainomics_require_executor
brainomics_require_file "$BRAINOMICS_RESULTS_DIR/r_checkpoints/objects_hvg.rds"

for stage in pca raw_umap harmony; do
  if [ "$BRAINOMICS_EXECUTOR" = "slurm" ]; then
    brainomics_run_sbatch hpc/parallel_integration.sbatch \
      "STAGE=$stage"
  else
    brainomics_log "running stage $stage"
    brainomics_require_command "$BRAINOMICS_RSCRIPT"
    "$BRAINOMICS_RSCRIPT" "$BRAINOMICS_REPO_ROOT/functions/run_parallel_integration.R" \
      "$BRAINOMICS_RUN_ROOT" "$stage"
  fi
done

brainomics_log "raw and Harmony reductions written to $BRAINOMICS_RESULTS_DIR"
