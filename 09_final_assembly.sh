#!/usr/bin/env bash
# Collect the completed R reductions, import scVI, validate the common-input
# comparison and evaluate the full-cell latent space.
#
# All integration methods must finish before this stage begins.
#
# Requirements: 100 GB memory; no GPU.
set -euo pipefail

BRAINOMICS_STAGE=09_final_assembly
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/functions/pipeline_lib.sh"

cd "$BRAINOMICS_REPO_ROOT"
brainomics_require_executor

overwrite="${1:-F}"
export BRAINOMICS_EVALUATION_OVERWRITE="$overwrite"

if [ "$BRAINOMICS_EXECUTOR" = "slurm" ]; then
  # hpc/integration_finalize.sbatch runs this whole stage, both the
  # finalize and the evaluation phase.
  brainomics_run_sbatch hpc/integration_finalize.sbatch \
    "OVERWRITE=$overwrite"
  exit 0
fi

brainomics_require_command "$BRAINOMICS_RSCRIPT"
for method in raw_umap harmony rpca; do
  brainomics_require_file "$BRAINOMICS_RESULTS_DIR/r_checkpoints/${method}_reductions.rds"
done
for artifact in scvi_output_audit.tsv cell_ids.txt scvi_latent.float32.bin scvi_versions.tsv; do
  brainomics_require_file "$BRAINOMICS_RESULTS_DIR/scvi_output/$artifact"
done
"$BRAINOMICS_RSCRIPT" functions/run_parallel_integration.R "$BRAINOMICS_RUN_ROOT" collect
"$BRAINOMICS_RSCRIPT" integration/datasets_integration_04.R
"$BRAINOMICS_RSCRIPT" integration/datasets_integration_05.R

brainomics_log "integrated object and latent-space evaluation written to $BRAINOMICS_RESULTS_DIR"
