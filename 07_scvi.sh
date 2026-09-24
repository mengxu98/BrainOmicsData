#!/usr/bin/env bash
# Fit scVI on every raw-count reference cell and export the latent space.
#
# Input and output are read and written under BRAINOMICS_RESULTS_DIR.
#
# Requirements: 16 GB VRAM; 32 GB host memory.
set -euo pipefail

BRAINOMICS_STAGE=07_scvi
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/functions/pipeline_lib.sh"

cd "$BRAINOMICS_REPO_ROOT"
brainomics_require_executor

overwrite="${1:-F}"

if [ "$BRAINOMICS_EXECUTOR" = "slurm" ]; then
  brainomics_run_sbatch hpc/scvi.sbatch "OVERWRITE=$overwrite"
  exit 0
fi

brainomics_require_command "$BRAINOMICS_SCVI_PYTHON"
input_dir="$BRAINOMICS_RESULTS_DIR/scvi_input"
output_dir="$BRAINOMICS_RESULTS_DIR/scvi_output"
brainomics_require_dir "$input_dir"
mkdir -p "$output_dir"

scvi_command=(
  "$BRAINOMICS_SCVI_PYTHON" integration/run_scvi.py
  --input-dir "$input_dir"
  --output-dir "$output_dir"
  --seed "${BRAINOMICS_INTEGRATION_SEED:-20260730}"
  --latent-dimensions "${BRAINOMICS_INTEGRATION_DIMS:-50}"
  --max-epochs "${BRAINOMICS_SCVI_MAX_EPOCHS:-100}"
  --batch-size "${BRAINOMICS_SCVI_BATCH_SIZE:-2048}"
  --n-hidden "${BRAINOMICS_SCVI_N_HIDDEN:-128}"
  --n-layers "${BRAINOMICS_SCVI_N_LAYERS:-2}"
  --dropout-rate "${BRAINOMICS_SCVI_DROPOUT_RATE:-0.1}"
  --dispersion "${BRAINOMICS_SCVI_DISPERSION:-gene}"
  --gene-likelihood "${BRAINOMICS_SCVI_GENE_LIKELIHOOD:-nb}"
  --train-size "${BRAINOMICS_SCVI_TRAIN_SIZE:-0.9}"
  --validation-size "${BRAINOMICS_SCVI_VALIDATION_SIZE:-0.1}"
  --early-stopping-patience "${BRAINOMICS_SCVI_EARLY_STOPPING_PATIENCE:-15}"
  --kl-warmup-steps "${BRAINOMICS_SCVI_KL_WARMUP_STEPS:-400}"
  --learning-rate "${BRAINOMICS_SCVI_LEARNING_RATE:-0.001}"
  --accelerator "${BRAINOMICS_SCVI_ACCELERATOR:-gpu}"
  --devices "${BRAINOMICS_SCVI_DEVICES:-1}"
)
if brainomics_is_true "$overwrite"; then
  scvi_command+=(--overwrite)
fi

"${scvi_command[@]}"
brainomics_log "scVI latent space written to $output_dir"
