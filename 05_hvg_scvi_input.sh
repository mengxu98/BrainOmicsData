#!/usr/bin/env bash
# Assemble the filtered and merged reference objects, select the 3,000 variable
# features and export the scVI input shards.
#
# Optional BRAINOMICS_RUN_ROOT: when set, the stage runs inside
# $BRAINOMICS_RUN_ROOT/pipeline instead of this checkout.
#
# Requirements: 512 GB memory; no GPU.
set -euo pipefail

BRAINOMICS_STAGE=05_hvg_scvi_input
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/functions/pipeline_lib.sh"
source "functions/utils.sh"

brainomics_require_executor

overwrite="${1:-F}"
export BRAINOMICS_SCVI_INPUT_OVERWRITE="$overwrite"

code_root="$BRAINOMICS_REPO_ROOT"
if [ -n "${BRAINOMICS_RUN_ROOT:-}" ]; then
  code_root="$BRAINOMICS_RUN_ROOT/pipeline"
  brainomics_require_dir "$code_root"
fi
cd "$code_root"

if [ "$BRAINOMICS_EXECUTOR" = "slurm" ]; then
  run_root_spec=()
  if [ -n "${BRAINOMICS_RUN_ROOT:-}" ]; then
    run_root_spec=("RUN_ROOT=$BRAINOMICS_RUN_ROOT")
  fi
  brainomics_run_sbatch hpc/integration.sbatch \
    "OVERWRITE=$overwrite" "${run_root_spec[@]+"${run_root_spec[@]}"}"
  exit 0
fi

check_command Rscript
Rscript integration/datasets_scvi_input.R

brainomics_log "variable features and scVI input exported under $BRAINOMICS_RESULTS_DIR"
