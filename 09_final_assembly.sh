#!/usr/bin/env bash
# Import scVI into the assembled object, validate the common-input comparison and
# evaluate the full-cell latent space.
#
# Optional BRAINOMICS_RUN_ROOT: when set, the stage runs inside
# $BRAINOMICS_RUN_ROOT/pipeline instead of this checkout.
#
# Requirements: 100 GB memory; no GPU.
set -euo pipefail

BRAINOMICS_STAGE=09_final_assembly
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/functions/pipeline_lib.sh"
source "functions/utils.sh"

brainomics_require_executor

overwrite="${1:-F}"
export BRAINOMICS_EVALUATION_OVERWRITE="$overwrite"

code_root="$BRAINOMICS_REPO_ROOT"
if [ -n "${BRAINOMICS_RUN_ROOT:-}" ]; then
  code_root="$BRAINOMICS_RUN_ROOT/pipeline"
  brainomics_require_dir "$code_root"
fi
cd "$code_root"

if [ "$BRAINOMICS_EXECUTOR" = "hpc" ]; then
  run_root_spec=()
  if [ -n "${BRAINOMICS_RUN_ROOT:-}" ]; then
    run_root_spec=("RUN_ROOT=$BRAINOMICS_RUN_ROOT")
  fi
  # hpc/hpc_integration_finalize.sbatch runs this whole stage, both the
  # finalize and the evaluation phase.
  brainomics_run_sbatch hpc/hpc_integration_finalize.sbatch \
    "OVERWRITE=$overwrite" "${run_root_spec[@]+"${run_root_spec[@]}"}"
  exit 0
fi

check_command Rscript
Rscript integration/datasets_integration_04.R
Rscript integration/datasets_integration_05.R

brainomics_log "integrated object and latent-space evaluation written under $code_root"
