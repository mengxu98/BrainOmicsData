#!/usr/bin/env bash
# Verify the pinned R and Python environments, or install them on the cluster.
#
# Requires the pinned R and Python environments in environment/ and network
# access when they are installed rather than verified.
set -euo pipefail

BRAINOMICS_STAGE=01_environment
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/functions/pipeline_lib.sh"
source "functions/utils.sh"

cd "$BRAINOMICS_REPO_ROOT"
brainomics_require_executor

if [ "$BRAINOMICS_EXECUTOR" = "slurm" ]; then
  brainomics_run_sbatch hpc/install_integration_dependencies.sbatch
  brainomics_run_sbatch hpc/install_scvi.sbatch
  exit 0
fi

check_command Rscript
check_command python3

Rscript environment/verify_environment.R "$BRAINOMICS_REPO_ROOT" preprocessing
Rscript environment/verify_environment.R "$BRAINOMICS_REPO_ROOT" integration
python3 environment/verify_python.py

brainomics_log "R and Python environments match the locked versions"
