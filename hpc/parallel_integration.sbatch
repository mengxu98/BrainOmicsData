#!/bin/bash
#SBATCH -J brainomics_parallel_integration
# brainomics-partition: cpu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=480G
#SBATCH -t 72:00:00

# One stage of the shared 3,000-HVG reductions. STAGE is pca, raw_umap,
# harmony, collect or rpca; RUN_ROOT is the run directory holding pipeline/.
set -euo pipefail

repo_dir="${BRAINOMICS_REPO_ROOT:?run through functions/pipeline_lib.sh or set BRAINOMICS_REPO_ROOT}"
stage="${STAGE:?Supply STAGE (pca, raw_umap, harmony, collect or rpca)}"
run_root="${RUN_ROOT:?Supply RUN_ROOT}"
overwrite="${OVERWRITE:-F}"

mkdir -p "$repo_dir/hpc/logs"

export OVERWRITE="$overwrite"
export BRAINOMICS_EXECUTOR=local
exec bash "$repo_dir/integration/parallel_integration_stage.sh" "$stage" "$run_root"
