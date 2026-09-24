#!/bin/bash
#SBATCH -J brainomics_parallel_integration
# brainomics-partition: cpu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=480G
#SBATCH -t 72:00:00

# One stage of the shared 3,000-HVG reductions. STAGE is pca, raw_umap,
# harmony, collect or rpca. Source code is read from the current checkout.
set -euo pipefail

repo_dir="${BRAINOMICS_REPO_ROOT:?run through functions/pipeline_lib.sh or set BRAINOMICS_REPO_ROOT}"
stage="${STAGE:?Supply STAGE (pca, raw_umap, harmony, collect or rpca)}"
source "$repo_dir/functions/pipeline_lib.sh"
export BRAINOMICS_EXECUTOR=local
cd "$repo_dir"
brainomics_require_command "$BRAINOMICS_RSCRIPT"
exec "$BRAINOMICS_RSCRIPT" "$repo_dir/functions/run_parallel_integration.R" \
  "$BRAINOMICS_RUN_ROOT" "$stage"
