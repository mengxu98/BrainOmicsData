#!/bin/bash
#SBATCH -J brainomics_scvi
# brainomics-partition: gpu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --mem=59G
#SBATCH -t 72:00:00

set -euo pipefail

repo_dir="${BRAINOMICS_REPO_ROOT:?run through functions/pipeline_lib.sh or set BRAINOMICS_REPO_ROOT}"
data_root="${BRAINOMICS_DATA_ROOT:?set BRAINOMICS_DATA_ROOT}"
integration_dir="$data_root/integration_25"
scvi_prefix="${BRAINOMICS_HPC_ROOT:?set BRAINOMICS_HPC_ROOT in hpc/local.env}/envs/brainomics-scvi"
runtime_lib_dir="${BRAINOMICS_HPC_ROOT:?set BRAINOMICS_HPC_ROOT in hpc/local.env}/software_sources/lib/brainomics-scvi-runtime"
overwrite="${OVERWRITE:-F}"
run_dir="$integration_dir/run"

export LD_LIBRARY_PATH="$runtime_lib_dir:${LD_LIBRARY_PATH:-}"
export PATH="$scvi_prefix/bin:$PATH"
export BRAINOMICS_EXECUTOR=local
if [ -n "${RUN_ROOT:-}" ]; then
  export BRAINOMICS_RUN_ROOT="$RUN_ROOT"
fi
export BRAINOMICS_SCVI_ACCELERATOR=gpu
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "$repo_dir/hpc/logs" "$run_dir"
cd "$repo_dir"

exit_file="$run_dir/scvi.exit"
temporary_exit="${exit_file}.tmp.$$"
rm -f "$exit_file"
set +e
bash 07_scvi.sh "$overwrite"
status="$?"
set -e
printf '%s\n' "$status" > "$temporary_exit"
mv "$temporary_exit" "$exit_file"
exit "$status"
