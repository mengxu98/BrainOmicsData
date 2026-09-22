#!/bin/bash
#SBATCH -J brainomics_scvi
#SBATCH -p gpu-partition
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --mem=59G
#SBATCH -t 72:00:00
#SBATCH -o /path/to/hpc/home/repositories/BrainOmicsData/hpc/logs/scvi_%j.out
#SBATCH -e /path/to/hpc/home/repositories/BrainOmicsData/hpc/logs/scvi_%j.err

set -euo pipefail

repo_dir="/path/to/hpc/home/repositories/BrainOmicsData"
data_root="/path/to/hpc/home/data/BrainOmicsData"
integration_dir="$data_root/integration_25"
scvi_prefix="/path/to/hpc/home/envs/brainomics-scvi"
runtime_lib_dir="/path/to/hpc/home/software_sources/lib/brainomics-scvi-runtime"
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
