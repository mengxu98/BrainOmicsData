#!/bin/bash
#SBATCH -J brainomics_finalize
#SBATCH -p cpu-partition
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=480G
#SBATCH -t 48:00:00
#SBATCH -o /path/to/hpc/home/repositories/BrainOmicsData/hpc/logs/finalize_%j.out
#SBATCH -e /path/to/hpc/home/repositories/BrainOmicsData/hpc/logs/finalize_%j.err

set -euo pipefail

unset LC_ALL LC_CTYPE LANGUAGE CPLUS_INCLUDE_PATH C_INCLUDE_PATH
export LANG=C
module load mathlib/hdf5/1.8.20-gcc

repo_dir="/path/to/hpc/home/repositories/BrainOmicsData"
r_prefix="/path/to/hpc/home/miniforge3"
data_root="/path/to/hpc/home/data/BrainOmicsData"
integration_dir="$data_root/integration_25"
overwrite="${OVERWRITE:-F}"
run_dir="$integration_dir/run"

export PATH="$r_prefix/bin:$PATH"
export R_LIBS_USER="${R_LIBS_USER:-/path/to/hpc/home/R/library/4.5}"
export BRAINOMICS_EXECUTOR=local
if [ -n "${RUN_ROOT:-}" ]; then
  export BRAINOMICS_RUN_ROOT="$RUN_ROOT"
fi
export BRAINOMICS_REMOVE_INTEGRATION_INTERMEDIATES=true

mkdir -p "$repo_dir/hpc/logs" "$run_dir"
cd "$repo_dir"

exit_file="$run_dir/finalize.exit"
temporary_exit="${exit_file}.tmp.$$"
rm -f "$exit_file"
set +e
bash 09_final_assembly.sh "$overwrite"
status="$?"
set -e
printf '%s\n' "$status" > "$temporary_exit"
mv "$temporary_exit" "$exit_file"
exit "$status"
