#!/bin/bash
#SBATCH -J brainomics_scvi_input
# brainomics-partition: cpu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=480G
#SBATCH -t 72:00:00

set -euo pipefail

unset LC_ALL LC_CTYPE LANGUAGE CPLUS_INCLUDE_PATH C_INCLUDE_PATH
export LANG=C
module load mathlib/hdf5/1.8.20-gcc

repo_dir="${BRAINOMICS_REPO_ROOT:?run through functions/pipeline_lib.sh or set BRAINOMICS_REPO_ROOT}"
r_prefix="${BRAINOMICS_HPC_ROOT:?set BRAINOMICS_HPC_ROOT in hpc/local.env}/miniforge3"
data_root="${BRAINOMICS_DATA_ROOT:?set BRAINOMICS_DATA_ROOT}"
integration_dir="$data_root/integration_25"
overwrite="${OVERWRITE:-F}"
run_dir="$integration_dir/run"

export PATH="$r_prefix/bin:$PATH"
export R_LIBS_USER="${R_LIBS_USER:-${BRAINOMICS_HPC_ROOT:?set BRAINOMICS_HPC_ROOT in hpc/local.env}/R/library/4.5}"
export BRAINOMICS_EXECUTOR=local
if [ -n "${RUN_ROOT:-}" ]; then
  export BRAINOMICS_RUN_ROOT="$RUN_ROOT"
fi
export BRAINOMICS_REMOVE_INTEGRATION_INTERMEDIATES=false

mkdir -p "$repo_dir/hpc/logs" "$run_dir"
cd "$repo_dir"

exit_file="$run_dir/scvi_input.exit"
temporary_exit="${exit_file}.tmp.$$"
rm -f "$exit_file"
set +e
bash 05_hvg_scvi_input.sh "$overwrite"
status="$?"
set -e
printf '%s\n' "$status" > "$temporary_exit"
mv "$temporary_exit" "$exit_file"
exit "$status"
