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

repo_dir="${BRAINOMICS_REPO_ROOT:?run through functions/pipeline_lib.sh or set BRAINOMICS_REPO_ROOT}"
source "$repo_dir/functions/pipeline_lib.sh"
overwrite="${OVERWRITE:-F}"
run_dir="$BRAINOMICS_RESULTS_DIR/run"

export BRAINOMICS_EXECUTOR=local
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
