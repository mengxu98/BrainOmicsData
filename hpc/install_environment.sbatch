#!/bin/bash
#SBATCH -J brainomics_environment
# brainomics-partition: cpu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -t 08:00:00

set -euo pipefail
repo_dir="${BRAINOMICS_REPO_ROOT:?run 01_environment.sh --install or set BRAINOMICS_REPO_ROOT}"
# Uses the same public installer and environment prefix as a local Linux run.
# Compilers and HDF5 come from the shared prefix; no site module is required.
exec bash "$repo_dir/environment/restore.sh"
