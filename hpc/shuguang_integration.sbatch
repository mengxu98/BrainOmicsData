#!/bin/bash
#SBATCH -J brainomics_integration
#SBATCH -p cpu-partition
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=480G
#SBATCH -t 72:00:00
#SBATCH -o /path/to/hpc/home/repositories/BrainOmicsData/hpc/logs/integration_%j.out
#SBATCH -e /path/to/hpc/home/repositories/BrainOmicsData/hpc/logs/integration_%j.err

set -euo pipefail

unset LC_ALL LC_CTYPE LANG LANGUAGE
export LANG=C

REPO_DIR="${REPO_DIR:-/path/to/hpc/home/repositories/BrainOmicsData}"
RSCRIPT_DIR="${RSCRIPT_DIR:-/path/to/hpc/home/miniforge3/bin}"
OVERWRITE="${OVERWRITE:-F}"

export PATH="$RSCRIPT_DIR:$PATH"

mkdir -p "$REPO_DIR/hpc/logs"
cd "$REPO_DIR"

bash 03_datasets_integration.sh "$OVERWRITE"
