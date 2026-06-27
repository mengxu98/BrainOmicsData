#!/bin/bash
#SBATCH -J brainomics_integration
#SBATCH -p xhacnormalb01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=512G
#SBATCH -t 72:00:00
#SBATCH -o /work/home/mengxu1310/repositories/BrainOmicsData/hpc/logs/integration_%j.out
#SBATCH -e /work/home/mengxu1310/repositories/BrainOmicsData/hpc/logs/integration_%j.err

set -euo pipefail

unset LC_ALL LC_CTYPE LANG LANGUAGE
export LANG=C

REPO_DIR="${REPO_DIR:-/work/home/mengxu1310/repositories/BrainOmicsData}"
RSCRIPT_DIR="${RSCRIPT_DIR:-/work/home/mengxu1310/miniforge3/bin}"
OVERWRITE="${OVERWRITE:-F}"

export PATH="$RSCRIPT_DIR:$PATH"

mkdir -p "$REPO_DIR/hpc/logs"
cd "$REPO_DIR"

bash 03_datasets_integration.sh "$OVERWRITE"
