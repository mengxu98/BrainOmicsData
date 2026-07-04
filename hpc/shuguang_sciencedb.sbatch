#!/bin/bash
#SBATCH -J brainomics_sciencedb
#SBATCH -p cpu-partition
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=480G
#SBATCH -t 48:00:00
#SBATCH --nodelist=a01r05n03
#SBATCH -o /path/to/hpc/home/repositories/BrainOmicsData/hpc/logs/sciencedb_%j.out
#SBATCH -e /path/to/hpc/home/repositories/BrainOmicsData/hpc/logs/sciencedb_%j.err

set -euo pipefail

unset LC_ALL LC_CTYPE LANG LANGUAGE
export LANG=C

REPO_DIR="${REPO_DIR:-/path/to/hpc/home/repositories/BrainOmicsData}"
INTEGRATION_DIR="${INTEGRATION_DIR:-/path/to/hpc/home/data/BrainOmicsData/integration}"
OBJECT_FILE="${OBJECT_FILE:-$INTEGRATION_DIR/objects_celltypes.rds}"
METADATA_OBJECT_FILE="${METADATA_OBJECT_FILE:-$INTEGRATION_DIR/objects_celltype_plot.rds}"
OUT_DIR="${OUT_DIR:-/path/to/hpc/home/data/BrainOmicsData/ScienceDB}"
RSCRIPT_BIN="${RSCRIPT_BIN:-/path/to/hpc/home/miniforge3/bin/Rscript}"

mkdir -p "$REPO_DIR/hpc/logs"
cd "$REPO_DIR"

REPO_DIR="$REPO_DIR" \
INTEGRATION_DIR="$INTEGRATION_DIR" \
OBJECT_FILE="$OBJECT_FILE" \
METADATA_OBJECT_FILE="$METADATA_OBJECT_FILE" \
OUT_DIR="$OUT_DIR" \
RSCRIPT_BIN="$RSCRIPT_BIN" \
bash 05_sciencedb.sh
