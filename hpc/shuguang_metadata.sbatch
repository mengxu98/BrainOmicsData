#!/bin/bash
#SBATCH -J brainomics_metadata
#SBATCH -p cpu-partition
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=480G
#SBATCH -t 24:00:00
#SBATCH -o /path/to/hpc/home/repositories/BrainOmicsData/hpc/logs/metadata_%j.out
#SBATCH -e /path/to/hpc/home/repositories/BrainOmicsData/hpc/logs/metadata_%j.err

set -euo pipefail

unset LC_ALL LC_CTYPE LANG LANGUAGE
export LANG=C

REPO_DIR="${REPO_DIR:-/path/to/hpc/home/repositories/BrainOmicsData}"
INTEGRATION_DIR="${INTEGRATION_DIR:-/path/to/hpc/home/data/BrainOmicsData/integration}"
TARGETS="${TARGETS:-metadata_filtered.rds,objects_celltypes.rds,objects_celltype_plot.rds}"
RSCRIPT_BIN="${RSCRIPT_BIN:-/path/to/hpc/home/miniforge3/bin/Rscript}"
BACKUP_DIR="${BACKUP_DIR:-$INTEGRATION_DIR/schema_backups}"

mkdir -p "$REPO_DIR/hpc/logs"
cd "$REPO_DIR"

"$RSCRIPT_BIN" sciencedb/metadata.R \
  --integration-dir "$INTEGRATION_DIR" \
  --targets "$TARGETS" \
  --backup-dir "$BACKUP_DIR"
