#!/bin/bash
#SBATCH -J brainomics_age_schema
#SBATCH -p xhacnormalb01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=512G
#SBATCH -t 24:00:00
#SBATCH -o /work/home/mengxu1310/repositories/BrainOmicsData/hpc/logs/age_schema_%j.out
#SBATCH -e /work/home/mengxu1310/repositories/BrainOmicsData/hpc/logs/age_schema_%j.err

set -euo pipefail

unset LC_ALL LC_CTYPE LANG LANGUAGE
export LANG=C

REPO_DIR="${REPO_DIR:-/work/home/mengxu1310/repositories/BrainOmicsData}"
INTEGRATION_DIR="${INTEGRATION_DIR:-/work/home/mengxu1310/data/BrainOmicsData/integration}"
TARGETS="${TARGETS:-metadata_filtered.rds,objects_celltypes.rds,objects_celltype_plot.rds}"
RSCRIPT_BIN="${RSCRIPT_BIN:-/work/home/mengxu1310/miniforge3/bin/Rscript}"

mkdir -p "$REPO_DIR/hpc/logs"
cd "$REPO_DIR"

"$RSCRIPT_BIN" sciencedb/add_age_interval_schema_to_rds.R \
  --integration-dir "$INTEGRATION_DIR" \
  --targets "$TARGETS"
