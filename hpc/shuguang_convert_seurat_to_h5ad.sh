#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${ENV_NAME:-brainomics-h5ad}"
CONDA_BIN="${CONDA_BIN:-/work/home/mengxu1310/miniforge3/bin/conda}"
REPO_DIR="${REPO_DIR:-/work/home/mengxu1310/repositories/BrainOmicsData}"
OBJECT_FILE="${OBJECT_FILE:-/work/home/mengxu1310/data/BrainOmicsData/integration/objects_celltypes.rds}"
OUT_FILE="${OUT_FILE:-/work/home/mengxu1310/data/ScienceDB/HumanBrain_sc_snRNAseq_age_resource/04_processed_objects/human_brain_integrated_full.h5ad}"

ENV_PREFIX="$("$CONDA_BIN" env list | awk -v env="$ENV_NAME" '$1 == env {print $NF}')"
if [[ -z "$ENV_PREFIX" || ! -d "$ENV_PREFIX" ]]; then
  echo "Missing conda environment: $ENV_NAME" >&2
  echo "Run hpc/shuguang_setup_h5ad_env.sh first." >&2
  exit 1
fi

export RETICULATE_PYTHON="$ENV_PREFIX/bin/python"

cd "$REPO_DIR"
"$ENV_PREFIX/bin/Rscript" hpc/shuguang_convert_seurat_to_h5ad.R \
  --object-file "$OBJECT_FILE" \
  --out-file "$OUT_FILE"
