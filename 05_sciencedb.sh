#!/usr/bin/env bash
set -euo pipefail

source "functions/utils.sh"

check_command Rscript
check_command python3

REPO_DIR="$(pwd)"
INTEGRATION_DIR="${INTEGRATION_DIR:-../../data/BrainOmicsData/integration}"
OUT_DIR="${OUT_DIR:-../../data/ScienceDB/HumanBrain_sc_snRNAseq_age_resource}"
OBJECT_FILE="${OBJECT_FILE:-$INTEGRATION_DIR/objects_celltypes.rds}"
LISI_RESULTS_FILE="${LISI_RESULTS_FILE:-results/lisi/lisi_results.rds}"
SKIP_HEAVY="${SKIP_HEAVY:-0}"

if [[ -n "${R_LIB_PATH:-}" ]]; then
  export LD_LIBRARY_PATH="$R_LIB_PATH:${LD_LIBRARY_PATH:-}"
fi

mkdir -p "$OUT_DIR"

log_message "Preparing ScienceDB package at {.file ${OUT_DIR}}..."

log_message "Exporting source, sample and lightweight cell metadata..."
Rscript sciencedb/export_sciencedb_metadata.R \
  --integration-dir "$INTEGRATION_DIR" \
  --out-dir "$OUT_DIR" \
  --repo-dir "$REPO_DIR"

if [[ "$SKIP_HEAVY" != "1" ]]; then
  log_message "Exporting heavy integrated-object tables from {.file ${OBJECT_FILE}}..."
  Rscript sciencedb/export_sciencedb_integrated_object.R \
    --object-file "$OBJECT_FILE" \
    --out-dir "$OUT_DIR"
  log_message "Staging processed Seurat objects..."
  bash sciencedb/stage_processed_objects.sh
else
  log_message "Skipping heavy integrated-object export because SKIP_HEAVY=1." --message-type warning
fi

log_message "Exporting LISI integration-quality summary..."
Rscript sciencedb/export_sciencedb_lisi_summary.R \
  --lisi-results-file "$LISI_RESULTS_FILE" \
  --out-dir "$OUT_DIR"

log_message "Writing figure and table source data..."
Rscript sciencedb/plot_sciencedb_resource_figures.R \
  --out-dir "$OUT_DIR"

workflow_dir="$OUT_DIR/08_reusable_workflow"
scripts_dir="$workflow_dir/scripts"
mkdir -p "$scripts_dir"

log_message "Copying reusable workflow scripts..."
rsync -a \
  --exclude '.DS_Store' \
  --exclude '__pycache__' \
  01_datasets_download.sh \
  02_datasets_preprocessing.sh \
  03_datasets_integration.sh \
  04_datasets_plotting.sh \
  05_sciencedb.sh \
  download \
  processing \
  integration \
  plotting \
  sciencedb \
  functions \
  "$scripts_dir/"

cat > "$workflow_dir/example_commands.md" <<'EOF'
# Example commands

Run from the BrainOmicsData repository root:

```sh
bash 01_datasets_download.sh
bash 02_datasets_preprocessing.sh
bash 03_datasets_integration.sh
bash 04_datasets_plotting.sh
bash 05_sciencedb.sh
```

For a lightweight metadata-only ScienceDB staging run:

```sh
SKIP_HEAVY=1 bash 05_sciencedb.sh
```

To force regeneration of upstream integration or plotting outputs:

```sh
bash 03_datasets_integration.sh T
bash 04_datasets_plotting.sh T
```
EOF

if command -v conda >/dev/null 2>&1; then
  conda env export > "$workflow_dir/conda_environment.yml" || true
fi

log_message "Writing file manifest and field dictionary..."
Rscript sciencedb/write_sciencedb_manifest.R \
  --out-dir "$OUT_DIR"

python3 sciencedb/write_xlsx_from_tsv.py \
  --input "$OUT_DIR/00_README/data_dictionary.tsv" \
  --output "$OUT_DIR/00_README/data_dictionary.xlsx"

Rscript sciencedb/write_sciencedb_manifest.R \
  --out-dir "$OUT_DIR"

log_message "ScienceDB package prepared at {.file ${OUT_DIR}}" --message-type success
