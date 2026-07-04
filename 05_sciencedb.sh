#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash 05_sciencedb.sh [--skip-export] [--skip-anonymize]

Options:
  --skip-export      Do not regenerate expression/metadata files; only run the
                     anonymization and manifest/MD5 refresh step on OUT_DIR.
  --skip-anonymize   Export the package without the public anonymization step.
  -h, --help         Show this help.

Environment:
  REPO_DIR, INTEGRATION_DIR, OBJECT_FILE, METADATA_OBJECT_FILE, OUT_DIR, RSCRIPT_BIN
EOF
}

REPO_DIR="${REPO_DIR:-$(pwd)}"
INTEGRATION_DIR="${INTEGRATION_DIR:-../../data/BrainOmicsData/integration}"
OBJECT_FILE="${OBJECT_FILE:-$INTEGRATION_DIR/objects_celltypes.rds}"
METADATA_OBJECT_FILE="${METADATA_OBJECT_FILE:-$INTEGRATION_DIR/objects_celltype_plot.rds}"
OUT_DIR="${OUT_DIR:-../../data/BrainOmicsData/ScienceDB}"
RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"
SKIP_EXPORT=0
SKIP_ANONYMIZE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --skip-export)
      SKIP_EXPORT=1
      shift
      ;;
    --skip-anonymize)
      SKIP_ANONYMIZE=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

source "functions/utils.sh"

check_command "$RSCRIPT_BIN"

if [[ "$SKIP_EXPORT" -eq 0 ]]; then
  log_message "Exporting compact ScienceDB package at {.file ${OUT_DIR}}..."

  "$RSCRIPT_BIN" sciencedb/export.R \
    --repo-dir "$REPO_DIR" \
    --object-file "$OBJECT_FILE" \
    --metadata-object-file "$METADATA_OBJECT_FILE" \
    --out-dir "$OUT_DIR" \
    --overwrite
else
  log_message "Skipping ScienceDB export; using existing package at {.file ${OUT_DIR}}"
fi

if [[ "$SKIP_ANONYMIZE" -eq 0 ]]; then
  log_message "Applying public-package anonymization at {.file ${OUT_DIR}}..."
  "$RSCRIPT_BIN" sciencedb/anonymize_public_package.R \
    --repo-dir "$REPO_DIR" \
    --package-dir "$OUT_DIR" \
    --source-object "$METADATA_OBJECT_FILE"
else
  log_message "Skipping public-package anonymization." --message-type warning
fi

log_message "Compact ScienceDB package prepared at {.file ${OUT_DIR}}" --message-type success
