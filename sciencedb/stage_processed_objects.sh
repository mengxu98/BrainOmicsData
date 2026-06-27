#!/usr/bin/env bash
set -euo pipefail

INTEGRATION_DIR="${INTEGRATION_DIR:-../../data/BrainOmicsData/integration}"
OUT_DIR="${OUT_DIR:-../../data/ScienceDB/HumanBrain_sc_snRNAseq_age_resource}"
OBJECT_FILE="${OBJECT_FILE:-$INTEGRATION_DIR/objects_celltypes.rds}"
PLOT_OBJECT_FILE="${PLOT_OBJECT_FILE:-$INTEGRATION_DIR/objects_celltype_plot.rds}"
OBJECTS_DIR="$OUT_DIR/04_processed_objects"

mkdir -p "$OBJECTS_DIR"

stage_file() {
  local source_file="$1"
  local target_name="$2"
  local description="$3"
  local target_file="$OBJECTS_DIR/$target_name"

  if [[ ! -f "$source_file" ]]; then
    printf 'Skipping missing processed object: %s\n' "$source_file" >&2
    return 0
  fi

  if [[ -f "$target_file" ]]; then
    rm -f "$target_file"
  fi

  if ln "$source_file" "$target_file" 2>/dev/null; then
    link_method="hardlink"
  else
    cp -p "$source_file" "$target_file"
    link_method="copy"
  fi

  cat > "$OBJECTS_DIR/${target_name}.note.txt" <<EOF
File: $target_name
Source: $source_file
Staging method: $link_method
Description: $description
EOF

  printf 'Staged processed object: %s -> %s (%s)\n' "$source_file" "$target_file" "$link_method"
}

stage_file \
  "$OBJECT_FILE" \
  "human_brain_integrated_full_seurat.rds" \
  "Full integrated Seurat object containing expression assays, harmonized metadata, reductions, clusters and cell-type annotations."

stage_file \
  "$PLOT_OBJECT_FILE" \
  "human_brain_marker_plot_seurat.rds" \
  "Seurat object used for marker-gene plotting and lightweight visualization workflows."
