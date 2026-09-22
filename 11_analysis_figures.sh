#!/usr/bin/env bash
# Refresh the annotation and integration tables, then draw every manuscript
# figure from the frozen analysis run.
#
# Phase 1 rewrites the annotation, resource-comparison, reuse-case and
# source-concordance tables and needs the frozen integration objects. Set
# BRAINOMICS_SKIP_TABLE_REFRESH=1 to draw the figures from the frozen tables
# only (for example when the integration objects are kept on another host).
#
# Requirements: 50 GB memory; no GPU.
set -euo pipefail

BRAINOMICS_STAGE=10_annotation_figures
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/functions/pipeline_lib.sh"
source "functions/utils.sh"

cd "$BRAINOMICS_REPO_ROOT"
brainomics_require_executor

if [ "$BRAINOMICS_EXECUTOR" = "slurm" ]; then
  brainomics_die "annotation and figure refresh runs on the storage host; use BRAINOMICS_EXECUTOR=local"
fi

check_command Rscript

if [ "${BRAINOMICS_SKIP_TABLE_REFRESH:-0}" != "1" ]; then
  Rscript annotation/export_annotation_inputs.R
  Rscript integration/export_resource_comparison.R
  Rscript analysis/reuse_case_dact1_expression_context.R
  Rscript integration/evaluate_source_concordance.R
fi

# Figures 1-4 and their panels, read from the frozen analysis run.
Rscript --vanilla plotting/fig1.R
Rscript --vanilla plotting/fig2_umap_panels.R
Rscript --vanilla plotting/fig2.R
Rscript --vanilla plotting/cluster_marker_balance.R
Rscript --vanilla plotting/reuse_gene_panel.R
Rscript --vanilla plotting/fig3_annotation_panels.R
Rscript --vanilla plotting/egad_reference_review.R
Rscript --vanilla plotting/assemble_figures.R

brainomics_log "annotation-dependent summaries and figures refreshed"
