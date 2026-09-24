#!/usr/bin/env Rscript
# Supplementary Figure S4: all 33 markers in the adopted formal list.
root <- Sys.getenv("BRAINOMICS_ANALYSIS_DIR", "results/analysis_run")
formal <- file.path(root, "07_downstream/revision_20260918/19_all_celltype_marker_review_20260918")
out <- file.path(root, "07_downstream/revision_20260918/21_figS4_compact_20260923")
if (!file.exists(file.path(formal, "marker_order.tsv"))) stop("Formal marker order is missing")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
if (!file.copy(file.path(formal, "marker_order.tsv"), file.path(out, "marker_order.tsv"), overwrite=TRUE))
  stop("Could not stage formal marker order")
# SCOP_SOURCE_PATH is optional and explicit; environment/scop-plotting.lock.tsv
# identifies the source used for this figure. Never pick a sibling checkout.
status <- system2("Rscript", c("--vanilla", "functions/figS4_source.R", root, out))
if (status != 0L) stop("Supplementary Figure S4 source generation failed")
if (!file.copy(file.path(out, "figS4.pdf"), "figures/figS4.pdf", overwrite=TRUE))
  stop("Could not copy Supplementary Figure S4 PDF")
source("functions/export_png.R")
# 425 dpi at the 240 mm source canvas gives 600 dpi at 170 mm placement.
export_pdf_png("figures/figS4.pdf", "figures/figS4.png", 170/25.4, min_dpi=ceiling(600*170/240))
