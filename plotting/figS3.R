#!/usr/bin/env Rscript
# Supplementary Figure S3: formal cell-type location UMAP.
# Use --redraw only where the complete analysis inputs are available.
# The exact recovered original is analysis/revision_sources/plot_celltype_locations.R;
# its original output is byte-identical to the retained formal PDF.
if ("--redraw" %in% commandArgs(trailingOnly = TRUE)) {
  root <- Sys.getenv("BRAINOMICS_ANALYSIS_DIR", "results/analysis_run")
  input <- file.path(root, "07_downstream/revision_20260918/13_heatmap_ordered_umaps_20260918")
  out <- tempfile("figS3-redraw-")
  dir.create(out)
  for (name in c("cluster_annotation.tsv", "celltype_order.tsv")) {
    if (!file.copy(file.path(input, name), file.path(out, name))) stop("Missing S3 input: ", name)
  }
  status <- system2("Rscript", c("--vanilla", "functions/figS3_source.R", shQuote(root), shQuote(out)))
  if (status != 0L) stop("S3 redraw failed; formal figure was not overwritten")
  if (!file.copy(file.path(out, "celltype_location_umap.pdf"), "figures/figS3.pdf", overwrite = TRUE))
    stop("Could not save redrawn S3")
} else {
  message("S3: exporting PNG from retained PDF; use --redraw for the recovered drawing source")
}
source("functions/export_png.R")
export_pdf_png("figures/figS3.pdf",
               "figures/figS3.png", 6.69)
