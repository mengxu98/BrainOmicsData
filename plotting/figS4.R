#!/usr/bin/env Rscript
# Supplementary Figure S4: all 33 markers in the adopted list.
root <- Sys.getenv("BRAINOMICS_ANALYSIS_DIR", "results/analysis_run")
marker_file <- Sys.getenv(
  "BRAINOMICS_MARKER_ORDER_FILE",
  file.path("results", "annotation", "display_markers.tsv")
)
out <- file.path("results", "figure_sources", "figS4")
if (!file.exists(marker_file)) stop("Marker order is missing: ", marker_file)
dir.create(out, recursive=TRUE, showWarnings=FALSE)
if (!file.copy(marker_file, file.path(out, "marker_order.tsv"), overwrite=TRUE))
  stop("Could not stage marker order")
# All figures use the installed package versions in environment/r-packages.lock.tsv.
status <- system2(
  Sys.getenv("BRAINOMICS_RSCRIPT", "Rscript"),
  c("--vanilla", "functions/figS4_source.R", shQuote(root), shQuote(out))
)
if (status != 0L) stop("Supplementary Figure S4 source generation failed")
if (!file.copy(file.path(out, "figS4.pdf"), "figures/figS4.pdf", overwrite=TRUE))
  stop("Could not copy Supplementary Figure S4 PDF")
source("functions/export_png.R")
# Use the measured PDF width to guarantee 600 dpi at 170 mm placement.
export_pdf_png("figures/figS4.pdf", "figures/figS4.png", 170/25.4,
               min_dpi=ceiling(600*170/240), min_effective_dpi=600L)
