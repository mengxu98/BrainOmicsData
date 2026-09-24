#!/usr/bin/env Rscript
# Supplementary Figure S3: cell-type locations on the RPCA embedding.
if ("--redraw" %in% commandArgs(trailingOnly = TRUE)) {
  root <- Sys.getenv("BRAINOMICS_ANALYSIS_DIR", "results/analysis_run")
  input <- Sys.getenv(
    "BRAINOMICS_ANNOTATION_TABLE",
    file.path("results", "annotation", "cluster_annotation.tsv")
  )
  out <- file.path("results", "figure_sources", "figS3")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(input, file.path(out, "cluster_annotation.tsv"), overwrite = TRUE))
    stop("Missing S3 annotation input: ", input)
  status <- system2(Sys.getenv("BRAINOMICS_RSCRIPT", "Rscript"),
                    c("--vanilla", "functions/figS3_source.R", shQuote(root), shQuote(out)))
  if (status != 0L) stop("S3 drawing failed")
  if (!file.copy(file.path(out, "celltype_location_umap.pdf"), "figures/figS3.pdf", overwrite = TRUE))
    stop("Could not save Figure S3")
}
source("functions/export_png.R")
export_pdf_png("figures/figS3.pdf",
               "figures/figS3.png", 6.69)
