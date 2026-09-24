#!/usr/bin/env Rscript
# Supplementary Figure S1: reported age intervals on the RPCA embedding.
status <- system2("Rscript", c("--vanilla", "functions/supplementary_umap.R", "age"))
if (status != 0L) stop("Supplementary Figure S1 source generation failed")
source("functions/export_png.R")
export_pdf_png("figures/figS1.pdf", "figures/figS1.png", 6.69)
