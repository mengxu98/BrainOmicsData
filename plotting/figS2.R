#!/usr/bin/env Rscript
# Supplementary Figure S2: reported brain regions on the RPCA embedding.
status <- system2("Rscript", c("--vanilla", "functions/supplementary_umap.R", "region"))
if (status != 0L) stop("Supplementary Figure S2 source generation failed")
source("functions/export_png.R")
export_pdf_png("figures/figS2.pdf", "figures/figS2.png", 6.69)
