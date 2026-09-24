#!/usr/bin/env Rscript
# Supplementary Figure S2: reported brain regions on the RPCA embedding.
rscript <- Sys.getenv("BRAINOMICS_RSCRIPT", file.path(R.home("bin"), "Rscript"))
status <- system2(rscript, c("--vanilla", "functions/supplementary_umap.R", "region"))
if (status != 0L) stop("Supplementary Figure S2 source generation failed")
source("functions/export_png.R")
export_pdf_png("figures/figS2.pdf", "figures/figS2.png", 6.69)
