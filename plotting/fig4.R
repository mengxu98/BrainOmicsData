#!/usr/bin/env Rscript
# EGAD transfer, source correspondence and donor consistency.
source("functions/fig4_source.R")
stopifnot(file.exists("figures/fig4.pdf"))
source("functions/export_png.R")
export_pdf_png("figures/fig4.pdf", "figures/fig4.png", 9.45)
