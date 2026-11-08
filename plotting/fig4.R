#!/usr/bin/env Rscript
# EGAD transfer, source correspondence and donor consistency.
source("functions/fig4_source.R")
figure_pdf <- file.path(out, "fig4.pdf")
stopifnot(file.exists(figure_pdf))
source("functions/export_png.R")
export_pdf_png(figure_pdf, file.path(out, "fig4.png"), 9.45)
