#!/usr/bin/env Rscript
# EGAD transfer, source correspondence and donor consistency.
previous_transfer_mode <- Sys.getenv("EGAD_TRANSFER_CURRENT", unset = "")
Sys.setenv(EGAD_TRANSFER_CURRENT = "1")
source("functions/egad_reference_review.R")
stopifnot(file.exists("figures/fig4.pdf"))
Sys.setenv(EGAD_TRANSFER_CURRENT = previous_transfer_mode)
source("functions/export_png.R")
export_pdf_png("figures/fig4.pdf", "figures/fig4.png", 9.45)
