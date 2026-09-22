#!/usr/bin/env Rscript
# EGAD transfer, source correspondence and donor consistency.
previous_transfer_mode <- Sys.getenv("EGAD_TRANSFER_CURRENT", unset = "")
Sys.setenv(EGAD_TRANSFER_CURRENT = "1")
source("plotting/egad_reference_review.R")
stopifnot(file.copy("figures/egad_transfer_review.pdf", "figures/fig4.pdf", overwrite = TRUE))
Sys.setenv(EGAD_TRANSFER_CURRENT = previous_transfer_mode)
