#!/usr/bin/env Rscript
source("functions/utils.R")
a <- read_celltype_assignments()
stopifnot(nrow(unique(a[c("Cluster", "CellType")])) == 75L)
index <- c(30L, 1L, 500L)
stopifnot(identical(read_celltype_assignments(a$Cells[index]), a[index, , drop = FALSE]))
cat("PASS: cell counts, cluster mapping, single 12-class annotation and requested cell order\n")
