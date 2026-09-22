suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})
set.seed(5)
counts <- Matrix(matrix(rpois(120, 4), nrow = 12), sparse = TRUE)
dimnames(counts) <- list(paste0("g", 1:12), paste0("c", 1:10))
selected <- counts[1:4, , drop = FALSE]
totals <- Matrix::colSums(counts)
observed <- selected
observed@x <- log1p(observed@x * rep.int(10000 / totals, diff(observed@p)))
expected <- NormalizeData(counts,
  normalization.method = "LogNormalize", scale.factor = 10000,
  verbose = FALSE
)[1:4, , drop = FALSE]
stopifnot(
  isTRUE(all.equal(as.matrix(observed), as.matrix(expected), tolerance = 1e-12)),
  all(totals >= Matrix::colSums(selected))
)
cat("PASS: query feature normalization agrees with full-RNA LogNormalize\n")
