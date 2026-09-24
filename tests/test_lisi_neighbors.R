#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(lisi))
source("functions/lisi_neighbors.R")

set.seed(2026)
embedding <- matrix(rnorm(400 * 50), nrow = 400)
embedding[201:220, ] <- embedding[1:20, ]
metadata <- data.frame(
  broad = rep(letters[1:4], 100),
  fine = rep(letters[1:8], 50)
)

expected <- lisi::compute_lisi(
  embedding, metadata, names(metadata), perplexity = 30, nn_eps = 0.1
)
observed <- compute_source_lisi(embedding, metadata)
stopifnot(isTRUE(all.equal(expected, observed, tolerance = 1e-12)))

hnsw <- compute_source_lisi_hnsw(embedding, metadata, metadata$broad)
stopifnot(max(abs(as.matrix(hnsw) - as.matrix(expected))) < 0.01)
cat("LISI neighbor checks passed\n")
