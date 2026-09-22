# Partition query rows only: every query still searches the entire common cohort.
compute_source_lisi <- function(X, meta_data, perplexity = 30, nn_eps = 0.1, workers = 2L) {
  encoded <- lapply(meta_data, function(x) as.integer(factor(x)) - 1L)
  chunks <- split(seq_len(nrow(X)), cut(seq_len(nrow(X)), breaks = workers, labels = FALSE))
  pieces <- parallel::mclapply(chunks, function(ii) {
    neighbors <- RANN::nn2(data = X, query = X[ii, , drop = FALSE],
      k = perplexity * 3L, eps = nn_eps)
    distances <- t(neighbors$nn.dists[, -1L, drop = FALSE])
    indices <- t(neighbors$nn.idx[, -1L, drop = FALSE]) - 1L
    as.data.frame(lapply(encoded, function(labels) {
      1 / lisi:::compute_simpson_index(distances, indices, labels,
        length(unique(labels)), perplexity)
    }))
  }, mc.cores = workers)
  if (any(vapply(pieces, inherits, logical(1), "try-error"))) stop("A LISI query worker failed")
  result <- do.call(rbind, pieces)
  rownames(result) <- rownames(meta_data)
  result
}

# Fixed HNSW approximation; exact full-cohort query checks precede full scoring.
compute_source_lisi_hnsw <- function(X, meta_data, datasets, perplexity = 30) {
  set.seed(2026)
  check_rows <- unlist(lapply(split(seq_len(nrow(X)), datasets), function(ii) {
    sample(ii, min(12L, length(ii)))
  }), use.names = FALSE)
  encoded <- lapply(meta_data, function(x) as.integer(factor(x)) - 1L)
  score_neighbors <- function(distance, index) {
    as.data.frame(lapply(encoded, function(label) {
      1 / lisi:::compute_simpson_index(t(distance[, -1L, drop = FALSE]),
        t(index[, -1L, drop = FALSE]) - 1L, label, length(unique(label)), perplexity)
    }))
  }
  message("Building HNSW: M=24, ef_construction=200, ef_search=400, threads=6, seed=2026")
  index <- RcppHNSW::hnsw_build(X, distance = "euclidean", M = 24, ef = 200,
    n_threads = 6, random_seed = 2026)
  k <- as.integer(3 * perplexity)
  approx <- RcppHNSW::hnsw_search(X[check_rows, , drop = FALSE], index,
    k = k, ef = 400, n_threads = 6)
  message("Checking ", length(check_rows), " stratified query cells against exact full-cohort RANN")
  exact <- RANN::nn2(X, X[check_rows, , drop = FALSE], k = k, eps = 0)
  recall <- vapply(seq_along(check_rows), function(i) {
    length(intersect(approx$idx[i, ], exact$nn.idx[i, ])) / k
  }, numeric(1))
  error <- abs(as.matrix(score_neighbors(approx$dist, approx$idx)) -
    as.matrix(score_neighbors(exact$nn.dists, exact$nn.idx)))
  audit <- data.frame(Check_Cells = length(check_rows), Neighbor_Recall = mean(recall),
    Mean_Absolute_cLISI_Error = mean(error), P95_Absolute_cLISI_Error = unname(quantile(error, 0.95)),
    Max_Absolute_cLISI_Error = max(error))
  print(audit)
  if (mean(recall) < 0.95 || any(colMeans(error) > 0.005) ||
      any(apply(error, 2, quantile, 0.95) > 0.05)) {
    stop("HNSW did not pass the fixed exact-search accuracy gate; no full scores published")
  }
  result <- matrix(NA_real_, nrow(X), length(encoded), dimnames = list(rownames(meta_data), names(encoded)))
  for (start in seq.int(1L, nrow(X), by = 100000L)) {
    rows <- start:min(start + 99999L, nrow(X))
    nn <- RcppHNSW::hnsw_search(X[rows, , drop = FALSE], index,
      k = k, ef = 400, n_threads = 6)
    result[rows, ] <- as.matrix(score_neighbors(nn$dist, nn$idx))
    message("Scored ", max(rows), "/", nrow(X), " cells")
  }
  result <- as.data.frame(result)
  attr(result, "neighbor_accuracy") <- audit
  result
}

