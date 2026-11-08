source("functions/utils.R")
set.seed(42)
reference <- matrix(rnorm(400), 80, dimnames = list(paste0("r", 1:80), paste0("PC", 1:5)))
query <- matrix(rnorm(45), 9, dimnames = list(paste0("q", 1:9), colnames(reference)))
coordinates <- reference[, 1:2]
labels <- rep(c("A", "B", "C", "D"), each = 20)
set.seed(100)
result <- reference_knn_projection(query, reference, coordinates, labels, k = 7L)
neighbors <- Seurat::FindNeighbors(
  object = reference, query = query, k.param = 7L,
  nn.method = "annoy", annoy.metric = "cosine", return.neighbor = TRUE,
  index = result$Index, verbose = FALSE
)
idx <- neighbors@nn.idx
expected <- t(vapply(seq_len(nrow(query)), function(i) {
  colMeans(coordinates[idx[i, ], , drop = FALSE])
}, numeric(2)))
stopifnot(isTRUE(all.equal(unname(result$Projection), unname(expected))))
votes <- scop:::knn_vote_labels(matrix(labels[idx], nrow = nrow(idx)), levels = unique(labels))
stopifnot(
  identical(result$Candidate_Prediction, unname(votes$best)),
  all(is.na(result$Prediction)), all(result$Mapping_Status == "scope_not_confirmed"),
  isTRUE(all.equal(result$Support, unname(apply(votes$probability, 1L, max)))),
  all(result$Support >= 0 & result$Support <= 1), all(result$Margin >= 0)
)
second <- reference_knn_projection(query, reference, coordinates, labels, index = result$Index, k = 7L)
stopifnot(identical(second$Candidate_Prediction, result$Candidate_Prediction), identical(second$Projection, result$Projection))
candidate <- c("Inhibitory neurons", "Neural progenitors", "Microglia", "Fibroblasts", "Excitatory neurons")
restricted <- restrict_reference_labels(candidate, c(FALSE, TRUE, TRUE, NA, TRUE))
stopifnot(
  identical(restricted$Prediction, c(NA_character_, NA_character_, "Immune", NA_character_, "Neuron")),
  identical(restricted$Mapping_Status, c("outside_demonstrated_scope", "unsupported_reference_assignment",
    "within_demonstrated_scope", "scope_not_confirmed", "within_demonstrated_scope")),
  inherits(try(restrict_reference_labels(candidate, TRUE), silent = TRUE), "try-error"),
  inherits(try(restrict_reference_labels(candidate, rep("TRUE", 5)), silent = TRUE), "try-error")
)
cat("PASS: projection, candidate votes, default withholding and explicit scope restriction\n")
