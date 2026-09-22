suppressPackageStartupMessages({
  library(data.table)
})

source("functions/data_paths.R")
source("functions/dataset_metadata.R")
source("functions/processed_object.R")
source("functions/integration.R")


source("functions/lisi_neighbors.R")

if ("--check-query-partition" %in% commandArgs(trailingOnly = TRUE)) {
  .libPaths(c(brainomics_data_path("environment", "lisi-double", "library"), .libPaths()))
  set.seed(2026)
  x <- matrix(rnorm(400 * 50), nrow = 400)
  x[201:220, ] <- x[1:20, ] # Include tied coordinates, not just generic random points.
  m <- data.frame(broad = rep(letters[1:4], 100), fine = rep(letters[1:8], 50))
  expected <- lisi::compute_lisi(x, m, names(m), perplexity = 30, nn_eps = 0.1)
  observed <- compute_source_lisi(x, m)
  stopifnot(isTRUE(all.equal(expected, observed, tolerance = 1e-12)))
  fast <- compute_source_lisi_hnsw(x, m, m$broad)
  stopifnot(max(abs(as.matrix(fast) - as.matrix(expected))) < 0.01)
  message("HNSW full scoring agrees with original LISI on the tied-coordinate fixture")
  message("Partitioned queries match the original LISI function, including tied coordinates")
  quit(save = "no")
}

method <- trimws(Sys.getenv("BRAINOMICS_ANNOTATION_METHOD", unset = ""))
if (!method %in% method_levels) {
  stop(
    "BRAINOMICS_ANNOTATION_METHOD must be one of: ",
    paste(method_levels, collapse = ", ")
  )
}
overwrite <- tolower(Sys.getenv(
  "BRAINOMICS_EVALUATION_OVERWRITE",
  unset = "false"
)) %in% c("true", "t", "1")
perplexity <- as.numeric(Sys.getenv(
  "BRAINOMICS_LISI_PERPLEXITY",
  unset = "30"
))
nn_eps <- as.numeric(Sys.getenv(
  "BRAINOMICS_LISI_NN_EPS",
  unset = "0.1"
))
if (!is.finite(perplexity) || perplexity <= 0L ||
  !is.finite(nn_eps) || nn_eps < 0) {
  stop("Invalid LISI perplexity or nearest-neighbor error bound")
}
task_library <- brainomics_data_path("environment", "lisi-double", "library")
if (dir.exists(task_library)) .libPaths(c(task_library, .libPaths()))
if (!requireNamespace("lisi", quietly = TRUE)) {
  stop("The pinned lisi package is required")
}
description <- packageDescription("lisi")
if (!identical(description$BrainOmicsNumericPrecision, "double")) {
  stop("Use the reviewed double-precision LISI build for source-label evaluation")
}

integration_dir <- brainomics_data_path("integration_25")
annotation_dir <- file.path(integration_dir, "annotation")
reduction_dir <- file.path(annotation_dir, "reductions")
checkpoint_dir <- file.path(annotation_dir, "lisi_checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

method_key <- tolower(method)
embedding_file <- file.path(
  reduction_dir,
  paste0(method_key, "_latent.rds")
)
label_file <- file.path(annotation_dir, "source_labels_harmonized.rds")
reduction_audit_file <- file.path(
  annotation_dir,
  "annotation_reduction_audit.tsv"
)
label_contract_file <- file.path(
  annotation_dir,
  "source_label_cell_contract.tsv"
)
output_file <- file.path(
  checkpoint_dir,
  paste0(method_key, "_source_label_lisi.rds")
)
required_files <- c(
  embedding_file, label_file, reduction_audit_file, label_contract_file
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Missing source-label LISI input: ", paste(missing_files, collapse = ", "))
}

labels <- readRDS(label_file)
required_label_fields <- c(
  "Cells", "Dataset", "Donor_ID", "Harmonized_CellClass",
  "Harmonized_CellType", "Mapping_Status",
  "Source_Label_Evaluation_Eligible", "Row_Order"
)
missing_label_fields <- setdiff(required_label_fields, names(labels))
if (length(missing_label_fields) > 0L) {
  stop(
    "Per-cell source-label mapping is missing: ",
    paste(missing_label_fields, collapse = ", ")
  )
}
if (anyDuplicated(labels$Cells) ||
  !identical(labels$Row_Order, seq_len(nrow(labels)))) {
  stop("Per-cell source-label mapping order is invalid")
}

reduction_audit <- read.delim(
  reduction_audit_file,
  sep = "\t",
  quote = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
audit_row <- reduction_audit[reduction_audit$Method == method, , drop = FALSE]
if (nrow(audit_row) != 1L) stop("Missing method in reduction audit")
label_contract <- read.delim(
  label_contract_file,
  sep = "\t",
  quote = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
contract_values <- stats::setNames(
  as.character(label_contract$Value),
  label_contract$Field
)
expected_eligible_cells <- as.integer(
  contract_values[["Mapped_Evaluation_Cells"]]
)
if (is.na(expected_eligible_cells) || expected_eligible_cells <= 0L) {
  stop("Invalid mapped-evaluation-cell count in source-label contract")
}

thisutils::log_message("[source-label-lisi] ", paste(
  "Loading",
  method,
  audit_row$Dimensions,
  "dimensional latent representation"
))
embedding <- readRDS(embedding_file)
if (!is.matrix(embedding) ||
  nrow(embedding) != nrow(labels) ||
  ncol(embedding) != audit_row$Dimensions ||
  !identical(rownames(embedding), as.character(labels$Cells)) ||
  any(!is.finite(embedding))) {
  stop(method, " latent representation failed dimension, order or value checks")
}

eligible <- labels$Source_Label_Evaluation_Eligible %in% TRUE &
  labels$Mapping_Status == "mapped_for_common_comparison" &
  !is.na(labels$Harmonized_CellClass) &
  !is.na(labels$Harmonized_CellType)
if (!any(eligible) || sum(eligible) != expected_eligible_cells) {
  stop("Common source-label evaluation cohort differs from its contract")
}
evaluation_embedding <- embedding[eligible, , drop = FALSE]
evaluation_metadata <- data.frame(
  Source_CellClass = labels$Harmonized_CellClass[eligible],
  Source_CellType = labels$Harmonized_CellType[eligible],
  stringsAsFactors = FALSE,
  row.names = labels$Cells[eligible]
)
if (!identical(rownames(evaluation_embedding), rownames(evaluation_metadata))) {
  stop("Source-label evaluation embedding and metadata orders differ")
}
rm(embedding)
gc()

if (file.exists(output_file) && !overwrite) {
  thisutils::log_message("[source-label-lisi] ", paste("Validating existing", method, "source-label LISI checkpoint"))
  result <- readRDS(output_file)
} else {
  thisutils::log_message("[source-label-lisi] ", paste(
    "Computing source CellClass and CellType cLISI for",
    method,
    "on",
    format(nrow(evaluation_embedding), big.mark = ","),
    "eligible cells"
  ))
  neighbor_method <- Sys.getenv("BRAINOMICS_LISI_NEIGHBORS", "rann")
  if (!neighbor_method %in% c("rann", "hnsw")) stop("Unknown neighbor method")
  result <- if (neighbor_method == "hnsw") {
    compute_source_lisi_hnsw(evaluation_embedding, evaluation_metadata,
      labels$Dataset[eligible], perplexity)
  } else {
    compute_source_lisi(evaluation_embedding, evaluation_metadata, perplexity, nn_eps)
  }
  attr(result, "neighbor_method") <- neighbor_method
  attr(result, "neighbor_parameters") <- if (neighbor_method == "hnsw") {
    "RcppHNSW; Euclidean; M=24; ef_construction=200; ef_search=400; threads=6; seed=2026; exact stratified check"
  } else paste("RANN; eps=", nn_eps)
  attr(result, "neighbor_package_version") <- as.character(packageVersion(
    if (neighbor_method == "hnsw") "RcppHNSW" else "RANN"))
  result$Source_CellClass <- canonicalize_lisi_lower_bound(
    result$Source_CellClass
  )
  result$Source_CellType <- canonicalize_lisi_lower_bound(
    result$Source_CellType
  )
  rownames(result) <- rownames(evaluation_metadata)
  attr(result, "method") <- method
  attr(result, "space") <- "method-specific 50-dimensional latent"
  attr(result, "labels") <- c("Source_CellClass", "Source_CellType")
  attr(result, "perplexity") <- perplexity
  attr(result, "nn_eps") <- if (neighbor_method == "rann") nn_eps else NA_real_
  attr(result, "eligible_cells") <- nrow(evaluation_metadata)
  attr(result, "source_labels") <- evaluation_metadata
  attr(result, "numeric_precision") <- "double"
  attr(result, "lisi_source_commit") <- unname(
    description$RemoteSha
  )
  temporary_output <- paste0(output_file, ".tmp.", Sys.getpid())
  saveRDS(result, temporary_output, compress = TRUE)
  if (!file.rename(temporary_output, output_file)) {
    unlink(temporary_output)
    stop("Could not publish ", method, " source-label LISI checkpoint")
  }
}

if (!identical(
  names(result),
  c("Source_CellClass", "Source_CellType")
) || nrow(result) != nrow(evaluation_metadata) ||
  !identical(attr(result, "method", exact = TRUE), method) ||
  !identical(rownames(result), rownames(evaluation_metadata)) ||
  !identical(attr(result, "source_labels"), evaluation_metadata) ||
  any(!is.finite(as.matrix(result))) ||
  any(as.matrix(result) < 1)) {
  stop(method, " source-label LISI checkpoint failed validation")
}
thisutils::log_message("[source-label-lisi] ", paste("Completed", method, "source-label cLISI checkpoint"))
