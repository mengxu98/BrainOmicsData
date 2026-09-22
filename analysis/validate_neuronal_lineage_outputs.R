#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

source("functions/data_paths.R")
source("functions/processed_object.R")
source("functions/integration.R")

lineage_id <- Sys.getenv("BRAINOMICS_LINEAGE_ID", unset = "")
expected_celltype <- Sys.getenv(
  "BRAINOMICS_LINEAGE_CELLTYPE",
  unset = ""
)
expected_cells <- as.numeric(Sys.getenv(
  "BRAINOMICS_LINEAGE_EXPECTED_CELLS",
  unset = "0"
))
expected_datasets <- as.integer(Sys.getenv(
  "BRAINOMICS_LINEAGE_EXPECTED_DATASETS",
  unset = "25"
))
expected_role <- Sys.getenv(
  "BRAINOMICS_LINEAGE_ANALYSIS_ROLE",
  unset = ""
)
lineage_dir <- brainomics_lineage_path(lineage_id, expected_role)
if (any(!nzchar(c(
  lineage_id, expected_celltype, expected_role
))) || !dir.exists(lineage_dir) || !is.finite(expected_cells) ||
  expected_cells <= 0 || !is.finite(expected_datasets) ||
  expected_datasets <= 0L) {
  stop("Complete lineage validation environment is required")
}

paths <- list(
  exit = file.path(lineage_dir, "run", "pipeline.exit"),
  cells = file.path(lineage_dir, "lineage_cells.tsv.gz"),
  datasets = file.path(lineage_dir, "lineage_dataset_counts.tsv"),
  parameters = file.path(lineage_dir, "lineage_integration_parameters.tsv"),
  integration_manifest = file.path(
    lineage_dir, "lineage_output_manifest.tsv"
  ),
  raw_pca = file.path(lineage_dir, "lineage_raw_pca_latent.rds"),
  rpca = file.path(lineage_dir, "lineage_rpca_latent.rds"),
  umap = file.path(lineage_dir, "lineage_rpca_umap.tsv.gz"),
  assignments = file.path(
    lineage_dir, "clustering", "lineage_cluster_assignments.tsv.gz"
  ),
  cluster_summary = file.path(
    lineage_dir, "clustering", "lineage_cluster_summary.tsv"
  ),
  cluster_stability = file.path(
    lineage_dir, "clustering", "lineage_cluster_stability.tsv"
  ),
  cluster_parameters = file.path(
    lineage_dir, "clustering", "lineage_clustering_parameters.tsv"
  ),
  cluster_manifest = file.path(
    lineage_dir, "clustering", "lineage_clustering_output_manifest.tsv"
  )
)
missing <- unlist(paths)[!file.exists(unlist(paths))]
if (length(missing) > 0L) {
  stop("Lineage outputs are incomplete: ", paste(missing, collapse = ", "))
}
if (!identical(trimws(readLines(paths$exit, warn = FALSE)), "0")) {
  stop("Lineage pipeline exit status is not zero")
}

validate_manifest <- function(path) {
  # Read byte sizes as text so validation does not depend on the optional
  # bit64 package. Without bit64, coercing data.table's integer64 storage with
  # base as.numeric() returns the underlying bit pattern instead of the size.
  manifest <- fread(
    path,
    sep = "\t",
    colClasses = c(
      File = "character", Size_Bytes = "character"
    )
  )
  required <- c("File", "Size_Bytes")
  expected_size <- suppressWarnings(as.numeric(manifest$Size_Bytes))
  if (any(!required %in% names(manifest)) || nrow(manifest) == 0L ||
    anyNA(manifest[, ..required]) || anyDuplicated(manifest$File) ||
    any(!is.finite(expected_size)) || any(expected_size < 0) ||
    any(expected_size != floor(expected_size)) ||
    any(!file.exists(manifest$File))) {
    stop("Invalid output manifest: ", path)
  }
  actual_size <- as.numeric(file.info(manifest$File)$size)
  if (!identical(actual_size, expected_size)) {
    stop("Size mismatch in manifest: ", path)
  }
  manifest
}
integration_manifest <- validate_manifest(paths$integration_manifest)
cluster_manifest <- validate_manifest(paths$cluster_manifest)

cells <- fread(paths$cells, sep = "\t")
required_cell_columns <- c("Cells", "Global_Cluster", "Main_CellType")
if (nrow(cells) != expected_cells ||
  any(!required_cell_columns %in% names(cells)) ||
  anyNA(cells[, ..required_cell_columns]) || anyDuplicated(cells$Cells) ||
  !identical(unique(as.character(cells$Main_CellType)), expected_celltype)) {
  stop(
    "Lineage cell sidecar failed its complete-cell contract: observed rows=",
    nrow(cells), "; expected rows=", expected_cells,
    "; observed Main_CellType=",
    paste(unique(as.character(cells$Main_CellType)), collapse = ";"),
    "; expected Main_CellType=", expected_celltype,
    "; missing required values=",
    if (all(required_cell_columns %in% names(cells))) {
      sum(is.na(cells[, ..required_cell_columns]))
    } else {
      "columns_missing"
    },
    "; duplicated cells=", anyDuplicated(cells$Cells)
  )
}

raw_pca <- readRDS(paths$raw_pca)
rpca <- readRDS(paths$rpca)
umap <- fread(paths$umap, sep = "\t")
if (!is.matrix(raw_pca) || !is.matrix(rpca) ||
  !identical(dim(raw_pca), c(as.integer(expected_cells), 50L)) ||
  !identical(dim(rpca), c(as.integer(expected_cells), 50L)) ||
  !identical(rownames(raw_pca), cells$Cells) ||
  !identical(rownames(rpca), cells$Cells) ||
  any(!is.finite(raw_pca)) || any(!is.finite(rpca)) ||
  nrow(umap) != expected_cells ||
  !identical(as.character(umap$Cells), cells$Cells) ||
  any(!is.finite(as.matrix(umap[, .(UMAP_1, UMAP_2)])))) {
  stop("Lineage latent-space or UMAP validation failed")
}
rm(raw_pca, rpca, umap)
gc()

dataset_counts <- fread(paths$datasets, sep = "\t")
if (!identical(names(dataset_counts), c("Dataset", "Cells")) ||
  nrow(dataset_counts) != expected_datasets ||
  anyNA(dataset_counts) || anyDuplicated(dataset_counts$Dataset) ||
  sum(as.numeric(dataset_counts$Cells)) != expected_cells ||
  any(as.numeric(dataset_counts$Cells) <= 0)) {
  stop("Lineage dataset coverage is incomplete")
}

parameters <- fread(paths$parameters, sep = "\t")
parameter <- setNames(as.character(parameters$Value), parameters$Parameter)
required_parameter_values <- c(
  lineage_id = lineage_id,
  lineage_celltype = expected_celltype,
  cells = as.character(expected_cells),
  datasets = as.character(expected_datasets),
  dimensions = "1-50",
  batch_model = "Integration_Batch_ID: dataset",
  analysis_role = expected_role,
  global_cluster_preserved = "true; stored as Global_Cluster"
)
if (any(!names(required_parameter_values) %in% names(parameter)) ||
  !identical(
    unname(parameter[names(required_parameter_values)]),
    unname(required_parameter_values)
  )) {
  stop("Lineage integration parameter contract changed")
}

assignments <- fread(paths$assignments, sep = "\t")
required_assignment_columns <- c(
  required_cell_columns, "Lineage_Cluster"
)
seed_columns <- grep("^Seed[0-9]+_Resolution", names(assignments), value = TRUE)
if (nrow(assignments) != expected_cells ||
  any(!required_assignment_columns %in% names(assignments)) ||
  length(seed_columns) < 6L ||
  anyNA(assignments[, c(required_assignment_columns, seed_columns), with = FALSE]) ||
  anyDuplicated(assignments$Cells) ||
  !identical(as.character(assignments$Cells), cells$Cells) ||
  !identical(
    as.character(assignments$Global_Cluster),
    as.character(cells$Global_Cluster)
  ) || !identical(
  as.character(assignments$Main_CellType),
  as.character(cells$Main_CellType)
) || !all(grepl(
  paste0("^", lineage_id, "_[0-9]{3}$"),
  assignments$Lineage_Cluster
))) {
  stop("Lineage clustering assignments failed validation")
}

cluster_summary <- fread(paths$cluster_summary, sep = "\t")
cluster_stability <- fread(paths$cluster_stability, sep = "\t")
cluster_parameters <- fread(paths$cluster_parameters, sep = "\t")
cluster_parameter <- setNames(
  as.character(cluster_parameters$Value), cluster_parameters$Parameter
)
if (nrow(cluster_summary) != length(seed_columns) ||
  anyNA(cluster_summary) || any(cluster_summary$Clusters < 2L) ||
  any(cluster_summary$Smallest_Cluster_Cells <= 0L) ||
  nrow(cluster_stability) < 9L ||
  anyNA(cluster_stability$Adjusted_Rand_Index) ||
  any(cluster_stability$Adjusted_Rand_Index < 0 |
    cluster_stability$Adjusted_Rand_Index > 1) ||
  !identical(cluster_parameter[["cells"]], as.character(expected_cells)) ||
  !identical(cluster_parameter[["dimensions"]], "1-50") ||
  !identical(cluster_parameter[["global_cluster_preserved"]], "true")) {
  stop("Lineage clustering summary or stability contract failed")
}

validation <- data.table(
  Field = c(
    "Lineage_ID", "Main_CellType", "Analysis_Role", "Cells", "Datasets",
    "Raw_PCA_Dimensions", "RPCA_Dimensions", "UMAP_Dimensions",
    "Global_Clusters_Preserved", "Clustering_Assignments",
    "Primary_Lineage_Clusters", "Seed_Stability_Comparisons",
    "Resolution_Comparisons",
    "Cell_Reduction"
  ),
  Value = c(
    lineage_id, expected_celltype, expected_role, expected_cells,
    expected_datasets, 50L, 50L, 2L, TRUE, length(seed_columns),
    uniqueN(assignments$Lineage_Cluster),
    sum(cluster_stability$Comparison_Type == "seed_at_primary_resolution"),
    sum(cluster_stability$Comparison_Type == "resolution_at_primary_seed"),
    "none; every input lineage cell retained exactly once"
  )
)
output <- file.path(lineage_dir, "lineage_validation_summary.tsv")
fwrite(validation, output, sep = "\t", quote = FALSE)
success <- file.path(lineage_dir, "_VALIDATED")
writeLines(
  paste(
    "validated", lineage_id, expected_role, expected_cells, "cells"
  ),
  success
)
message(
  "[validate-neuronal-lineage] validated ", lineage_id, ": ",
  format(expected_cells, big.mark = ","), " cells, ", expected_datasets,
  " datasets, ", uniqueN(assignments$Lineage_Cluster),
  " primary clusters"
)
