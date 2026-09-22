#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(SeuratObject)
})

source("functions/data_paths.R")
source("functions/processed_object.R")
source("functions/integration.R")

# FindNeighbors uses future-aware apply code even under a sequential plan. The
# complete ExN 50-dimensional matrix exceeds future's 500-MiB default export
# guard, so set an explicit finite ceiling without creating parallel copies.
future::plan(future::sequential)
options(future.globals.maxSize = 4 * 1024^3)


lineage_id <- Sys.getenv("BRAINOMICS_LINEAGE_ID", unset = "")
analysis_role <- Sys.getenv("BRAINOMICS_LINEAGE_ANALYSIS_ROLE", unset = "primary")
lineage_dir <- brainomics_lineage_path(lineage_id, analysis_role)
if (!dir.exists(lineage_dir)) {
  stop("Completed lineage directory is required: ", lineage_dir)
}
embedding_file <- file.path(lineage_dir, "lineage_rpca_latent.rds")
cell_file <- file.path(lineage_dir, "lineage_cells.tsv.gz")
required_files <- c(embedding_file, cell_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Missing lineage-clustering input: ", paste(missing_files, collapse = ", "))
}

seed <- as.integer(Sys.getenv("BRAINOMICS_LINEAGE_CLUSTER_SEED", unset = "20260730"))
dims <- seq_len(as.integer(Sys.getenv(
  "BRAINOMICS_LINEAGE_CLUSTER_DIMENSIONS",
  unset = "50"
)))
k_neighbors <- as.integer(Sys.getenv(
  "BRAINOMICS_LINEAGE_CLUSTER_NEIGHBORS",
  unset = "20"
))
parse_numeric_vector <- function(value, field) {
  parsed <- suppressWarnings(as.numeric(strsplit(value, ",", fixed = TRUE)[[1L]]))
  if (length(parsed) == 0L || any(!is.finite(parsed)) || any(parsed <= 0)) {
    stop(field, " must contain comma-separated positive numbers")
  }
  unique(parsed)
}
resolutions <- parse_numeric_vector(
  Sys.getenv(
    "BRAINOMICS_LINEAGE_CLUSTER_RESOLUTIONS",
    unset = "0.5,1.0,1.5,2.0"
  ),
  "BRAINOMICS_LINEAGE_CLUSTER_RESOLUTIONS"
)
primary_resolution <- as.numeric(Sys.getenv(
  "BRAINOMICS_LINEAGE_CLUSTER_PRIMARY_RESOLUTION",
  unset = "1.0"
))
stability_seeds <- unique(c(seed, seed + 1L, seed + 2L))
if (!is.finite(seed) || seed < 1L || length(dims) < 10L ||
  max(dims) > 100L || !is.finite(k_neighbors) || k_neighbors < 5L ||
  length(resolutions) < 2L ||
  !any(abs(resolutions - primary_resolution) < .Machine$double.eps^0.5)) {
  stop("Invalid lineage clustering seed, dimensions, neighbors or resolution")
}

cluster_dir <- file.path(lineage_dir, "clustering")
checkpoint_dir <- file.path(cluster_dir, "checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
neighbor_file <- file.path(checkpoint_dir, "lineage_rpca_knn.rds")
graph_file <- file.path(checkpoint_dir, "lineage_rpca_graph.rds")
assignment_file <- file.path(cluster_dir, "lineage_cluster_assignments.tsv.gz")
summary_file <- file.path(cluster_dir, "lineage_cluster_summary.tsv")
stability_file <- file.path(cluster_dir, "lineage_cluster_stability.tsv")
parameter_file <- file.path(cluster_dir, "lineage_clustering_parameters.tsv")
input_manifest_file <- file.path(cluster_dir, "lineage_clustering_input_manifest.tsv")
output_manifest_file <- file.path(cluster_dir, "lineage_clustering_output_manifest.tsv")
material_outputs <- c(
  assignment_file, summary_file, stability_file, parameter_file,
  input_manifest_file, output_manifest_file
)
if (any(file.exists(material_outputs))) {
  stop("Lineage clustering refuses to overwrite completed outputs")
}

cells <- fread(cell_file, sep = "\t")
required_cell_columns <- c("Cells", "Global_Cluster", "Main_CellType")
if (!all(required_cell_columns %in% names(cells)) ||
  anyNA(cells[, ..required_cell_columns]) || anyDuplicated(cells$Cells)) {
  stop("Lineage cell sidecar failed its schema contract")
}
thisutils::log_message("[cluster-neuronal-lineage] ", "Loading the complete lineage RPCA representation")
embedding <- readRDS(embedding_file)
if (!is.matrix(embedding) || nrow(embedding) != nrow(cells) ||
  ncol(embedding) < max(dims) || anyDuplicated(rownames(embedding)) ||
  !identical(rownames(embedding), cells$Cells) || any(!is.finite(embedding))) {
  stop("Lineage RPCA representation failed dimension, order or value checks")
}
embedding <- embedding[, dims, drop = FALSE]

save_checkpoint <- function(object, path) {
  temporary <- paste0(path, ".tmp.", Sys.getpid())
  saveRDS(object, temporary, compress = FALSE)
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    stop("Could not publish lineage clustering checkpoint: ", path)
  }
}

validate_neighbor <- function(neighbor) {
  is(neighbor, "Neighbor") &&
    identical(dim(neighbor@nn.idx), c(nrow(cells), k_neighbors)) &&
    identical(dim(neighbor@nn.dist), c(nrow(cells), k_neighbors)) &&
    identical(as.character(neighbor@cell.names), cells$Cells) &&
    identical(attr(neighbor, "k_neighbors", exact = TRUE), k_neighbors) &&
    all(neighbor@nn.idx >= 1L) && all(neighbor@nn.idx <= nrow(cells)) &&
    all(is.finite(neighbor@nn.dist))
}

if (file.exists(neighbor_file)) {
  thisutils::log_message("[cluster-neuronal-lineage] ", "Loading the lineage kNN checkpoint")
  neighbor <- readRDS(neighbor_file)
  if (!validate_neighbor(neighbor)) {
    stop("Existing lineage kNN checkpoint failed validation")
  }
} else {
  thisutils::log_message("[cluster-neuronal-lineage] ", "Computing the complete lineage kNN checkpoint")
  set.seed(seed)
  neighbor <- FindNeighbors(
    object = embedding, k.param = k_neighbors, return.neighbor = TRUE,
    compute.SNN = FALSE, nn.method = "annoy", annoy.metric = "euclidean",
    n.trees = 50L, verbose = TRUE
  )
  attr(neighbor, "k_neighbors") <- k_neighbors
  if (!validate_neighbor(neighbor)) {
    stop("New lineage kNN checkpoint failed validation")
  }
  save_checkpoint(neighbor, neighbor_file)
}
rm(embedding)
gc()

validate_graph <- function(graph) {
  is(graph, "Graph") && identical(dim(graph), c(nrow(cells), nrow(cells))) &&
    identical(rownames(graph), cells$Cells) &&
    identical(colnames(graph), cells$Cells) &&
    identical(attr(graph, "k_neighbors", exact = TRUE), k_neighbors) &&
    all(is.finite(graph@x)) && all(graph@x > 0) && isTRUE(isSymmetric(graph))
}

if (file.exists(graph_file)) {
  thisutils::log_message("[cluster-neuronal-lineage] ", "Loading the lineage graph checkpoint")
  clustering_graph <- readRDS(graph_file)
  if (!validate_graph(clustering_graph)) {
    stop("Existing lineage graph checkpoint failed validation")
  }
} else {
  thisutils::log_message("[cluster-neuronal-lineage] ", "Constructing the complete symmetric lineage graph")
  neighbor_index <- neighbor@nn.idx
  neighbor_distance <- neighbor@nn.dist
  row_scale <- pmax(neighbor_distance[, ncol(neighbor_distance)], .Machine$double.eps)
  source_index <- rep(seq_len(nrow(neighbor_index)), each = k_neighbors)
  target_index <- as.integer(t(neighbor_index))
  edge_distance <- as.numeric(t(neighbor_distance))
  edge_scale <- rep(row_scale, each = k_neighbors)
  keep <- source_index != target_index
  directed_graph <- sparseMatrix(
    i = source_index[keep], j = target_index[keep],
    x = exp(-((edge_distance[keep] / edge_scale[keep])^2)),
    dims = c(nrow(cells), nrow(cells)),
    dimnames = list(cells$Cells, cells$Cells), giveCsparse = TRUE
  )
  clustering_graph <- drop0(directed_graph + t(directed_graph))
  clustering_graph <- as(clustering_graph, "Graph")
  attr(clustering_graph, "k_neighbors") <- k_neighbors
  if (!validate_graph(clustering_graph)) {
    stop("New lineage graph checkpoint failed validation")
  }
  save_checkpoint(clustering_graph, graph_file)
}
rm(neighbor)
gc()

resolution_key <- function(value) {
  gsub("[.]", "_", format(value, nsmall = 1L, trim = TRUE))
}
cluster_column <- function(run_seed, resolution) {
  paste0("Seed", run_seed, "_Resolution", resolution_key(resolution))
}
run_clustering <- function(run_seed, resolution) {
  thisutils::log_message("[cluster-neuronal-lineage] ", paste("Clustering seed", run_seed, "resolution", resolution))
  set.seed(run_seed)
  result <- FindClusters(
    clustering_graph,
    algorithm = 1L, resolution = resolution,
    n.start = 10L, n.iter = 10L, random.seed = run_seed,
    group.singletons = TRUE, verbose = TRUE
  )
  as.character(result[[1L]])
}

assignments <- copy(cells)
for (resolution in resolutions) {
  column <- cluster_column(seed, resolution)
  assignments[, (column) := run_clustering(seed, resolution)]
}
for (run_seed in setdiff(stability_seeds, seed)) {
  column <- cluster_column(run_seed, primary_resolution)
  assignments[, (column) := run_clustering(run_seed, primary_resolution)]
}
primary_column <- cluster_column(seed, primary_resolution)
primary_numeric <- suppressWarnings(as.integer(assignments[[primary_column]]))
if (anyNA(primary_numeric)) {
  stop("Primary lineage cluster labels are not integer-valued")
}
assignments[, Lineage_Cluster := paste0(
  lineage_id, "_", sprintf("%03d", primary_numeric)
)]

choose2 <- function(value) value * (value - 1) / 2
adjusted_rand_index <- function(left, right) {
  contingency <- table(left, right)
  total <- sum(contingency)
  cell_sum <- sum(choose2(as.numeric(contingency)))
  left_sum <- sum(choose2(rowSums(contingency)))
  right_sum <- sum(choose2(colSums(contingency)))
  expected <- left_sum * right_sum / choose2(total)
  maximum <- (left_sum + right_sum) / 2
  if (identical(maximum, expected)) 1 else (cell_sum - expected) / (maximum - expected)
}

assignment_columns <- grep("^Seed", names(assignments), value = TRUE)
summary <- rbindlist(lapply(assignment_columns, function(column) {
  sizes <- table(assignments[[column]])
  data.table(
    Assignment = column,
    Seed = as.integer(sub("^Seed([0-9]+)_.*$", "\\1", column)),
    Resolution = as.numeric(gsub("_", ".", sub("^.*_Resolution", "", column))),
    Clusters = length(sizes),
    Smallest_Cluster_Cells = min(sizes),
    Median_Cluster_Cells = stats::median(as.numeric(sizes)),
    Largest_Cluster_Cells = max(sizes)
  )
}))
seed_columns <- vapply(
  stability_seeds, cluster_column, character(1L),
  resolution = primary_resolution
)
pairs <- utils::combn(seed_columns, 2L, simplify = FALSE)
seed_stability <- rbindlist(lapply(pairs, function(columns) {
  data.table(
    Comparison_Type = "seed_at_primary_resolution",
    Assignment_A = columns[[1L]], Assignment_B = columns[[2L]],
    Adjusted_Rand_Index = adjusted_rand_index(
      assignments[[columns[[1L]]]], assignments[[columns[[2L]]]]
    ),
    Cells = nrow(assignments)
  )
}))
resolution_columns <- vapply(
  resolutions, cluster_column, character(1L),
  run_seed = seed
)
resolution_pairs <- utils::combn(
  resolution_columns, 2L,
  simplify = FALSE
)
resolution_stability <- rbindlist(lapply(resolution_pairs, function(columns) {
  data.table(
    Comparison_Type = "resolution_at_primary_seed",
    Assignment_A = columns[[1L]], Assignment_B = columns[[2L]],
    Adjusted_Rand_Index = adjusted_rand_index(
      assignments[[columns[[1L]]]], assignments[[columns[[2L]]]]
    ),
    Cells = nrow(assignments)
  )
}))
stability <- rbindlist(list(seed_stability, resolution_stability))

fwrite(assignments, assignment_file, sep = "\t", quote = FALSE)
fwrite(summary, summary_file, sep = "\t", quote = FALSE)
fwrite(stability, stability_file, sep = "\t", quote = FALSE)
parameters <- data.table(
  Parameter = c(
    "lineage_id", "cells", "dimensions", "neighbors", "random_seed",
    "resolutions", "primary_resolution", "stability_seeds",
    "clustering_algorithm", "global_cluster_preserved"
  ),
  Value = c(
    lineage_id, nrow(cells), paste(range(dims), collapse = "-"),
    k_neighbors, seed, paste(resolutions, collapse = ","),
    primary_resolution, paste(stability_seeds, collapse = ","),
    "Seurat FindClusters algorithm 1", "true"
  )
)
fwrite(parameters, parameter_file, sep = "\t", quote = FALSE)
input_manifest <- data.table(
  File = required_files, Size_Bytes = file.info(required_files)$size,
)
fwrite(input_manifest, input_manifest_file, sep = "\t", quote = FALSE)
outputs_before_manifest <- setdiff(material_outputs, output_manifest_file)
output_manifest <- data.table(
  File = outputs_before_manifest,
  Size_Bytes = file.info(outputs_before_manifest)$size,,
  Lineage_ID = lineage_id, Cells = nrow(cells)
)
fwrite(output_manifest, output_manifest_file, sep = "\t", quote = FALSE)
primary_clusters <- summary[Assignment == primary_column, Clusters]
if (length(primary_clusters) != 1L) {
  stop("Primary clustering summary must contain exactly one assignment")
}
thisutils::log_message("[cluster-neuronal-lineage] ", paste(
  "Completed", lineage_id, "lineage clustering with",
  primary_clusters,
  "primary clusters"
))
