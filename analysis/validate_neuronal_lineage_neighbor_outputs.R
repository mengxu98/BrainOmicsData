#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

source("functions/data_paths.R")
source("functions/processed_object.R")

lineage_id <- Sys.getenv("BRAINOMICS_LINEAGE_ID", unset = "")
lineage_dir <- brainomics_lineage_path(lineage_id)
if (!dir.exists(lineage_dir) ||
  !lineage_id %chin% c("excitatory", "inhibitory")) {
  stop("A complete excitatory or inhibitory lineage directory is required")
}

output_dir <- file.path(lineage_dir, "neighborhood_evidence")
paths <- list(
  assignments = file.path(
    lineage_dir, "clustering", "lineage_cluster_assignments.tsv.gz"
  ),
  graph = file.path(
    lineage_dir, "clustering", "checkpoints", "lineage_rpca_graph.rds"
  ),
  rpca = file.path(lineage_dir, "lineage_rpca_latent.rds"),
  umap = file.path(lineage_dir, "lineage_rpca_umap.tsv.gz"),
  neighborhoods = file.path(output_dir, "lineage_cluster_neighborhoods.tsv"),
  centroids = file.path(output_dir, "lineage_cluster_centroids.tsv"),
  contract = file.path(output_dir, "lineage_neighborhood_audit_contract.tsv"),
  manifest = file.path(output_dir, "lineage_neighborhood_audit_manifest.tsv")
)
missing <- unlist(paths)[!file.exists(unlist(paths))]
if (length(missing) > 0L ||
  !file.exists(file.path(lineage_dir, "_VALIDATED"))) {
  stop(
    "Lineage-neighborhood validation inputs are missing: ",
    paste(missing, collapse = ", ")
  )
}

assignments <- fread(paths$assignments, sep = "\t")
neighborhoods <- fread(paths$neighborhoods, sep = "\t")
centroids <- fread(paths$centroids, sep = "\t")
contract <- fread(paths$contract, sep = "\t")
manifest <- fread(paths$manifest, sep = "\t")
if (!identical(names(contract), c("Field", "Value")) ||
  anyDuplicated(contract$Field) ||
  !identical(names(manifest), c("File", "Size_Bytes")) ||
  anyDuplicated(manifest$File)) {
  stop("Lineage-neighborhood contract or manifest schema changed")
}
contract_values <- setNames(contract$Value, contract$Field)
expected_manifest_files <- c(
  "lineage_cluster_neighborhoods.tsv",
  "lineage_cluster_centroids.tsv",
  "lineage_neighborhood_audit_contract.tsv"
)
if (!identical(sort(manifest$File), sort(expected_manifest_files))) {
  stop("Lineage-neighborhood manifest membership changed")
}
manifest_paths <- file.path(output_dir, manifest$File)
if (any(!file.exists(manifest_paths)) ||
  any(as.numeric(manifest$Size_Bytes) != file.info(manifest_paths)$size)) {
  stop("Lineage-neighborhood manifest size or SHA-256 validation failed")
}

required_assignment <- c("Cells", "Lineage_Cluster")
required_neighbor <- c(
  "Lineage_Cluster", "Rank", "Graph_Neighbor_Cluster",
  "Graph_Weight_Fraction", "RPCA_Neighbor_Cluster",
  "RPCA_Centroid_Distance", "UMAP_Neighbor_Cluster",
  "UMAP_Centroid_Distance"
)
if (any(!required_assignment %in% names(assignments)) ||
  any(!required_neighbor %in% names(neighborhoods)) ||
  anyNA(assignments[, ..required_assignment]) ||
  anyNA(neighborhoods[, ..required_neighbor]) ||
  anyDuplicated(assignments$Cells) ||
  anyDuplicated(neighborhoods[, .(Lineage_Cluster, Rank)])) {
  stop("Lineage-neighborhood cell or neighbor schema changed")
}
cluster_ids <- sort(unique(assignments$Lineage_Cluster))
expected_ranks <- seq_len(min(5L, length(cluster_ids) - 1L))
rank_contract <- neighborhoods[
  , .(Ranks = paste(sort(unique(Rank)), collapse = ",")),
  by = Lineage_Cluster
]
if (!identical(sort(unique(neighborhoods$Lineage_Cluster)), cluster_ids) ||
  nrow(neighborhoods) != length(cluster_ids) * length(expected_ranks) ||
  any(rank_contract$Ranks != paste(expected_ranks, collapse = ",")) ||
  any(!neighborhoods$Graph_Neighbor_Cluster %chin% cluster_ids) ||
  any(!neighborhoods$RPCA_Neighbor_Cluster %chin% cluster_ids) ||
  any(!neighborhoods$UMAP_Neighbor_Cluster %chin% cluster_ids) ||
  any(neighborhoods$Lineage_Cluster == neighborhoods$Graph_Neighbor_Cluster) ||
  any(neighborhoods$Lineage_Cluster == neighborhoods$RPCA_Neighbor_Cluster) ||
  any(neighborhoods$Lineage_Cluster == neighborhoods$UMAP_Neighbor_Cluster) ||
  any(!is.finite(as.matrix(neighborhoods[, -c(1L, 3L, 5L, 7L)])))) {
  stop("Lineage-neighborhood rank, identity or numeric contract failed")
}
if (!all(c("Lineage_Cluster", "Representation", "Coordinate", "Value") %in%
  names(centroids)) || anyNA(centroids) ||
  nrow(centroids) != length(cluster_ids) * 52L ||
  !identical(sort(unique(centroids$Lineage_Cluster)), cluster_ids) ||
  !identical(sort(unique(centroids$Representation)), c("RPCA", "UMAP")) ||
  any(!is.finite(centroids$Value))) {
  stop("Lineage-neighborhood centroid contract failed")
}

expected_contract <- c(
  Lineage_ID = lineage_id,
  Cells = as.character(nrow(assignments)),
  Clusters = as.character(length(cluster_ids)),
  RPCA_Dimensions = "50",
  UMAP_Dimensions = "2",
)
if (any(!names(expected_contract) %chin% names(contract_values)) ||
  any(contract_values[names(expected_contract)] != expected_contract)) {
  stop("Lineage-neighborhood source identity contract failed")
}

message(
  "[lineage-neighborhood-validation] passed for ", nrow(assignments), " ",
  lineage_id, " cells and ", length(cluster_ids), " clusters"
)
