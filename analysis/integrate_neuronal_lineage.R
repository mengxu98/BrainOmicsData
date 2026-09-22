#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(future)
  library(Seurat)
  library(SeuratObject)
})

source("functions/data_paths.R")
source("functions/processed_object.R")
source("functions/integration.R")


integration_dir <- brainomics_data_path("integration_25")
main_celltype_dir <- file.path(
  integration_dir, "annotation"
)
assignment_file <- file.path(
  main_celltype_dir, "celltype_assignments.rds"
)
audit_file <- file.path(integration_dir, "lineage_analysis_20260825_v2", "lineage_selection.tsv")
lineage_id <- Sys.getenv("BRAINOMICS_LINEAGE_ID", unset = "")
lineage_celltype <- Sys.getenv("BRAINOMICS_LINEAGE_CELLTYPE", unset = "")
analysis_role <- Sys.getenv("BRAINOMICS_LINEAGE_ANALYSIS_ROLE", unset = "")
if (any(!nzchar(c(
  lineage_id, lineage_celltype, analysis_role
)))) {
  stop(
    "Lineage identity, cell type, and BRAINOMICS_LINEAGE_ANALYSIS_ROLE ",
    "are required"
  )
}
if (!analysis_role %chin% c(
  "primary", "complete_sensitivity", "excluded_sensitivity"
)) {
  stop(
    "BRAINOMICS_LINEAGE_ANALYSIS_ROLE must be primary, ",
    "complete_sensitivity or excluded_sensitivity"
  )
}
if (!grepl("^[a-z][a-z0-9_]*$", lineage_id)) {
  stop("BRAINOMICS_LINEAGE_ID must be a lowercase filesystem identifier")
}
output_dir <- brainomics_lineage_path(lineage_id, analysis_role)
excluded_global_clusters <- trimws(strsplit(
  Sys.getenv("BRAINOMICS_LINEAGE_EXCLUDED_GLOBAL_CLUSTERS", unset = ""),
  ",",
  fixed = TRUE
)[[1L]])
excluded_global_clusters <- unique(excluded_global_clusters[
  nzchar(excluded_global_clusters)
])
if (length(excluded_global_clusters) > 0L &&
  any(!grepl("^C[0-9]{2}$", excluded_global_clusters))) {
  stop(
    "BRAINOMICS_LINEAGE_EXCLUDED_GLOBAL_CLUSTERS must contain ",
    "comma-separated resolution-2 cluster IDs"
  )
}

input_candidates <- c(
  file.path(integration_dir, "objects_filtered.rds"),
  file.path(integration_dir, "r_checkpoints", "objects_hvg.rds")
)
existing_input_candidates <- input_candidates[file.exists(input_candidates)]
if (length(existing_input_candidates) == 0L) {
  stop(
    "No global filtered-count checkpoint is available: ",
    paste(input_candidates, collapse = ", ")
  )
}
input_file <- existing_input_candidates[[1L]]
required_files <- c(input_file, assignment_file, audit_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Missing lineage-integration input: ", paste(missing_files, collapse = ", "))
}

seed <- as.integer(Sys.getenv("BRAINOMICS_LINEAGE_SEED", unset = "20260730"))
dims <- seq_len(as.integer(Sys.getenv(
  "BRAINOMICS_LINEAGE_DIMENSIONS",
  unset = "50"
)))
n_features <- as.integer(Sys.getenv(
  "BRAINOMICS_LINEAGE_VARIABLE_FEATURES",
  unset = "3000"
))
rpca_k_weight <- as.integer(Sys.getenv(
  "BRAINOMICS_LINEAGE_RPCA_K_WEIGHT",
  unset = "100"
))
expected_cells <- as.numeric(Sys.getenv(
  "BRAINOMICS_LINEAGE_EXPECTED_CELLS",
  unset = "0"
))
if (!is.finite(seed) || seed < 1L || length(dims) < 10L ||
  max(dims) > 100L || !is.finite(n_features) || n_features < 1000L ||
  !is.finite(rpca_k_weight) || rpca_k_weight < 10L ||
  !is.finite(expected_cells) || expected_cells < 0) {
  stop("Invalid lineage seed, dimensions, variable features or expected cells")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
checkpoint_dir <- file.path(output_dir, "checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
pca_file <- file.path(checkpoint_dir, "objects_lineage_pca.rds")
rpca_file <- file.path(output_dir, "objects_lineage_integrated.rds")
cell_file <- file.path(output_dir, "lineage_cells.tsv.gz")
feature_file <- file.path(output_dir, "lineage_variable_features.tsv")
parameter_file <- file.path(output_dir, "lineage_integration_parameters.tsv")
version_file <- file.path(output_dir, "lineage_integration_versions.tsv")
dataset_file <- file.path(output_dir, "lineage_dataset_counts.tsv")
timing_file <- file.path(output_dir, "lineage_stage_timing.tsv")
input_manifest_file <- file.path(output_dir, "lineage_input_manifest.tsv")
output_manifest_file <- file.path(output_dir, "lineage_output_manifest.tsv")
raw_pca_file <- file.path(output_dir, "lineage_raw_pca_latent.rds")
rpca_latent_file <- file.path(output_dir, "lineage_rpca_latent.rds")
umap_file <- file.path(output_dir, "lineage_rpca_umap.tsv.gz")

future::plan(future::sequential)
if (!identical(as.integer(future::nbrOfWorkers()), 1L)) {
  stop("Neuronal-lineage integration requires one in-process future worker")
}
options(future.globals.maxSize = Inf)

run_id <- format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC")
record_timing <- function(stage, event, started = NULL) {
  row <- data.frame(
    Run_ID = run_id,
    Lineage_ID = lineage_id,
    Stage = stage,
    Event = event,
    Timestamp_UTC = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    Elapsed_Seconds = if (is.null(started)) {
      NA_real_
    } else {
      as.numeric(
        difftime(Sys.time(), started, units = "secs")
      )
    },
    stringsAsFactors = FALSE
  )
  write.table(
    row, timing_file,
    sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = !file.exists(timing_file), append = file.exists(timing_file)
  )
  invisible(row)
}

assignment <- as.data.table(readRDS(assignment_file))
audit <- fread(audit_file, sep = "\t")
required_assignment_columns <- c("Cells", "Cluster", "CellType")
if (!all(required_assignment_columns %in% names(assignment)) ||
  anyNA(assignment[, ..required_assignment_columns]) ||
  anyDuplicated(assignment$Cells)) {
  stop("Main-cell-type assignments failed the complete sidecar contract")
}
required_audit_columns <- c(
  "Cluster", "Cells", "Main_CellType",
  "Exclude_From_Primary_Lineage_Analysis", "Exclude_From_Subtype_Inference",
  "Sensitivity_Set"
)
if (!all(required_audit_columns %in% names(audit)) || nrow(audit) != 75L ||
  anyDuplicated(audit$Cluster) || anyNA(audit[, ..required_audit_columns]) ||
  !is.logical(audit$Exclude_From_Primary_Lineage_Analysis)) {
  stop("Main-cell-type audit failed its 75-cluster schema contract")
}
assignment_audit <- assignment[, .(Cells = .N, Main_CellType = unique(CellType)),
  by = Cluster
]
audit_index <- match(assignment_audit$Cluster, audit$Cluster)
if (anyNA(audit_index) ||
  !identical(as.integer(assignment_audit$Cells), as.integer(audit$Cells[audit_index])) ||
  !identical(
    as.character(assignment_audit$Main_CellType),
    as.character(audit$Main_CellType[audit_index])
  )) {
  stop("Complete-cell assignment differs from the formal cluster audit")
}
required_excluded_clusters <- audit[
  Main_CellType == lineage_celltype &
    Exclude_From_Primary_Lineage_Analysis == TRUE,
  Cluster
]
required_subtype_sensitivity_exclusions <- audit[
  Main_CellType == lineage_celltype &
    Exclude_From_Subtype_Inference == TRUE &
    Sensitivity_Set != "none",
  Cluster
]
if (analysis_role == "primary" &&
  !setequal(excluded_global_clusters, required_excluded_clusters)) {
  stop(
    "Primary lineage exclusions differ from the formal audit: required ",
    paste(required_excluded_clusters, collapse = ","), "; received ",
    paste(excluded_global_clusters, collapse = ",")
  )
}
if (analysis_role == "complete_sensitivity" &&
  length(excluded_global_clusters) > 0L) {
  stop("Complete-lineage sensitivity runs must retain every assigned cell")
}
if (analysis_role == "excluded_sensitivity" &&
  (length(excluded_global_clusters) == 0L || !setequal(
    excluded_global_clusters, required_subtype_sensitivity_exclusions
  ))) {
  stop(
    "Excluded-lineage sensitivity clusters differ from the formal subtype ",
    "uncertainty audit: required ",
    paste(required_subtype_sensitivity_exclusions, collapse = ","),
    "; received ", paste(excluded_global_clusters, collapse = ",")
  )
}
lineage_assignment <- assignment[CellType == lineage_celltype]
unknown_excluded_clusters <- setdiff(
  excluded_global_clusters, unique(lineage_assignment$Cluster)
)
if (length(unknown_excluded_clusters) > 0L) {
  stop(
    "Requested excluded clusters are absent from the lineage: ",
    paste(unknown_excluded_clusters, collapse = ", ")
  )
}
if (length(excluded_global_clusters) > 0L) {
  lineage_assignment <- lineage_assignment[
    !Cluster %chin% excluded_global_clusters
  ]
}
if (nrow(lineage_assignment) == 0L) {
  stop("No cells were assigned to the requested lineage")
}
if (expected_cells > 0 && nrow(lineage_assignment) != expected_cells) {
  stop(
    "Requested lineage contains ", nrow(lineage_assignment),
    " cells, expected ", expected_cells
  )
}

save_lineage_checkpoint <- function(object, stage, path) {
  object@misc$BrainOmics_Lineage_Integration <- list(
    Lineage_ID = lineage_id,
    Lineage_CellType = lineage_celltype,
    Stage = stage,
    Seed = seed,
    Dimensions = dims,
    Variable_Features = n_features,
    RPCA_K_Weight = rpca_k_weight,
    Batch_Model = "Integration_Batch_ID: dataset",
    Analysis_Role = analysis_role,
    Excluded_Global_Clusters = excluded_global_clusters,
    Cells = ncol(object),
    Features = nrow(object)
  )
  save_processed_object(object, path, validate_reload = FALSE)
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Lineage checkpoint was not written: ", path)
  }
  invisible(path)
}

validate_checkpoint <- function(object, stage) {
  contract <- object@misc$BrainOmics_Lineage_Integration
  required_reductions <- switch(stage,
    pca = c("pca.lineage", "umap.lineage.unintegrated"),
    rpca = c(
      "pca.lineage", "umap.lineage.unintegrated",
      "integrated.lineage.rpca", "umap.lineage.rpca"
    ),
    stop("Unknown lineage checkpoint stage: ", stage)
  )
  weight_contract_valid <- if (identical(stage, "pca")) {
    # k.weight is not used until RPCA anchor weighting, so a completed PCA
    # checkpoint remains valid when a smaller lineage layer requires a lower
    # k.weight after the explicit layer-size preflight.
    TRUE
  } else {
    identical(as.integer(contract$RPCA_K_Weight), rpca_k_weight)
  }
  checkpoint_checks <- c(
    contract_is_list = is.list(contract),
    lineage_id = identical(as.character(contract$Lineage_ID), lineage_id),
    lineage_celltype = identical(
      as.character(contract$Lineage_CellType), lineage_celltype
    ),
    stage = identical(as.character(contract$Stage), stage),
    seed = identical(as.integer(contract$Seed), seed),
    dimensions = identical(
      as.integer(contract$Dimensions), as.integer(dims)
    ),
    variable_feature_contract = identical(
      as.integer(contract$Variable_Features), n_features
    ),
    rpca_k_weight = weight_contract_valid,
    batch_model = identical(
      as.character(contract$Batch_Model),
      "Integration_Batch_ID: dataset"
    ),
    analysis_role = identical(
      as.character(contract$Analysis_Role), analysis_role
    ),
    excluded_global_clusters = setequal(
      as.character(contract$Excluded_Global_Clusters),
      as.character(excluded_global_clusters)
    ),
    cells = identical(
      as.numeric(contract$Cells), as.numeric(nrow(lineage_assignment))
    ),
    variable_features = length(VariableFeatures(object)) == n_features,
    reductions = all(required_reductions %in% names(object@reductions))
  )
  if (!all(checkpoint_checks)) {
    stop(
      "Lineage checkpoint differs from the requested contract at stage ",
      stage, "; failed fields: ",
      paste(names(checkpoint_checks)[!checkpoint_checks], collapse = ", ")
    )
  }
  invisible(TRUE)
}

if (file.exists(rpca_file)) {
  thisutils::log_message("[integrate-neuronal-lineage] ", "Resuming from the complete RPCA checkpoint")
  object <- load_processed_object(rpca_file)
  validate_checkpoint(object, "rpca")
  completed_stage <- 2L
} else if (file.exists(pca_file)) {
  thisutils::log_message("[integrate-neuronal-lineage] ", "Resuming from the complete PCA checkpoint")
  object <- load_processed_object(pca_file)
  validate_checkpoint(object, "pca")
  completed_stage <- 1L
} else {
  thisutils::log_message("[integrate-neuronal-lineage] ", "Loading the globally filtered complete count object")
  object <- load_processed_object(input_file)
  if (!setequal(lineage_assignment$Cells, intersect(
    lineage_assignment$Cells, colnames(object)
  ))) {
    stop("One or more lineage cells are absent from the filtered count object")
  }
  lineage_cell_order <- colnames(object)[
    colnames(object) %chin% lineage_assignment$Cells
  ]
  lineage_assignment <- lineage_assignment[match(lineage_cell_order, Cells)]
  object <- subset(object, cells = lineage_assignment$Cells)
  lineage_assignment <- lineage_assignment[match(colnames(object), Cells)]
  if (!identical(colnames(object), lineage_assignment$Cells) ||
    ncol(object) != nrow(lineage_assignment)) {
    stop("Lineage subset and assignment order differ")
  }
  object$Global_Cluster <- lineage_assignment$Cluster
  object$Main_CellType <- lineage_assignment$CellType
  completed_stage <- 0L
}

lineage_assignment <- lineage_assignment[match(colnames(object), Cells)]
if (ncol(object) != nrow(lineage_assignment) ||
  anyNA(lineage_assignment$Cells) ||
  !identical(colnames(object), lineage_assignment$Cells) ||
  !identical(unique(as.character(object$Main_CellType)), lineage_celltype)) {
  stop("Lineage object failed cell, order or identity validation")
}
required_batch_columns <- c(
  "Dataset", "Integration_Batch_ID", "Integration_Batch_Model"
)
missing_batch_columns <- setdiff(required_batch_columns, colnames(object[[]]))
if (length(missing_batch_columns) > 0L ||
  anyNA(object$Integration_Batch_ID) ||
  !identical(unique(as.character(object$Integration_Batch_Model)), "dataset") ||
  !identical(
    as.character(object$Integration_Batch_ID),
    paste0("dataset:", as.character(object$Dataset))
  )) {
  stop("Lineage object differs from the audited dataset-batch contract")
}
count_layers <- Layers(object, assay = "RNA", search = "^counts")
layer_cell_lists <- lapply(count_layers, function(layer) {
  Cells(object[["RNA"]], layer = layer)
})
layer_cells <- unlist(layer_cell_lists, use.names = FALSE)
layer_sizes <- lengths(layer_cell_lists)
if (length(count_layers) < 2L || length(layer_cells) != ncol(object) ||
  anyDuplicated(layer_cells) || !setequal(layer_cells, colnames(object))) {
  stop("Lineage count layers do not cover every cell exactly once")
}
if (rpca_k_weight > min(layer_sizes)) {
  stop(
    "RPCA k.weight (", rpca_k_weight,
    ") exceeds the smallest complete lineage dataset layer (",
    min(layer_sizes), " cells); lower BRAINOMICS_LINEAGE_RPCA_K_WEIGHT"
  )
}

if (completed_stage < 1L) {
  started <- Sys.time()
  record_timing("PCA", "start")
  thisutils::log_message("[integrate-neuronal-lineage] ", "Normalizing the complete lineage and selecting lineage HVGs")
  object <- NormalizeData(object, verbose = FALSE)
  object <- FindVariableFeatures(
    object,
    selection.method = "vst", nfeatures = n_features,
    verbose = FALSE
  )
  selected_features <- VariableFeatures(object)
  if (length(selected_features) != n_features || anyNA(selected_features) ||
    anyDuplicated(selected_features)) {
    stop("Lineage variable-feature selection failed its fixed contract")
  }
  object <- ScaleData(object, features = selected_features, verbose = FALSE)
  set.seed(seed)
  object <- RunPCA(
    object,
    features = selected_features, npcs = max(dims),
    reduction.name = "pca.lineage", reduction.key = "PCLINEAGE_",
    seed.use = seed, verbose = FALSE
  )
  set.seed(seed)
  object <- RunUMAP(
    object,
    reduction = "pca.lineage", dims = dims,
    reduction.name = "umap.lineage.unintegrated",
    reduction.key = "UMAPLINEAGEUNINTEGRATED_", seed.use = seed,
    verbose = FALSE
  )
  record_timing("PCA", "end", started)
  save_lineage_checkpoint(object, "pca", pca_file)
} else {
  selected_features <- VariableFeatures(object)
}

if (completed_stage < 2L) {
  started <- Sys.time()
  record_timing("RPCA", "start")
  thisutils::log_message("[integrate-neuronal-lineage] ", "Re-estimating RPCA anchors within the complete lineage")
  set.seed(seed)
  object <- IntegrateLayers(
    object = object,
    method = RPCAIntegration,
    orig.reduction = "pca.lineage",
    new.reduction = "integrated.lineage.rpca",
    dims = dims,
    k.weight = rpca_k_weight,
    verbose = FALSE
  )
  set.seed(seed)
  object <- RunUMAP(
    object,
    reduction = "integrated.lineage.rpca", dims = dims,
    reduction.name = "umap.lineage.rpca",
    reduction.key = "UMAPLINEAGERPCA_", seed.use = seed,
    verbose = FALSE
  )
  record_timing("RPCA", "end", started)
  save_lineage_checkpoint(object, "rpca", rpca_file)
}

fwrite(
  lineage_assignment[, .(Cells, Global_Cluster = Cluster, Main_CellType = CellType)],
  cell_file,
  sep = "\t", quote = FALSE
)
fwrite(
  data.table(Rank = seq_along(selected_features), Gene = selected_features),
  feature_file,
  sep = "\t", quote = FALSE
)
dataset_counts <- as.data.table(table(Dataset = as.character(object$Dataset)))
setnames(dataset_counts, "N", "Cells")
setorder(dataset_counts, -Cells, Dataset)
fwrite(dataset_counts, dataset_file, sep = "\t", quote = FALSE)

raw_pca <- Embeddings(object, reduction = "pca.lineage")[, dims, drop = FALSE]
rpca_latent <- Embeddings(
  object,
  reduction = "integrated.lineage.rpca"
)[, dims, drop = FALSE]
rpca_umap <- Embeddings(object, reduction = "umap.lineage.rpca")
if (!identical(rownames(raw_pca), lineage_assignment$Cells) ||
  !identical(rownames(rpca_latent), lineage_assignment$Cells) ||
  !identical(rownames(rpca_umap), lineage_assignment$Cells) ||
  any(!is.finite(raw_pca)) || any(!is.finite(rpca_latent)) ||
  any(!is.finite(rpca_umap))) {
  stop("Lineage embedding export failed its order or finite-value contract")
}
saveRDS(raw_pca, raw_pca_file, compress = FALSE)
saveRDS(rpca_latent, rpca_latent_file, compress = FALSE)
fwrite(
  data.table(
    Cells = rownames(rpca_umap),
    UMAP_1 = rpca_umap[, 1L],
    UMAP_2 = rpca_umap[, 2L]
  ),
  umap_file,
  sep = "\t", quote = FALSE
)

parameters <- data.table(
  Parameter = c(
    "lineage_id", "lineage_celltype", "cells", "datasets",
    "count_layers", "variable_features", "dimensions", "random_seed",
    "rpca_dimensions", "rpca_k_weight", "smallest_dataset_layer_cells",
    "batch_model", "analysis_role", "excluded_global_clusters", "input_feature_filter",
    "global_cluster_preserved"
  ),
  Value = c(
    lineage_id, lineage_celltype, ncol(object),
    length(unique(as.character(object$Dataset))), length(count_layers),
    n_features, paste(range(dims), collapse = "-"), seed,
    paste(range(dims), collapse = "-"), rpca_k_weight, min(layer_sizes),
    "Integration_Batch_ID: dataset",
    analysis_role,
    if (length(excluded_global_clusters) == 0L) {
      "none"
    } else {
      paste(excluded_global_clusters, collapse = ",")
    },
    paste0(
      "reuse global filtered-count checkpoint ", basename(input_file),
      "; no new source-level QC or feature filtering"
    ),
    "true; stored as Global_Cluster"
  )
)
fwrite(parameters, parameter_file, sep = "\t", quote = FALSE)
versions <- data.table(
  Package = c("R", "Seurat", "SeuratObject"),
  Version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    as.character(utils::packageVersion("Seurat")),
    as.character(utils::packageVersion("SeuratObject"))
  )
)
fwrite(versions, version_file, sep = "\t", quote = FALSE)

input_manifest <- data.table(
  File = required_files,
  Size_Bytes = file.info(required_files)$size,
)
fwrite(input_manifest, input_manifest_file, sep = "\t", quote = FALSE)
material_outputs <- c(
  pca_file, rpca_file, cell_file, feature_file,
  parameter_file, version_file, dataset_file, timing_file, input_manifest_file,
  raw_pca_file, rpca_latent_file, umap_file
)
if (any(!file.exists(material_outputs)) || any(file.info(material_outputs)$size <= 0)) {
  stop("One or more lineage integration outputs are absent or empty")
}
output_manifest <- data.table(
  File = material_outputs,
  Size_Bytes = file.info(material_outputs)$size,,
  Lineage_ID = lineage_id,
  Cells = ncol(object),
  Features = nrow(object)
)
fwrite(output_manifest, output_manifest_file, sep = "\t", quote = FALSE)
thisutils::log_message("[integrate-neuronal-lineage] ", paste(
  "Completed", lineage_id, "lineage integration:",
  format(ncol(object), big.mark = ","), "cells across",
  length(unique(as.character(object$Dataset))), "datasets"
))
