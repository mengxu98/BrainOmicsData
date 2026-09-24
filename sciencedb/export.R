#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) {
    return(default)
  }
  args[[idx + 1]]
}

repo_dir <- normalizePath(value_after("--repo-dir", "."), mustWork = TRUE)
source(file.path(repo_dir, "functions", "metadata_schema.R"))
source(file.path(repo_dir, "functions", "sample_schema.R"))
source(file.path(repo_dir, "functions", "utils.R"))
source(file.path(repo_dir, "sciencedb", "manifest.R"))
source(file.path(repo_dir, "sciencedb", "package_metadata.R"))
source(file.path(repo_dir, "sciencedb", "readers.R"))

object_file <- normalizePath(
  value_after("--object-file", "../../data/BrainOmicsData/integration_25/objects_merged.rds"),
  mustWork = TRUE
)
metadata_object_file <- normalizePath(
  value_after("--metadata-object-file", object_file),
  mustWork = TRUE
)
celltype_assignment_file <- value_after(
  "--celltype-assignment-file",
  "results/annotation/celltype_assignments.rds"
)
if (nzchar(celltype_assignment_file)) {
  celltype_assignment_file <- normalizePath(
    celltype_assignment_file,
    mustWork = TRUE
  )
}
out_dir <- value_after("--out-dir", "../../data/BrainOmicsData/ScienceDB")
assay <- value_after("--assay", "RNA")
reduction_pca <- value_after("--pca-reduction", "integrated.rpca")
reduction_umap <- value_after("--umap-reduction", "umap.rpca")
reduction_unintegrated_umap <- value_after("--unintegrated-umap-reduction", "umap.unintegrated")
lisi_file <- value_after(
  "--lisi-file",
  file.path(dirname(object_file), "lisi_results.rds")
)
rpca_latent_file <- value_after(
  "--rpca-latent-file",
  file.path(dirname(object_file), "annotation", "reductions", "rpca_latent.rds")
)
umap_plot_file <- value_after(
  "--umap-plot-file",
  file.path(dirname(object_file), "evaluation", "umap_plot_data.rds")
)
overwrite <- "--overwrite" %in% args
reuse_expression <- "--reuse-expression" %in% args

dirs <- list(
  expression = file.path(out_dir, "expression"),
  metadata = file.path(out_dir, "metadata"),
  embeddings = file.path(out_dir, "embeddings"),
  validation = file.path(out_dir, "validation"),
  objects = file.path(out_dir, "objects"),
  scripts = file.path(out_dir, "scripts"),
  provenance = file.path(out_dir, "provenance")
)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

write_tsv <- function(x, file, gzip = grepl("\\.gz$", file)) {
  con <- if (gzip) gzfile(file, "wt") else file(file, "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "NA")
}

write_tsv_no_header <- function(x, file, gzip = grepl("\\.gz$", file)) {
  con <- if (gzip) gzfile(file, "wt") else file(file, "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, na = "NA")
}

safe_remove <- function(files) {
  if (overwrite) unlink(files[file.exists(files)])
}

message("Loading expression Seurat object: ", object_file)
object <- readRDS(object_file)
if (!inherits(object, "Seurat")) stop("object is not a Seurat object")
if (!requireNamespace("SeuratObject", quietly = TRUE)) stop("SeuratObject is required")
if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix is required")
if (!assay %in% names(object@assays)) stop("missing assay: ", assay)

expression_meta <- object@meta.data
if (!"Dataset" %in% colnames(expression_meta)) {
  stop("expression object metadata is missing Dataset")
}

counts_layers <- SeuratObject::Layers(object[[assay]], search = "^counts")
if (length(counts_layers) == 0) stop("no counts layers found in assay ", assay)
message("Counts layers: ", paste(counts_layers, collapse = ", "))

layer_info <- lapply(counts_layers, function(layer) {
  mat <- SeuratObject::LayerData(object[[assay]], layer = layer, fast = FALSE)
  mat <- methods::as(mat, "dgCMatrix")
  layer_cells <- colnames(mat)
  datasets <- unique(as.character(expression_meta[layer_cells, "Dataset"]))
  datasets <- datasets[!is.na(datasets) & nzchar(datasets)]
  if (length(datasets) != 1L) {
    stop("Counts layer does not map to exactly one source dataset: ", layer)
  }
  list(
    layer = layer,
    dataset = datasets[[1L]],
    genes = rownames(mat),
    cells = layer_cells,
    nrow = nrow(mat),
    ncol = ncol(mat),
    nnz = length(mat@x)
  )
})

features <- unique(unlist(lapply(layer_info, `[[`, "genes"), use.names = FALSE))
cells <- unlist(lapply(layer_info, `[[`, "cells"), use.names = FALSE)
if (anyDuplicated(cells) || anyDuplicated(vapply(
  layer_info, `[[`, character(1), "dataset"
))) {
  stop("Expression layers must contain unique cells and one unique dataset each")
}
frozen_assignments <- NULL
if (nzchar(celltype_assignment_file)) {
  message(
    "Preflighting frozen main-cell-type assignments: ",
    celltype_assignment_file
  )
  frozen_assignments <- readRDS(celltype_assignment_file)
  required_assignment_columns <- c("Cells", "Cluster", "CellType")
  if (!is.data.frame(frozen_assignments) ||
    any(!required_assignment_columns %in% colnames(frozen_assignments)) ||
    nrow(frozen_assignments) != length(cells) ||
    anyDuplicated(frozen_assignments$Cells) ||
    !setequal(as.character(frozen_assignments$Cells), cells)) {
    stop("Frozen main-cell-type assignments fail the complete-cell contract")
  }
  assignment_index <- match(cells, as.character(frozen_assignments$Cells))
  frozen_assignments <- frozen_assignments[
    assignment_index, required_assignment_columns,
    drop = FALSE
  ]
  if (!identical(as.character(frozen_assignments$Cells), cells) ||
    anyNA(frozen_assignments$Cluster) || anyNA(frozen_assignments$CellType) ||
    any(!nzchar(as.character(frozen_assignments$Cluster))) ||
    any(!nzchar(as.character(frozen_assignments$CellType)))) {
    stop("Frozen main-cell-type assignments could not be aligned to expression cells")
  }
  rm(assignment_index)
}
layer_nnz <- vapply(layer_info, `[[`, numeric(1), "nnz")
total_nnz <- sum(as.double(layer_nnz))
if (any(layer_nnz > .Machine$integer.max)) {
  stop("At least one expression shard exceeds the dgCMatrix index limit")
}
message("Export dimensions: ", length(features), " genes x ", length(cells), " cells; nnz=", total_nnz)

existing_expression <- list.files(
  dirs$expression,
  recursive = TRUE, full.names = TRUE, all.files = TRUE,
  no.. = TRUE
)
if (!reuse_expression) {
  if (length(existing_expression) && !overwrite) {
    stop("expression files already exist; use --overwrite to replace them")
  }
  if (length(existing_expression)) {
    unlink(existing_expression, recursive = TRUE, force = TRUE)
  }
  dir.create(
    file.path(dirs$expression, "shards"),
    recursive = TRUE,
    showWarnings = FALSE
  )
} else if (!length(existing_expression)) {
  stop("--reuse-expression requires an existing expression package")
}

feature_table <- data.frame(
  gene_id = features,
  gene_name = features,
  feature_type = "Gene Expression",
  stringsAsFactors = FALSE
)

write_matrix_market_gz <- function(mat, file, block_columns = 500L) {
  if (!inherits(mat, "dgCMatrix")) stop("Matrix Market writer requires dgCMatrix")
  con <- gzfile(file, "wt", compression = 6L)
  on.exit(close(con), add = TRUE)
  writeLines("%%MatrixMarket matrix coordinate real general", con)
  writeLines("% BrainOmicsData complete dataset-resolved raw-count shard", con)
  writeLines(
    paste(nrow(mat), ncol(mat), format(length(mat@x), scientific = FALSE)),
    con
  )
  for (start in seq.int(1L, ncol(mat), by = block_columns)) {
    end <- min(ncol(mat), start + block_columns - 1L)
    first <- as.double(mat@p[[start]]) + 1
    last <- as.double(mat@p[[end + 1L]])
    if (last < first) next
    positions <- seq.int(first, last)
    columns <- rep.int(
      seq.int(start, end), diff(mat@p[start:(end + 1L)])
    )
    block <- data.frame(
      i = mat@i[positions] + 1L,
      j = columns,
      x = mat@x[positions]
    )
    utils::write.table(
      block, con,
      sep = " ", quote = FALSE, row.names = FALSE,
      col.names = FALSE
    )
    rm(block, positions, columns)
  }
  close(con)
  on.exit(NULL)
}

cell_offset <- 0
shard_rows <- vector("list", length(layer_info))
for (i in seq_along(layer_info)) {
  info <- layer_info[[i]]
  safe_dataset <- gsub("[^A-Za-z0-9_.-]+", "_", info$dataset)
  relative_dir <- file.path(
    "shards", sprintf("%02d_%s", i, safe_dataset)
  )
  shard_rows[[i]] <- data.frame(
    Order = i,
    Dataset = info$dataset,
    Source_Layer = info$layer,
    Relative_Directory = relative_dir,
    Features = info$nrow,
    Cells = info$ncol,
    Nonzero_Values = as.double(info$nnz),
    Cell_Start = cell_offset + 1,
    Cell_End = cell_offset + info$ncol,
    Expression_Layout_Version = sciencedb_expression_layout_version(),
    stringsAsFactors = FALSE
  )
  cell_offset <- cell_offset + info$ncol
}
expected_shard_manifest <- do.call(rbind, shard_rows)
if (cell_offset != length(cells) ||
  sum(expected_shard_manifest$Nonzero_Values) != total_nnz ||
  !identical(expected_shard_manifest$Dataset, vapply(
    layer_info, `[[`, character(1), "dataset"
  ))) {
  stop("Expected expression shard manifest failed the complete-cell contract")
}

if (reuse_expression) {
  manifest_file <- file.path(dirs$expression, "shard_manifest.tsv")
  if (!file.exists(manifest_file)) {
    stop("--reuse-expression requires expression/shard_manifest.tsv")
  }
  shard_manifest <- utils::read.delim(
    manifest_file,
    check.names = FALSE, stringsAsFactors = FALSE
  )
  if (!identical(names(shard_manifest), names(expected_shard_manifest)) ||
    nrow(shard_manifest) != nrow(expected_shard_manifest)) {
    stop("Existing expression shard manifest has an incompatible schema")
  }
  numeric_columns <- c(
    "Order", "Features", "Cells", "Nonzero_Values", "Cell_Start", "Cell_End"
  )
  for (column in numeric_columns) {
    shard_manifest[[column]] <- as.numeric(shard_manifest[[column]])
    expected_shard_manifest[[column]] <- as.numeric(
      expected_shard_manifest[[column]]
    )
  }
  if (!identical(shard_manifest, expected_shard_manifest)) {
    stop("Existing expression shard manifest differs from the source object")
  }
  for (i in seq_along(layer_info)) {
    info <- layer_info[[i]]
    shard_dir <- file.path(
      dirs$expression, shard_manifest$Relative_Directory[[i]]
    )
    shard_files <- file.path(
      shard_dir, c("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz")
    )
    if (any(!file.exists(shard_files)) || any(file.info(shard_files)$size <= 0)) {
      stop("Existing expression shard is incomplete: ", info$dataset)
    }
    gzip_status <- system2("gzip", c("-t", shard_files))
    if (!identical(gzip_status, 0L)) {
      stop("Existing expression shard failed gzip integrity: ", info$dataset)
    }
    matrix_con <- gzfile(shard_files[[1L]], "rt")
    matrix_header <- readLines(matrix_con, n = 3L, warn = FALSE)
    close(matrix_con)
    expected_header <- paste(info$nrow, info$ncol, info$nnz)
    if (length(matrix_header) != 3L ||
      !identical(matrix_header[[3L]], expected_header)) {
      stop("Existing Matrix Market header differs: ", info$dataset)
    }
    feature_con <- gzfile(shard_files[[2L]], "rt")
    feature_lines <- readLines(feature_con, warn = FALSE)
    close(feature_con)
    barcode_con <- gzfile(shard_files[[3L]], "rt")
    barcode_lines <- readLines(barcode_con, warn = FALSE)
    close(barcode_con)
    feature_ids <- sub("\\t.*$", "", feature_lines)
    if (!identical(feature_ids, info$genes) ||
      !identical(barcode_lines, info$cells)) {
      stop("Existing expression feature or cell order differs: ", info$dataset)
    }
    message(
      "Reusing validated complete expression shard ", i, "/",
      length(layer_info), ": ", info$dataset
    )
  }
} else {
  shard_manifest <- expected_shard_manifest
  for (i in seq_along(layer_info)) {
    info <- layer_info[[i]]
    shard_dir <- file.path(
      dirs$expression, shard_manifest$Relative_Directory[[i]]
    )
    dir.create(shard_dir, recursive = TRUE, showWarnings = FALSE)
    message(
      "Writing complete expression shard ", i, "/", length(layer_info),
      ": ", info$dataset
    )
    mat <- SeuratObject::LayerData(
      object[[assay]],
      layer = info$layer, fast = FALSE
    )
    mat <- methods::as(mat, "dgCMatrix")
    if (!identical(rownames(mat), info$genes) ||
      !identical(colnames(mat), info$cells) || length(mat@x) != info$nnz) {
      stop("Expression shard changed after the pre-export audit: ", info$dataset)
    }
    write_matrix_market_gz(mat, file.path(shard_dir, "matrix.mtx.gz"))
    shard_features <- data.frame(
      gene_id = rownames(mat), gene_name = rownames(mat),
      feature_type = "Gene Expression", stringsAsFactors = FALSE
    )
    write_tsv_no_header(
      shard_features, file.path(shard_dir, "features.tsv.gz")
    )
    barcodes_con <- gzfile(file.path(shard_dir, "barcodes.tsv.gz"), "wt")
    writeLines(colnames(mat), barcodes_con)
    close(barcodes_con)
    rm(mat, shard_features)
    gc(FALSE)
  }
  write_tsv(shard_manifest, file.path(dirs$expression, "shard_manifest.tsv"))
}

rm(object, expression_meta, layer_info)
gc(FALSE)

message("Loading metadata Seurat object: ", metadata_object_file)
metadata_object <- readRDS(metadata_object_file)
if (!inherits(metadata_object, "Seurat")) {
  stop("metadata object is not a Seurat object")
}
if (ncol(metadata_object) != length(cells) ||
  !setequal(colnames(metadata_object), cells)) {
  stop("metadata object cells differ from the expression-cell contract")
}
meta <- metadata_object@meta.data
meta[] <- lapply(meta, as.character)
meta <- add_metadata_schema(meta)
meta <- add_sample_schema(meta)

if (!is.null(frozen_assignments)) {
  meta[cells, "Cells"] <- cells
  meta[cells, "Cluster"] <- as.character(frozen_assignments$Cluster)
  meta[cells, "CellType"] <- as.character(frozen_assignments$CellType)
  rm(frozen_assignments)
  gc(FALSE)
}

metadata <- build_sciencedb_cell_metadata(meta, cells)
write_tsv(metadata, file.path(dirs$metadata, "metadata.tsv.gz"))
sample_audit <- sample_schema_audit(meta)
write_tsv(sample_audit, file.path(dirs$metadata, "sample_schema_audit.tsv"))

object_meta <- metadata
rownames(object_meta) <- object_meta$Cells
object_meta <- object_meta[colnames(metadata_object), , drop = FALSE]
if (!identical(rownames(object_meta), colnames(metadata_object))) {
  stop("The lightweight object metadata is not aligned to object cells")
}
metadata_object@meta.data <- object_meta
saveRDS(metadata_object, file.path(dirs$objects, "objects_celltype_plot.rds"), compress = TRUE)

dataset_summary <- build_sciencedb_dataset_summary(meta, metadata, repo_dir)
write_tsv(dataset_summary, file.path(dirs$metadata, "dataset_summary.tsv"))
attrition_summary <- build_sciencedb_attrition_summary(
  repo_dir, dataset_summary
)
write_tsv(
  attrition_summary,
  file.path(dirs$metadata, "dataset_attrition_summary.tsv")
)
verification_dictionary <- brainomics_verification_status_dictionary()
write_tsv(
  verification_dictionary,
  file.path(dirs$metadata, "verification_status_dictionary.tsv")
)
write_tsv(feature_table, file.path(dirs$metadata, "feature_metadata.tsv"))

if (!file.exists(rpca_latent_file) || !file.exists(umap_plot_file)) {
  stop(
    "Validated embedding sidecars are missing: ",
    paste(c(rpca_latent_file, umap_plot_file)[
      !file.exists(c(rpca_latent_file, umap_plot_file))
    ], collapse = ", ")
  )
}
rpca_latent <- readRDS(rpca_latent_file)
if (!is.matrix(rpca_latent) || nrow(rpca_latent) != length(cells) ||
  ncol(rpca_latent) != 50L || !identical(rownames(rpca_latent), cells) ||
  any(!is.finite(rpca_latent))) {
  stop("Validated RPCA latent sidecar differs from the expression-cell contract")
}
rpca_table <- data.frame(
  cell_id = rownames(rpca_latent), rpca_latent,
  check.names = FALSE, stringsAsFactors = FALSE
)
rpca_output_columns <- ncol(rpca_table)
write_tsv(rpca_table, file.path(dirs$embeddings, "integrated_pca.tsv.gz"))
rm(rpca_latent, rpca_table)
gc(FALSE)

umap_plot <- readRDS(umap_plot_file)
required_umap_columns <- c("Cell", "Raw_1", "Raw_2", "RPCA_1", "RPCA_2")
if (!is.data.frame(umap_plot) || nrow(umap_plot) != length(cells) ||
  any(!required_umap_columns %in% colnames(umap_plot)) ||
  !identical(as.character(umap_plot$Cell), cells) ||
  any(!is.finite(as.matrix(umap_plot[, required_umap_columns[-1L]])))) {
  stop("Validated UMAP plotting sidecar differs from the expression-cell contract")
}
write_tsv(
  data.frame(
    cell_id = umap_plot$Cell,
    UMAP_1 = umap_plot$RPCA_1,
    UMAP_2 = umap_plot$RPCA_2,
    stringsAsFactors = FALSE
  ),
  file.path(dirs$embeddings, "integrated_umap.tsv.gz")
)
write_tsv(
  data.frame(
    cell_id = umap_plot$Cell,
    UMAP_1 = umap_plot$Raw_1,
    UMAP_2 = umap_plot$Raw_2,
    stringsAsFactors = FALSE
  ),
  file.path(dirs$embeddings, "unintegrated_umap.tsv.gz")
)
rm(umap_plot)
gc(FALSE)

if (!file.exists(lisi_file)) {
  stop("formal integration-space LISI results are missing: ", lisi_file)
}
lisi <- readRDS(lisi_file)
if (ncol(lisi) != 4L || !setequal(names(lisi), method_levels) ||
  !identical(attr(lisi, "space"), "latent") ||
  any(!cells %in% rownames(lisi))) {
  stop("LISI results are not the validated latent-space four-method output")
}
lisi <- as.data.frame(lisi)[cells, method_levels, drop = FALSE]
lisi <- cbind(cell_id = rownames(lisi), lisi)
write_tsv(lisi, file.path(dirs$validation, "lisi.tsv.gz"))

write_sciencedb_readers(out_dir)

write_tsv(
  build_sciencedb_references(dataset_summary),
  file.path(dirs$provenance, "dataset_manifest.tsv")
)

write_sciencedb_readme(out_dir)

dimension_audit <- data.frame(
  path = c(
    "expression/shard_manifest.tsv", "metadata/metadata.tsv.gz",
    "metadata/dataset_summary.tsv", "metadata/dataset_attrition_summary.tsv",
    "metadata/verification_status_dictionary.tsv",
    "metadata/sample_schema_audit.tsv", "metadata/feature_metadata.tsv",
    "validation/lisi.tsv.gz",
    "objects/objects_celltype_plot.rds"
  ),
  rows = c(
    nrow(shard_manifest), nrow(metadata), nrow(dataset_summary),
    nrow(attrition_summary), nrow(verification_dictionary),
    nrow(sample_audit), nrow(feature_table), nrow(lisi),
    nrow(metadata_object)
  ),
  columns = c(
    ncol(shard_manifest), ncol(metadata), ncol(dataset_summary),
    ncol(attrition_summary), ncol(verification_dictionary),
    ncol(sample_audit), ncol(feature_table), ncol(lisi),
    ncol(metadata_object)
  ),
  nonzero_values = rep(NA_real_, 9L),
  stringsAsFactors = FALSE
)
shard_dimension_rows <- do.call(rbind, lapply(
  seq_len(nrow(shard_manifest)),
  function(i) {
    prefix <- file.path(
      "expression", shard_manifest$Relative_Directory[[i]]
    )
    data.frame(
      path = file.path(
        prefix, c("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz")
      ),
      rows = c(
        shard_manifest$Features[[i]], shard_manifest$Features[[i]],
        shard_manifest$Cells[[i]]
      ),
      columns = c(shard_manifest$Cells[[i]], 3, 1),
      nonzero_values = c(
        shard_manifest$Nonzero_Values[[i]], NA_real_, NA_real_
      ),
      stringsAsFactors = FALSE
    )
  }
))
dimension_audit <- rbind(dimension_audit, shard_dimension_rows)
dimension_audit <- rbind(
  dimension_audit,
  data.frame(
    path = c(
      "embeddings/integrated_pca.tsv.gz",
      "embeddings/integrated_umap.tsv.gz",
      "embeddings/unintegrated_umap.tsv.gz"
    ),
    rows = rep(length(cells), 3L),
    columns = c(rpca_output_columns, 3L, 3L),
    nonzero_values = rep(NA_real_, 3L),
    stringsAsFactors = FALSE
  )
)
write_release_manifest(out_dir, dimension_audit)

message("Clean ScienceDB package written to: ", out_dir)
