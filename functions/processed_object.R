open_h5ad_counts <- function(h5ad_file, matrix_group) {
  if (!requireNamespace("BPCells", quietly = TRUE)) {
    stop(
      "BPCells is required only to stream the complete H5AD input during ",
      "preprocessing"
    )
  }
  if (!file.exists(h5ad_file)) {
    stop("H5AD input does not exist: ", h5ad_file)
  }
  BPCells::open_matrix_anndata_hdf5(
    path = h5ad_file,
    group = matrix_group
  )
}

materialize_sparse_counts <- function(
  matrix,
  expected_features,
  expected_cells,
  expected_feature_names = NULL,
  expected_cell_names = NULL
) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix is required to materialize sparse counts")
  }

  observed <- dim(matrix)
  expected <- c(expected_features, expected_cells)
  if (!identical(as.numeric(observed), as.numeric(expected))) {
    stop(
      "input matrix dimension mismatch: expected ",
      paste(expected, collapse = " x "),
      ", observed ",
      paste(observed, collapse = " x ")
    )
  }
  if (!is.null(expected_feature_names) &&
    !identical(rownames(matrix), expected_feature_names)) {
    stop("input matrix feature order differs from source metadata")
  }
  if (!is.null(expected_cell_names) &&
    !identical(colnames(matrix), expected_cell_names)) {
    stop("input matrix cell order differs from source metadata")
  }

  counts <- methods::as(matrix, "dgCMatrix")
  if (!inherits(counts, "dgCMatrix")) {
    stop("failed to materialize a standard dgCMatrix")
  }
  if (!identical(as.numeric(dim(counts)), as.numeric(expected))) {
    stop("materialized sparse matrix dimensions changed")
  }
  if (!is.null(expected_feature_names) &&
    !identical(rownames(counts), expected_feature_names)) {
    stop("materialized sparse matrix feature order changed")
  }
  if (!is.null(expected_cell_names) &&
    !identical(colnames(counts), expected_cell_names)) {
    stop("materialized sparse matrix cell order changed")
  }
  if (anyNA(counts@x) || any(counts@x < 0)) {
    stop("materialized sparse count matrix contains invalid values")
  }
  if (any(counts@x != floor(counts@x))) {
    stop("materialized sparse count matrix contains non-integer values")
  }
  counts
}

verify_materialized_matrix_content <- function(
  source_matrix,
  materialized_matrix,
  block_cells = 5000L
) {
  if (!inherits(materialized_matrix, "dgCMatrix")) {
    stop("materialized matrix must be a dgCMatrix")
  }
  block_cells <- as.integer(block_cells)
  if (length(block_cells) != 1L || is.na(block_cells) ||
    block_cells < 1L) {
    stop("matrix verification block size must be a positive integer")
  }
  if (!identical(dim(source_matrix), dim(materialized_matrix)) ||
    !identical(dimnames(source_matrix), dimnames(materialized_matrix))) {
    stop("source and materialized matrix dimensions or names differ")
  }

  starts <- seq.int(1L, ncol(materialized_matrix), by = block_cells)
  for (start in starts) {
    end <- min(start + block_cells - 1L, ncol(materialized_matrix))
    columns <- seq.int(start, end)
    source_block <- methods::as(
      source_matrix[, columns, drop = FALSE],
      "dgCMatrix"
    )
    materialized_block <- materialized_matrix[
      , columns,
      drop = FALSE
    ]
    if (!inherits(materialized_block, "dgCMatrix")) {
      materialized_block <- methods::as(
        materialized_block,
        "dgCMatrix"
      )
    }
    identical_content <- identical(
      dimnames(source_block),
      dimnames(materialized_block)
    ) && identical(source_block@p, materialized_block@p) &&
      identical(source_block@i, materialized_block@i) &&
      identical(source_block@x, materialized_block@x)
    if (!identical_content) {
      stop(
        "source and materialized matrix differ in cell columns ",
        start,
        "-",
        end
      )
    }
  }
  list(
    verified = TRUE,
    block_cells = block_cells,
    blocks = length(starts)
  )
}

matrix_content_summary <- function(matrix) {
  data.frame(
    features = nrow(matrix),
    cells = ncol(matrix),
    nonzero_values = as.numeric(Matrix::nnzero(matrix)),
    total_counts = as.numeric(sum(matrix)),
    stringsAsFactors = FALSE
  )
}

create_processed_object <- function(
  counts,
  metadata,
  dataset,
  input_cells,
  input_features
) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("SeuratObject is required")
  }
  if (!inherits(counts, "dgCMatrix")) {
    stop("processed counts must be a standard dgCMatrix")
  }
  if (!identical(colnames(counts), metadata$Cells)) {
    stop("sparse matrix and metadata cell order differ")
  }
  if (!identical(colnames(counts), input_cells)) {
    stop("sparse matrix does not contain all requested cells in order")
  }
  if (!identical(rownames(counts), input_features)) {
    stop("sparse matrix does not contain all input features in order")
  }
  SeuratObject::CreateSeuratObject(
    counts = counts,
    meta.data = metadata,
    project = dataset
  )
}

create_layered_processed_object <- function(
  counts_layers,
  metadata,
  dataset,
  input_cells,
  input_features
) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("SeuratObject is required")
  }
  if (!is.list(counts_layers) ||
    length(counts_layers) < 2L ||
    is.null(names(counts_layers)) ||
    any(names(counts_layers) == "") ||
    anyDuplicated(names(counts_layers))) {
    stop("layered counts must be a named list with at least two layers")
  }
  valid_layers <- vapply(
    counts_layers,
    inherits,
    logical(1),
    what = "dgCMatrix"
  )
  if (any(!valid_layers)) {
    stop("every processed count layer must be a standard dgCMatrix")
  }
  if (any(vapply(
    counts_layers,
    function(layer) !identical(rownames(layer), input_features),
    logical(1)
  ))) {
    stop("layered sparse matrices do not share the input feature order")
  }
  layer_cells <- unlist(
    lapply(counts_layers, colnames),
    use.names = FALSE
  )
  if (anyDuplicated(layer_cells) ||
    !setequal(layer_cells, input_cells)) {
    stop("layered sparse matrices do not cover every input cell once")
  }
  if (!identical(rownames(metadata), metadata$Cells) ||
    !identical(metadata$Cells, input_cells)) {
    stop("layered metadata does not preserve the input cell order")
  }

  object <- SeuratObject::CreateSeuratObject(
    counts = counts_layers,
    meta.data = metadata,
    project = dataset
  )
  if (!setequal(colnames(object), input_cells)) {
    stop("layered Seurat object changed the input cell membership")
  }
  object
}

processed_count_layers <- function(object) {
  assay <- SeuratObject::DefaultAssay(object)
  layer_names <- SeuratObject::Layers(object[[assay]])
  count_layer_names <- layer_names[grepl(
    "^counts(?:[.]|$)",
    layer_names
  )]
  if (length(count_layer_names) == 0L) {
    stop("processed object has no count layer")
  }
  result <- lapply(
    count_layer_names,
    function(layer_name) {
      SeuratObject::LayerData(
        object,
        assay = assay,
        layer = layer_name
      )
    }
  )
  names(result) <- count_layer_names
  result
}

processed_counts <- function(object) {
  layers <- processed_count_layers(object)
  if (length(layers) != 1L) {
    stop(
      "processed object contains ",
      length(layers),
      " count layers; use processed_count_layers()"
    )
  }
  counts <- layers[[1L]]
  if (!inherits(counts, "dgCMatrix")) {
    stop(
      "processed object counts are not an in-memory standard dgCMatrix: ",
      paste(class(counts), collapse = "/")
    )
  }
  counts
}

processed_object_matrix_class <- function(object) {
  layers <- processed_count_layers(object)
  classes <- unique(vapply(
    layers,
    function(layer) class(layer)[[1L]],
    character(1)
  ))
  if (length(classes) != 1L || classes != "dgCMatrix") {
    return(paste(classes, collapse = ";"))
  }
  if (length(layers) == 1L) "dgCMatrix" else "dgCMatrix_layers"
}

validate_processed_object <- function(
  object,
  expected_features = NULL,
  expected_cells = NULL
) {
  if (!inherits(object, "Seurat")) {
    stop("processed RDS does not contain a Seurat object")
  }
  layers <- processed_count_layers(object)
  valid_layers <- vapply(
    layers,
    inherits,
    logical(1),
    what = "dgCMatrix"
  )
  if (any(!valid_layers)) {
    stop("processed object contains a non-dgCMatrix count layer")
  }
  if (!is.null(expected_features) && nrow(object) != expected_features) {
    stop("processed object feature count mismatch")
  }
  if (!is.null(expected_cells) && ncol(object) != expected_cells) {
    stop("processed object cell count mismatch")
  }
  if (!identical(colnames(object), rownames(object@meta.data))) {
    stop("processed object metadata and count matrix are not aligned")
  }
  if (any(vapply(
    layers,
    function(layer) {
      any(grepl(
        "BPCells|IterableMatrix|Delayed",
        class(layer),
        ignore.case = TRUE
      ))
    },
    logical(1)
  ))) {
    stop("processed object retains an on-disk or delayed matrix dependency")
  }
  object_features <- rownames(object)
  layer_features <- lapply(layers, rownames)
  if (any(vapply(layer_features, anyDuplicated, integer(1)) > 0L)) {
    stop("processed count layer contains duplicated features")
  }
  invalid_feature_order <- vapply(
    layer_features,
    function(features) {
      positions <- match(features, object_features)
      anyNA(positions) || is.unsorted(positions, strictly = TRUE)
    },
    logical(1)
  )
  if (any(invalid_feature_order)) {
    stop(
      "processed count layers are not ordered feature subsets of the object"
    )
  }
  if (!setequal(
    unique(unlist(layer_features, use.names = FALSE)),
    object_features
  )) {
    stop("processed count layers do not cover every object feature")
  }
  layer_cells <- unlist(lapply(layers, colnames), use.names = FALSE)
  if (anyDuplicated(layer_cells) ||
    !setequal(layer_cells, colnames(object))) {
    stop("processed count layers do not cover every object cell once")
  }
  invisible(TRUE)
}

processed_object_content_summary <- function(object) {
  validate_processed_object(object)
  layers <- processed_count_layers(object)
  data.frame(
    features = nrow(object),
    cells = ncol(object),
    count_layers = length(layers),
    nonzero_values = sum(vapply(
      layers,
      function(layer) as.numeric(Matrix::nnzero(layer)),
      numeric(1)
    )),
    total_counts = sum(vapply(
      layers,
      function(layer) as.numeric(sum(layer)),
      numeric(1)
    )),
    stringsAsFactors = FALSE
  )
}

processed_file_sha256 <- function(file) {
  command <- if (nzchar(Sys.which("sha256sum"))) {
    "sha256sum"
  } else if (nzchar(Sys.which("shasum"))) {
    "shasum"
  } else {
    stop("sha256sum or shasum is required")
  }
  args <- if (command == "shasum") {
    c("-a", "256", shQuote(file))
  } else {
    shQuote(file)
  }
  output <- system2(command, args, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop("SHA-256 calculation failed for ", file)
  }
  checksum <- strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
  if (!grepl("^[0-9a-f]{64}$", checksum)) {
    stop("invalid SHA-256 output for ", file)
  }
  checksum
}

refresh_processed_bundle_sha256 <- function(
  processed_dir,
  changed_files
) {
  manifest_file <- file.path(
    processed_dir,
    ".processed_bundle.sha256"
  )
  if (file.exists(manifest_file)) {
    lines <- readLines(manifest_file, warn = FALSE)
    valid <- grepl(
      "^[0-9a-f]{64}[[:space:]]+[*]?.+$",
      lines
    )
    if (length(lines) == 0L || any(!valid)) {
      stop("processed bundle checksum manifest is invalid")
    }
    checksums <- sub(
      "^([0-9a-f]{64})[[:space:]]+[*]?.+$",
      "\\1",
      lines
    )
    paths <- sub(
      "^[0-9a-f]{64}[[:space:]]+[*]?",
      "",
      lines
    )
  } else {
    checksums <- character()
    paths <- character()
  }
  changed_files <- unique(as.character(changed_files))
  all_paths <- unique(c(paths, changed_files))
  unsafe <- is.na(all_paths) |
    all_paths == "" |
    grepl("^/|(^|/)[.][.](/|$)", all_paths)
  if (any(unsafe) || anyDuplicated(paths)) {
    stop("processed bundle checksum manifest contains unsafe paths")
  }
  full_paths <- file.path(processed_dir, all_paths)
  if (any(!file.exists(full_paths))) {
    stop("processed bundle checksum refresh found a missing file")
  }

  checksum_by_path <- stats::setNames(checksums, paths)
  for (path in changed_files) {
    checksum_by_path[[path]] <- processed_file_sha256(file.path(
      processed_dir,
      path
    ))
  }
  ordered_paths <- sort(names(checksum_by_path))
  output <- paste(
    unname(checksum_by_path[ordered_paths]),
    ordered_paths,
    sep = "  "
  )
  temporary <- paste0(manifest_file, ".tmp.", Sys.getpid())
  on.exit(unlink(temporary), add = TRUE)
  writeLines(output, temporary, useBytes = TRUE)
  if (!file.rename(temporary, manifest_file)) {
    stop("failed to atomically refresh processed bundle checksums")
  }
  invisible(TRUE)
}

read_processed_shard_manifest <- function(processed_dir) {
  manifest_file <- file.path(
    processed_dir,
    "count_shard_manifest.tsv"
  )
  if (!file.exists(manifest_file)) {
    stop("processed count-shard manifest is missing: ", manifest_file)
  }
  manifest <- read.delim(
    manifest_file,
    sep = "\t",
    quote = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required <- c(
    "Dataset", "Shard_ID", "Source_Layer_Value", "RDS_File",
    "RDS_Bytes", "RDS_SHA256", "Cells", "Analysis_Cells",
    "Features", "Nonzero_Values", "Total_Counts", "Matrix_Class"
  )
  missing <- setdiff(required, names(manifest))
  if (length(missing) > 0L) {
    stop(
      "processed count-shard manifest lacks columns: ",
      paste(missing, collapse = ", ")
    )
  }
  if (nrow(manifest) == 0L ||
    anyDuplicated(manifest$Shard_ID) ||
    anyDuplicated(manifest$RDS_File) ||
    any(is.na(manifest$RDS_File)) ||
    any(grepl("^/|(^|/)[.][.](/|$)", manifest$RDS_File))) {
    stop("processed count-shard manifest contains unsafe shard paths")
  }
  manifest
}

validate_processed_shard_bundle <- function(
  processed_dir,
  expected_dataset = NULL,
  verify_sha256 = TRUE,
  verify_canonical_metadata = TRUE
) {
  manifest <- read_processed_shard_manifest(processed_dir)
  canonical_metadata <- NULL
  if (isTRUE(verify_canonical_metadata)) {
    canonical_file <- file.path(
      processed_dir,
      "metadata_canonical.tsv.gz"
    )
    if (!file.exists(canonical_file)) {
      stop("canonical metadata is missing for count-shard validation")
    }
    canonical_metadata <- if (
      requireNamespace("data.table", quietly = TRUE)
    ) {
      data.table::fread(canonical_file, data.table = FALSE)
    } else {
      read.delim(
        gzfile(canonical_file),
        sep = "\t",
        quote = "",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
    required_canonical <- c(
      "Cells", "Dataset", "Source_Cell_Index", "Analysis_Include"
    )
    missing_canonical <- setdiff(
      required_canonical,
      names(canonical_metadata)
    )
    if (length(missing_canonical) > 0L ||
      nrow(canonical_metadata) != sum(manifest$Cells) ||
      anyNA(canonical_metadata$Cells) ||
      anyDuplicated(canonical_metadata$Cells) ||
      !identical(
        as.integer(canonical_metadata$Source_Cell_Index),
        seq_len(nrow(canonical_metadata))
      )) {
      stop("canonical metadata is incomplete or nondeterministic")
    }
  }
  metadata_values_identical <- function(observed, expected) {
    if (is.factor(observed) || is.factor(expected) ||
      is.character(observed) || is.character(expected)) {
      return(identical(as.character(observed), as.character(expected)))
    }
    if (is.numeric(observed) && is.numeric(expected)) {
      observed <- as.numeric(observed)
      expected <- as.numeric(expected)
      if (!identical(is.na(observed), is.na(expected))) {
        return(FALSE)
      }
      retained <- !is.na(observed)
      return(isTRUE(all.equal(
        observed[retained],
        expected[retained],
        tolerance = sqrt(.Machine$double.eps),
        check.attributes = FALSE
      )))
    }
    identical(observed, expected)
  }
  source_indices <- vector("list", nrow(manifest))
  observed_total_counts <- numeric(nrow(manifest))

  for (index in seq_len(nrow(manifest))) {
    shard_file <- file.path(processed_dir, manifest$RDS_File[[index]])
    if (!file.exists(shard_file)) {
      stop("processed count shard is missing: ", shard_file)
    }
    if (as.numeric(file.info(shard_file)$size) !=
      as.numeric(manifest$RDS_Bytes[[index]])) {
      stop("processed count-shard byte size changed: ", shard_file)
    }
    if (isTRUE(verify_sha256) &&
      processed_file_sha256(shard_file) !=
        manifest$RDS_SHA256[[index]]) {
      stop("processed count-shard SHA-256 changed: ", shard_file)
    }

    object <- load_processed_object(shard_file)
    counts <- processed_counts(object)
    required_metadata <- c(
      "Dataset", "Count_Shard", "Source_Cell_Index",
      "Analysis_Include"
    )
    if (any(!required_metadata %in% names(object[[]]))) {
      stop("processed count shard lacks required metadata: ", shard_file)
    }
    if (ncol(object) != manifest$Cells[[index]] ||
      nrow(object) != manifest$Features[[index]] ||
      as.numeric(Matrix::nnzero(counts)) !=
        as.numeric(manifest$Nonzero_Values[[index]]) ||
      sum(object$Analysis_Include) !=
        manifest$Analysis_Cells[[index]] ||
      processed_object_matrix_class(object) !=
        manifest$Matrix_Class[[index]] ||
      any(object$Count_Shard != manifest$Shard_ID[[index]]) ||
      any(object$Dataset != manifest$Dataset[[index]])) {
      stop("processed count-shard contents changed: ", shard_file)
    }
    if (isTRUE(verify_canonical_metadata)) {
      expected_metadata <- canonical_metadata[
        as.integer(object$Source_Cell_Index), ,
        drop = FALSE
      ]
      observed_metadata <- object[[]]
      missing_embedded <- setdiff(
        names(expected_metadata),
        names(observed_metadata)
      )
      mismatched_columns <- if (length(missing_embedded) == 0L) {
        names(expected_metadata)[!vapply(
          names(expected_metadata),
          function(column) {
            metadata_values_identical(
              observed_metadata[[column]],
              expected_metadata[[column]]
            )
          },
          logical(1)
        )]
      } else {
        character()
      }
      if (length(missing_embedded) > 0L ||
        length(mismatched_columns) > 0L ||
        !identical(
          colnames(object),
          as.character(expected_metadata$Cells)
        )) {
        stop(
          "processed count-shard metadata differs from canonical sidecar: ",
          shard_file,
          if (length(missing_embedded) > 0L) {
            paste0(
              "; missing columns=",
              paste(missing_embedded, collapse = ",")
            )
          } else {
            ""
          },
          if (length(mismatched_columns) > 0L) {
            paste0(
              "; mismatched columns=",
              paste(mismatched_columns, collapse = ",")
            )
          } else {
            ""
          }
        )
      }
    }
    observed_total_counts[[index]] <- as.numeric(sum(counts))
    if (!isTRUE(all.equal(
      observed_total_counts[[index]],
      as.numeric(manifest$Total_Counts[[index]]),
      tolerance = 0
    ))) {
      stop("processed count-shard total counts changed: ", shard_file)
    }
    source_indices[[index]] <- as.integer(object$Source_Cell_Index)
    message(
      "[",
      index,
      "/",
      nrow(manifest),
      "] processed count shard and canonical metadata validated: ",
      basename(shard_file)
    )
    rm(object, counts)
    gc()
  }

  dataset_values <- unique(manifest$Dataset)
  if (length(dataset_values) != 1L ||
    (!is.null(expected_dataset) &&
      dataset_values != expected_dataset)) {
    stop("processed count-shard dataset identity is inconsistent")
  }
  combined_indices <- unlist(source_indices, use.names = FALSE)
  if (anyNA(combined_indices) ||
    anyDuplicated(combined_indices) ||
    !identical(sort(combined_indices), seq_along(combined_indices))) {
    stop("processed count shards do not cover every source cell once")
  }
  feature_values <- unique(manifest$Features)
  if (length(feature_values) != 1L) {
    stop("processed count shards do not share one feature dimension")
  }
  rm(canonical_metadata)
  gc()

  data.frame(
    Dataset = dataset_values,
    Shards = nrow(manifest),
    Full_Cells = sum(manifest$Cells),
    Analysis_Cells = sum(manifest$Analysis_Cells),
    Features = feature_values,
    Nonzero_Values = sum(manifest$Nonzero_Values),
    Total_Counts = sum(observed_total_counts),
    Matrix_Class = "dgCMatrix_shards",
    stringsAsFactors = FALSE
  )
}

refresh_processed_shard_metadata <- function(
  processed_dir,
  metadata,
  expected_dataset = NULL,
  verify_sha256 = TRUE
) {
  manifest <- read_processed_shard_manifest(processed_dir)
  required_metadata <- c(
    "Cells", "Dataset", "Source_Cell_Index", "Analysis_Include"
  )
  missing_metadata <- setdiff(required_metadata, names(metadata))
  if (length(missing_metadata) > 0L) {
    stop(
      "canonical metadata lacks shard-refresh columns: ",
      paste(missing_metadata, collapse = ", ")
    )
  }
  metadata$Cells <- as.character(metadata$Cells)
  if (anyNA(metadata$Cells) || anyDuplicated(metadata$Cells) ||
    anyNA(metadata$Source_Cell_Index) ||
    !identical(
      as.integer(metadata$Source_Cell_Index),
      seq_len(nrow(metadata))
    ) ||
    anyNA(metadata$Analysis_Include)) {
    stop("canonical metadata is not a complete deterministic cell table")
  }
  rownames(metadata) <- metadata$Cells
  dataset_values <- unique(as.character(metadata$Dataset))
  if (length(dataset_values) != 1L ||
    (!is.null(expected_dataset) &&
      !identical(dataset_values, expected_dataset)) ||
    any(as.character(manifest$Dataset) != dataset_values)) {
    stop("canonical metadata and shard manifest dataset identities differ")
  }

  shard_dir <- file.path(processed_dir, "count_shards")
  manifest_file <- file.path(processed_dir, "count_shard_manifest.tsv")
  stage_suffix <- paste0(
    ".metadata_stage.",
    format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"),
    ".",
    Sys.getpid()
  )
  stage_dir <- paste0(shard_dir, stage_suffix)
  stage_manifest <- paste0(manifest_file, stage_suffix)
  backup_dir <- paste0(shard_dir, ".metadata_previous.", Sys.getpid())
  backup_manifest <- paste0(
    manifest_file,
    ".metadata_previous.",
    Sys.getpid()
  )
  stale_paths <- c(
    stage_dir, stage_manifest, backup_dir, backup_manifest
  )
  if (any(file.exists(stale_paths) | dir.exists(stale_paths))) {
    stop("shard metadata refresh found a conflicting staging path")
  }
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  installed <- FALSE
  on.exit(
    {
      if (!installed) {
        if (dir.exists(stage_dir)) {
          unlink(stage_dir, recursive = TRUE)
        }
        if (file.exists(stage_manifest)) {
          unlink(stage_manifest)
        }
      }
    },
    add = TRUE
  )

  updated_manifest <- manifest
  covered_indices <- vector("list", nrow(manifest))
  for (index in seq_len(nrow(manifest))) {
    source_file <- file.path(
      processed_dir,
      manifest$RDS_File[[index]]
    )
    if (!file.exists(source_file) ||
      as.numeric(file.info(source_file)$size) !=
        as.numeric(manifest$RDS_Bytes[[index]]) ||
      (isTRUE(verify_sha256) &&
        processed_file_sha256(source_file) !=
          manifest$RDS_SHA256[[index]])) {
      stop("source count shard differs from its manifest: ", source_file)
    }
    object <- load_processed_object(source_file)
    counts <- processed_counts(object)
    source_indices <- as.integer(object$Source_Cell_Index)
    if (anyNA(source_indices) || anyDuplicated(source_indices) ||
      !identical(
        colnames(object),
        metadata$Cells[source_indices]
      )) {
      stop("count shard cell indices differ from canonical metadata")
    }
    shard_metadata <- metadata[source_indices, , drop = FALSE]
    shard_metadata$Count_Shard <- manifest$Shard_ID[[index]]
    rownames(shard_metadata) <- shard_metadata$Cells
    original_summary <- processed_object_content_summary(object)
    refreshed <- create_processed_object(
      counts = counts,
      metadata = shard_metadata,
      dataset = dataset_values,
      input_cells = shard_metadata$Cells,
      input_features = rownames(counts)
    )
    refreshed_summary <- processed_object_content_summary(refreshed)
    if (!identical(original_summary, refreshed_summary) ||
      sum(refreshed$Analysis_Include) !=
        sum(shard_metadata$Analysis_Include)) {
      stop("metadata refresh changed a count-shard matrix summary")
    }
    stage_file <- file.path(
      stage_dir,
      basename(manifest$RDS_File[[index]])
    )
    save_processed_object(refreshed, stage_file)
    updated_manifest$RDS_Bytes[[index]] <- as.numeric(
      file.info(stage_file)$size
    )
    updated_manifest$RDS_SHA256[[index]] <-
      processed_file_sha256(stage_file)
    updated_manifest$Analysis_Cells[[index]] <-
      sum(shard_metadata$Analysis_Include)
    updated_manifest$Excluded_Cells[[index]] <-
      sum(!shard_metadata$Analysis_Include)
    covered_indices[[index]] <- source_indices
    message(
      "[",
      index,
      "/",
      nrow(manifest),
      "] metadata-refreshed count shard validated: ",
      basename(stage_file)
    )
    rm(object, counts, shard_metadata, refreshed)
    gc()
  }
  combined_indices <- unlist(covered_indices, use.names = FALSE)
  if (anyNA(combined_indices) || anyDuplicated(combined_indices) ||
    !identical(sort(combined_indices), seq_len(nrow(metadata))) ||
    sum(updated_manifest$Cells) != nrow(metadata) ||
    sum(updated_manifest$Analysis_Cells) !=
      sum(metadata$Analysis_Include)) {
    stop("refreshed shards do not cover every canonical cell once")
  }
  write.table(
    updated_manifest,
    stage_manifest,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )

  if (!file.rename(shard_dir, backup_dir)) {
    stop("failed to preserve the previous count-shard directory")
  }
  if (!file.rename(stage_dir, shard_dir)) {
    file.rename(backup_dir, shard_dir)
    stop("failed to install metadata-refreshed count shards")
  }
  if (!file.rename(manifest_file, backup_manifest)) {
    file.rename(shard_dir, stage_dir)
    file.rename(backup_dir, shard_dir)
    stop("failed to preserve the previous count-shard manifest")
  }
  if (!file.rename(stage_manifest, manifest_file)) {
    file.rename(backup_manifest, manifest_file)
    file.rename(shard_dir, stage_dir)
    file.rename(backup_dir, shard_dir)
    stop("failed to install the metadata-refreshed shard manifest")
  }
  unlink(backup_dir, recursive = TRUE)
  unlink(backup_manifest)
  installed <- TRUE

  format_file <- file.path(processed_dir, "processed_object_format.tsv")
  refreshed_utc <- format(
    Sys.time(),
    "%Y-%m-%dT%H:%M:%SZ",
    tz = "UTC"
  )
  changed_files <- c(
    "count_shard_manifest.tsv",
    updated_manifest$RDS_File
  )
  if (file.exists(format_file)) {
    format_audit <- read.delim(
      format_file,
      sep = "\t",
      quote = "",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    if (nrow(format_audit) != 1L ||
      as.character(format_audit$Dataset[[1L]]) != dataset_values) {
      stop("processed-object format differs from refreshed shards")
    }
    format_audit$Full_Object_SHA256 <-
      processed_file_sha256(manifest_file)
    format_audit$Shard_Metadata_Refreshed_UTC <- refreshed_utc
    format_audit$Shard_Matrix_Recomputed <- FALSE
    temporary_format <- paste0(format_file, ".tmp.", Sys.getpid())
    write.table(
      format_audit,
      temporary_format,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      na = "NA"
    )
    if (!file.rename(temporary_format, format_file)) {
      stop("failed to update the processed-object format after shard refresh")
    }
    changed_files <- c(changed_files, "processed_object_format.tsv")
  }

  list(
    Dataset = dataset_values,
    Shards = nrow(updated_manifest),
    Full_Cells = sum(updated_manifest$Cells),
    Analysis_Cells = sum(updated_manifest$Analysis_Cells),
    Matrix_Recomputed = FALSE,
    Refreshed_UTC = refreshed_utc,
    Changed_Files = changed_files
  )
}

read_processed_object_format <- function(processed_dir) {
  format_file <- file.path(
    processed_dir,
    "processed_object_format.tsv"
  )
  if (!file.exists(format_file)) {
    stop("processed-object format manifest is missing: ", format_file)
  }
  format_audit <- read.delim(
    format_file,
    sep = "\t",
    quote = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required <- c(
    "Dataset", "Full_Object", "Analysis_Object",
    "Full_Cells", "Analysis_Cells", "Matrix_Class",
    "Self_Contained", "Technical_Batch_Design_Audit",
    "Source_Technical_Field_Inventory"
  )
  missing <- setdiff(required, names(format_audit))
  if (length(missing) > 0L) {
    stop(
      "processed-object format manifest lacks columns: ",
      paste(missing, collapse = ", ")
    )
  }
  if (!"Format_Type" %in% names(format_audit)) {
    format_audit$Format_Type <- "paired_full_analysis_objects"
  }
  audit_files <- unlist(
    format_audit[
      1L,
      c(
        "Technical_Batch_Design_Audit",
        "Source_Technical_Field_Inventory"
      )
    ],
    use.names = FALSE
  )
  if (anyNA(audit_files) || any(audit_files == "") || any(grepl(
    "^/|(^|/)[.][.](/|$)",
    audit_files
  ))) {
    stop("processed-object format has unsafe metadata-audit paths")
  }
  if (nrow(format_audit) != 1L ||
    !format_audit$Format_Type[[1L]] %in% c(
      "paired_full_analysis_objects",
      "single_retained_object_with_canonical_metadata",
      "single_full_object_without_analysis",
      "sharded_full_object_with_analysis_flag"
    )) {
    stop("processed-object format manifest is invalid or unsafe")
  }
  full_object <- format_audit$Full_Object[[1L]]
  if (is.na(full_object) ||
    full_object == "" ||
    grepl("^/|(^|/)[.][.](/|$)", full_object)) {
    stop("processed-object format manifest has an unsafe full object")
  }
  if (format_audit$Format_Type[[1L]] ==
    "paired_full_analysis_objects") {
    analysis_object <- format_audit$Analysis_Object[[1L]]
    if (is.na(analysis_object) ||
      analysis_object == "" ||
      grepl("^/|(^|/)[.][.](/|$)", analysis_object)) {
      stop(
        "paired processed-object format has an unsafe analysis object"
      )
    }
  } else if (format_audit$Format_Type[[1L]] %in% c(
    "single_retained_object_with_canonical_metadata",
    "single_full_object_without_analysis"
  )) {
    required_single <- c(
      "Canonical_Metadata", "Feature_Metadata",
      "Full_Object_SHA256",
      "Matrix_Recomputed"
    )
    missing_single <- setdiff(
      required_single,
      names(format_audit)
    )
    if (length(missing_single) > 0L) {
      stop(
        "single retained-object format lacks columns: ",
        paste(missing_single, collapse = ", ")
      )
    }
    if (format_audit$Format_Type[[1L]] ==
      "single_retained_object_with_canonical_metadata") {
      required_crosswalks <- c(
        "Age_Crosswalk", "Region_Crosswalk",
        "Source_Cell_Type_Crosswalk",
        "Technical_Batch_Crosswalk",
        "Biological_Unit_Crosswalk"
      )
      missing_crosswalks <- setdiff(
        required_crosswalks,
        names(format_audit)
      )
      if (length(missing_crosswalks) > 0L) {
        stop(
          "retained-object format lacks crosswalk columns: ",
          paste(missing_crosswalks, collapse = ", ")
        )
      }
      crosswalk_files <- unlist(
        format_audit[1L, required_crosswalks],
        use.names = FALSE
      )
      if (anyNA(crosswalk_files) ||
        any(crosswalk_files == "") ||
        any(grepl(
          "^/|(^|/)[.][.](/|$)",
          crosswalk_files
        ))) {
        stop("retained-object format has unsafe crosswalk paths")
      }
    }
    canonical_metadata <- format_audit$Canonical_Metadata[[1L]]
    feature_metadata <- format_audit$Feature_Metadata[[1L]]
    if (is.na(canonical_metadata) ||
      canonical_metadata == "" ||
      grepl(
        "^/|(^|/)[.][.](/|$)",
        canonical_metadata
      ) ||
      is.na(feature_metadata) ||
      feature_metadata == "" ||
      grepl(
        "^/|(^|/)[.][.](/|$)",
        feature_metadata
      ) ||
      (format_audit$Format_Type[[1L]] ==
        "single_retained_object_with_canonical_metadata" &&
        isTRUE(as.logical(
          format_audit$Matrix_Recomputed[[1L]]
        )))) {
      stop("single retained-object format manifest is invalid")
    }
  }
  format_audit
}

validate_metadata_audit_sidecars <- function(processed_dir, format_audit) {
  read_audit <- function(column) {
    file <- file.path(processed_dir, format_audit[[column]][[1L]])
    read.delim(
      file,
      sep = "\t",
      quote = "",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  dataset <- as.character(format_audit$Dataset[[1L]])
  full_cells <- as.numeric(format_audit$Full_Cells[[1L]])
  analysis_cells <- as.numeric(format_audit$Analysis_Cells[[1L]])
  design <- read_audit("Technical_Batch_Design_Audit")
  inventory <- read_audit("Source_Technical_Field_Inventory")
  required_design <- c(
    "Dataset", "Record_Type", "Cells", "Analysis_Cells",
    "Technical_Batch_Known_Cells",
    "Technical_Batch_Levels", "Technical_Batch_Design_Gate_Passed",
    "Analysis_Technical_Batch_Known_Cells",
    "Analysis_Technical_Batch_Levels",
    "All_Analysis_Batch_Levels_Span_Multiple_Donors_Or_Specimens",
    "Declared_Eligibility_Consistent_With_Design_Gate",
    "Technical_Batch_Derivation_Rule"
  )
  required_inventory <- c(
    "Dataset", "Record_Type", "Source_Metadata_Rows",
    "Candidate_Source_Fields", "Exact_Audited_Technical_Batch_Source"
  )
  missing_design <- setdiff(required_design, names(design))
  missing_inventory <- setdiff(required_inventory, names(inventory))
  if (length(missing_design) > 0L || length(missing_inventory) > 0L) {
    stop(
      "metadata audit sidecars lack required columns: ",
      paste(c(missing_design, missing_inventory), collapse = ", ")
    )
  }
  design_summary <- design[
    design$Record_Type == "dataset_summary", ,
    drop = FALSE
  ]
  inventory_summary <- inventory[
    inventory$Record_Type == "dataset_summary", ,
    drop = FALSE
  ]
  if (nrow(design_summary) != 1L || nrow(inventory_summary) != 1L ||
    !identical(as.character(design_summary$Dataset), dataset) ||
    !identical(as.character(inventory_summary$Dataset), dataset) ||
    as.numeric(design_summary$Cells[[1L]]) != full_cells ||
    as.numeric(design_summary$Analysis_Cells[[1L]]) != analysis_cells ||
    as.numeric(inventory_summary$Source_Metadata_Rows[[1L]]) !=
      full_cells) {
    stop("metadata audit summaries differ from the object manifest")
  }
  if (!isTRUE(as.logical(
    design_summary$Declared_Eligibility_Consistent_With_Design_Gate[[1L]]
  ))) {
    stop(
      "technical-batch eligibility is inconsistent with the verified ",
      "coverage and donor/specimen design gate"
    )
  }
  batch_levels <- design[
    design$Record_Type == "batch_level", ,
    drop = FALSE
  ]
  audit_count <- function(value, label) {
    value <- suppressWarnings(as.numeric(value))
    if (length(value) != 1L || is.na(value) || !is.finite(value) ||
      value < 0 || value != floor(value)) {
      stop("technical-batch audit has invalid ", label)
    }
    value
  }
  expected_levels <- audit_count(
    design_summary$Technical_Batch_Levels[[1L]],
    "dataset-level batch count"
  )
  expected_known_cells <- audit_count(
    design_summary$Technical_Batch_Known_Cells[[1L]],
    "dataset-level known-cell count"
  )
  expected_analysis_known_cells <- audit_count(
    design_summary$Analysis_Technical_Batch_Known_Cells[[1L]],
    "analysis known-cell count"
  )
  level_cells <- vapply(
    batch_levels$Cells,
    audit_count,
    numeric(1),
    label = "batch-level cell count"
  )
  level_analysis_cells <- vapply(
    batch_levels$Analysis_Cells,
    audit_count,
    numeric(1),
    label = "batch-level analysis-cell count"
  )
  if (nrow(batch_levels) != expected_levels ||
    sum(level_cells) != expected_known_cells ||
    sum(level_analysis_cells) != expected_analysis_known_cells) {
    stop("technical-batch level audit does not cover all known batch cells")
  }
  candidate_rows <- inventory[
    inventory$Record_Type == "candidate_field", ,
    drop = FALSE
  ]
  expected_candidates <- as.integer(
    inventory_summary$Candidate_Source_Fields[[1L]]
  )
  if (nrow(candidate_rows) != expected_candidates) {
    stop("source technical-field inventory candidate count is inconsistent")
  }
  exact_source <- as.logical(
    candidate_rows$Exact_Audited_Technical_Batch_Source
  )
  exact_source[is.na(exact_source)] <- FALSE
  if (expected_levels > 0L &&
    !any(exact_source) &&
    all(is.na(normalize_missing_metadata(
      design_summary$Technical_Batch_Derivation_Rule
    )))) {
    stop(
      "verified technical batches lack both an exact source field and a ",
      "recorded derivation rule"
    )
  }
  invisible(TRUE)
}

validate_processed_object_bundle <- function(
  processed_dir,
  expected_dataset = NULL
) {
  format_audit <- read_processed_object_format(processed_dir)
  dataset <- format_audit$Dataset[[1L]]
  if (!is.null(expected_dataset) &&
    !identical(dataset, expected_dataset)) {
    stop("processed-object format dataset identity is inconsistent")
  }
  if (!isTRUE(as.logical(format_audit$Self_Contained[[1L]]))) {
    stop("processed-object format is not self-contained")
  }

  audit_files <- file.path(
    processed_dir,
    unlist(
      format_audit[
        1L,
        c(
          "Technical_Batch_Design_Audit",
          "Source_Technical_Field_Inventory"
        )
      ],
      use.names = FALSE
    )
  )
  if (any(!file.exists(audit_files))) {
    stop("one or more metadata audit sidecars are missing")
  }
  validate_metadata_audit_sidecars(processed_dir, format_audit)

  format_type <- format_audit$Format_Type[[1L]]
  if (format_type == "sharded_full_object_with_analysis_flag") {
    return(validate_processed_shard_bundle(
      processed_dir,
      expected_dataset = expected_dataset,
      verify_sha256 = TRUE
    ))
  }
  full_file <- file.path(
    processed_dir,
    format_audit$Full_Object[[1L]]
  )
  if ("Full_Object_SHA256" %in% names(format_audit)) {
    recorded_sha256 <- as.character(
      format_audit$Full_Object_SHA256[[1L]]
    )
    if (!is.na(recorded_sha256) && nzchar(recorded_sha256)) {
      if (!grepl("^[0-9a-f]{64}$", recorded_sha256) ||
        processed_file_sha256(full_file) != recorded_sha256) {
        stop("processed full object SHA-256 changed")
      }
    }
  }
  full_object <- load_processed_object(full_file)
  if (format_type %in% c(
    "single_retained_object_with_canonical_metadata",
    "single_full_object_without_analysis"
  )) {
    canonical_file <- file.path(
      processed_dir,
      format_audit$Canonical_Metadata[[1L]]
    )
    feature_file <- file.path(
      processed_dir,
      format_audit$Feature_Metadata[[1L]]
    )
    if (!file.exists(canonical_file)) {
      stop("canonical metadata sidecar is missing: ", canonical_file)
    }
    if (!file.exists(feature_file)) {
      stop("feature metadata sidecar is missing: ", feature_file)
    }
    if (format_type ==
      "single_retained_object_with_canonical_metadata") {
      crosswalk_files <- file.path(
        processed_dir,
        unlist(
          format_audit[
            1L,
            c(
              "Age_Crosswalk", "Region_Crosswalk",
              "Source_Cell_Type_Crosswalk",
              "Technical_Batch_Crosswalk",
              "Biological_Unit_Crosswalk"
            )
          ],
          use.names = FALSE
        )
      )
      if (any(!file.exists(crosswalk_files))) {
        stop("one or more retained-object crosswalks are missing")
      }
    }
    connection <- if (grepl(
      "[.]gz$",
      canonical_file,
      ignore.case = TRUE
    )) {
      gzfile(canonical_file, "rt")
    } else {
      base::file(canonical_file, "rt")
    }
    canonical <- tryCatch(
      read.delim(
        connection,
        sep = "\t",
        quote = "",
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),
      finally = close(connection)
    )
    if (!all(c(
      "Cells", "Dataset", "Analysis_Include"
    ) %in% names(canonical)) ||
      !identical(as.character(canonical$Cells), colnames(full_object)) ||
      any(as.character(canonical$Dataset) != dataset)) {
      stop(
        "canonical metadata sidecar differs from the retained object"
      )
    }
    analysis_include <- as.logical(canonical$Analysis_Include)
    feature_connection <- if (grepl(
      "[.]gz$",
      feature_file,
      ignore.case = TRUE
    )) {
      gzfile(feature_file, "rt")
    } else {
      base::file(feature_file, "rt")
    }
    feature_metadata <- tryCatch(
      read.delim(
        feature_connection,
        sep = "\t",
        quote = "",
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),
      finally = close(feature_connection)
    )
    if (anyNA(analysis_include) ||
      nrow(canonical) != format_audit$Full_Cells[[1L]] ||
      sum(analysis_include) !=
        format_audit$Analysis_Cells[[1L]] ||
      ncol(full_object) != format_audit$Full_Cells[[1L]] ||
      !"Original_Feature_ID" %in% names(feature_metadata) ||
      !identical(
        as.character(feature_metadata$Original_Feature_ID),
        rownames(full_object)
      ) ||
      !identical(
        processed_object_matrix_class(full_object),
        format_audit$Matrix_Class[[1L]]
      )) {
      stop(
        "single retained-object contents differ from its manifest"
      )
    }
    result <- data.frame(
      Dataset = dataset,
      Objects = 1L,
      Full_Cells = ncol(full_object),
      Analysis_Cells = sum(analysis_include),
      Features = nrow(full_object),
      Matrix_Class = processed_object_matrix_class(full_object),
      stringsAsFactors = FALSE
    )
    rm(full_object, canonical, feature_metadata)
    gc()
    return(result)
  }

  if (!"Analysis_Include" %in% names(full_object[[]])) {
    stop("complete processed object lacks Analysis_Include")
  }
  full_dataset_values <- unique(as.character(full_object$Dataset))
  expected_analysis_cells <- colnames(full_object)[
    as.logical(full_object$Analysis_Include)
  ]
  expected_features <- rownames(full_object)
  full_cells <- ncol(full_object)
  full_features <- nrow(full_object)
  matrix_class <- processed_object_matrix_class(full_object)
  if (length(full_dataset_values) != 1L ||
    !identical(full_dataset_values, dataset) ||
    !identical(matrix_class, format_audit$Matrix_Class[[1L]]) ||
    full_cells != format_audit$Full_Cells[[1L]] ||
    length(expected_analysis_cells) !=
      format_audit$Analysis_Cells[[1L]]) {
    stop("complete processed object differs from its format manifest")
  }
  rm(full_object)
  gc()

  analysis_object <- load_processed_object(file.path(
    processed_dir,
    format_audit$Analysis_Object[[1L]]
  ))
  if (!"Analysis_Include" %in% names(analysis_object[[]])) {
    stop("analysis processed object lacks Analysis_Include")
  }
  analysis_dataset_values <- unique(as.character(
    analysis_object$Dataset
  ))
  if (length(analysis_dataset_values) != 1L ||
    !identical(analysis_dataset_values, dataset) ||
    !identical(
      processed_object_matrix_class(analysis_object),
      matrix_class
    ) ||
    !identical(colnames(analysis_object), expected_analysis_cells) ||
    !identical(rownames(analysis_object), expected_features) ||
    ncol(analysis_object) != format_audit$Analysis_Cells[[1L]] ||
    nrow(analysis_object) != full_features ||
    !all(as.logical(analysis_object$Analysis_Include))) {
    stop(
      "analysis processed object is not the manifest-defined ",
      "complete-object subset"
    )
  }

  data.frame(
    Dataset = dataset,
    Objects = 2L,
    Full_Cells = full_cells,
    Analysis_Cells = ncol(analysis_object),
    Features = nrow(analysis_object),
    Matrix_Class = matrix_class,
    stringsAsFactors = FALSE
  )
}

load_processed_shards <- function(
  processed_dir,
  analysis_only = FALSE,
  verify_sha256 = TRUE
) {
  validate_processed_shard_bundle(
    processed_dir,
    verify_sha256 = verify_sha256
  )
  manifest <- read_processed_shard_manifest(processed_dir)
  objects <- vector("list", nrow(manifest))
  names(objects) <- manifest$Shard_ID
  for (index in seq_len(nrow(manifest))) {
    object <- load_processed_object(file.path(
      processed_dir,
      manifest$RDS_File[[index]]
    ))
    if (isTRUE(analysis_only)) {
      object <- object[, object$Analysis_Include]
      validate_processed_object(object)
    }
    objects[[index]] <- object
  }
  objects
}

save_processed_object <- function(
  object,
  file,
  compression = Sys.getenv(
    "BRAINOMICS_RDS_COMPRESSION",
    unset = "gzip"
  ),
  validate_reload = TRUE
) {
  compression <- tolower(compression)
  compression_value <- switch(compression,
    gzip = "gzip",
    bzip2 = "bzip2",
    xz = "xz",
    none = FALSE,
    stop(
      "BRAINOMICS_RDS_COMPRESSION must be gzip, bzip2, xz, or none"
    )
  )

  validate_processed_object(object)
  original_summary <- if (isTRUE(validate_reload)) {
    processed_object_content_summary(object)
  } else {
    NULL
  }
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  temporary_file <- paste0(file, ".tmp.", Sys.getpid())
  on.exit(unlink(temporary_file), add = TRUE)
  if (identical(compression_value, "gzip")) {
    gzip_level <- suppressWarnings(as.integer(Sys.getenv(
      "BRAINOMICS_RDS_GZIP_LEVEL",
      unset = "1"
    )))
    if (is.na(gzip_level) || gzip_level < 0L || gzip_level > 9L) {
      stop("BRAINOMICS_RDS_GZIP_LEVEL must be an integer from 0 to 9")
    }
    connection <- gzfile(
      temporary_file,
      open = "wb",
      compression = gzip_level
    )
    tryCatch(
      serialize(object, connection),
      finally = close(connection)
    )
  } else {
    saveRDS(object, temporary_file, compress = compression_value)
  }

  if (isTRUE(validate_reload)) {
    restored <- readRDS(temporary_file)
    validate_processed_object(
      restored,
      expected_features = nrow(object),
      expected_cells = ncol(object)
    )
    if (!identical(
      processed_object_content_summary(restored),
      original_summary
    )) {
      stop("reloaded processed object content summary changed")
    }
    rm(restored)
    gc()
  }

  if (file.exists(file)) {
    backup <- paste0(
      file,
      ".previous.",
      format(Sys.time(), "%Y%m%dT%H%M%S")
    )
    if (!file.rename(file, backup)) {
      stop("failed to preserve the previous processed object")
    }
    if (!file.rename(temporary_file, file)) {
      file.rename(backup, file)
      stop("failed to install the new processed object")
    }
    unlink(backup)
  } else if (!file.rename(temporary_file, file)) {
    stop("failed to install the new processed object")
  }
  invisible(file)
}

load_processed_object <- function(file) {
  if (!file.exists(file)) {
    stop("processed object does not exist: ", file)
  }
  object <- readRDS(file)
  validate_processed_object(object)
  object
}
