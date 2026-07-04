`%||%` <- function(x, y) if (is.null(x)) y else x

log_message <- function(..., message_type = NULL, verbose = TRUE) {
  text <- paste0(..., collapse = "")
  text <- gsub("\\{\\.arg ([^}]*)\\}", "\\1", text)
  text <- gsub("\\{\\.cls ([^}]*)\\}", "\\1", text)
  text <- gsub("\\{\\.file ([^}]*)\\}", "\\1", text)
  text <- gsub("\\{\\.val ([^}]*)\\}", "\\1", text)
  text <- gsub("\\{\\.code ([^}]*)\\}", "\\1", text)
  if (identical(message_type, "error")) stop(text, call. = FALSE)
  if (isTRUE(verbose)) message(text)
  invisible(text)
}

GetAssayData5 <- function(srt, assay, layer = "counts") {
  SeuratObject::LayerData(srt, assay = assay, layer = layer)
}

GetFeaturesData <- function(srt, assay) {
  assay_obj <- srt[[assay]]
  meta <- tryCatch(assay_obj[[]], error = function(e) NULL)
  if (is.null(meta)) {
    data.frame(row.names = rownames(assay_obj))
  } else {
    as.data.frame(meta)
  }
}

align_sparse_features <- function(x, features) {
  if (!inherits(x, "dgCMatrix")) {
    x <- methods::as(x, "dgCMatrix")
  }
  if (identical(rownames(x), features)) {
    return(x)
  }

  row_map <- match(rownames(x), features)
  keep_rows <- !is.na(row_map)
  if (!all(keep_rows)) {
    x <- x[keep_rows, , drop = FALSE]
    row_map <- row_map[keep_rows]
  }
  if (length(row_map) == length(features) && identical(row_map, seq_along(features))) {
    rownames(x) <- features
    return(x)
  }

  col_counts <- diff(x@p)
  j <- rep.int(seq_len(ncol(x)), col_counts)
  Matrix::sparseMatrix(
    i = row_map[x@i + 1L],
    j = j,
    x = x@x,
    dims = c(length(features), ncol(x)),
    dimnames = list(features, colnames(x))
  )
}

build_sparse_x <- function(srt, assay, layer, features, verbose = TRUE) {
  sp <- reticulate::import("scipy.sparse", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)
  layer_names <- SeuratObject::Layers(srt[[assay]], search = layer)
  if (length(layer_names) == 0L) {
    log_message("Layer not found in assay ", assay, ": ", layer, message_type = "error")
  }

  cell_order <- character(0)
  x_parts <- vector("list", length(layer_names))
  for (i in seq_along(layer_names)) {
    layer_matrix <- align_sparse_features(
      GetAssayData5(srt, assay = assay, layer = layer_names[[i]]),
      features = features
    )
    cell_order <- c(cell_order, colnames(layer_matrix))
    x_part <- Matrix::t(layer_matrix)
    x_parts[[i]] <- sp$csr_matrix(reticulate::r_to_py(x_part, convert = FALSE), dtype = np$float32)
    rm(layer_matrix, x_part)
    gc()
  }

  if (length(x_parts) == 1L) {
    x <- x_parts[[1]]
  } else {
    log_message(
      "Stacking ", length(x_parts), " sparse layer matrices matching ", layer,
      " without joining Seurat layers in R.",
      verbose = verbose
    )
    x <- sp$vstack(x_parts, format = "csr", dtype = np$float32)
  }
  list(x = x, cells = cell_order)
}

srt_to_adata <- function(
  srt,
  features = NULL,
  assay_x = "RNA",
  layer_x = "counts",
  assay_y = c("spliced", "unspliced"),
  layer_y = "counts",
  reductions = NULL,
  graphs = NULL,
  neighbors = NULL,
  convert_tools = FALSE,
  convert_misc = FALSE,
  verbose = TRUE
) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} is not a {.cls Seurat}", message_type = "error")
  }
  log_message("Converting {.cls Seurat} to {.cls AnnData} ...", verbose = verbose)
  if (is.null(features)) features <- rownames(srt[[assay_x]])
  if (length(layer_y) == 1) {
    layer_y <- rep(layer_y, length(assay_y))
    names(layer_y) <- assay_y
  } else if (length(layer_y) != length(assay_y)) {
    log_message("{.arg layer_y} must be one character or the same length of the {.arg assay_y}", message_type = "error")
  }

  ad <- reticulate::import("anndata", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)
  sparse_x <- build_sparse_x(srt, assay = assay_x, layer = layer_x, features = features, verbose = verbose)
  X <- sparse_x$x
  cell_order <- sparse_x$cells

  obs <- srt@meta.data[cell_order, , drop = FALSE]
  if (ncol(obs) > 0) {
    for (i in seq_len(ncol(obs))) {
      if (is.logical(obs[, i])) obs[, i] <- factor(as.character(obs[, i]), levels = c("TRUE", "FALSE"))
    }
  }

  var <- GetFeaturesData(srt, assay = assay_x)[features, , drop = FALSE]
  if (ncol(var) > 0) {
    for (i in seq_len(ncol(var))) {
      if (is.logical(var[, i]) && !identical(colnames(var)[i], "highly_variable")) {
        var[, i] <- factor(as.character(var[, i]), levels = c("TRUE", "FALSE"))
      }
    }
  }
  var_features <- SeuratObject::VariableFeatures(srt, assay = assay_x)
  if (length(var_features) > 0) {
    if ("highly_variable" %in% colnames(var)) var <- var[, colnames(var) != "highly_variable", drop = FALSE]
    var[["highly_variable"]] <- features %in% var_features
  }

  adata <- ad$AnnData(
    X = X,
    obs = obs,
    var = cbind(data.frame(features = features), var)
  )
  adata$var_names <- features

  layer_list <- list()
  for (assay in names(srt@assays)[names(srt@assays) != assay_x]) {
    if (assay %in% assay_y) {
      assay_layer_names <- SeuratObject::Layers(srt[[assay]], search = layer_y[assay])
      if (length(assay_layer_names) == 1L) {
        layer <- Matrix::t(GetAssayData5(srt, assay = assay, layer = assay_layer_names[[1]]))
        layer <- layer[cell_order, colnames(X), drop = FALSE]
        layer_list[[assay]] <- layer
      } else {
        log_message(assay, " has split layers and is not converted as an AnnData layer", message_type = "warning", verbose = verbose)
      }
    } else {
      log_message(assay, " is in the srt object but not converted", message_type = "warning", verbose = verbose)
    }
  }
  if (length(layer_list) > 0) adata$layers <- layer_list

  reduction_names <- reductions %||% names(srt@reductions)
  reduction_names <- intersect(reduction_names, names(srt@reductions))
  reduction_list <- list()
  for (reduction in reduction_names) {
    reduction_list[[reduction]] <- srt[[reduction]]@cell.embeddings[cell_order, , drop = FALSE]
  }
  if (length(reduction_list) > 0) adata$obsm <- reduction_list

  if (!identical(graphs, character(0)) || !identical(neighbors, character(0))) {
    log_message("graphs/neighbors are not converted by the lightweight h5ad exporter", message_type = "warning", verbose = verbose)
  }

  uns_list <- list()
  if (isTRUE(convert_misc)) {
    for (nm in names(srt@misc)) if (nm != "") uns_list[[nm]] <- srt@misc[[nm]]
  } else {
    log_message("misc slot is not converted", message_type = "warning", verbose = verbose)
  }
  if (isTRUE(convert_tools)) {
    for (nm in names(srt@tools)) if (nm != "") uns_list[[nm]] <- srt@tools[[nm]]
  } else {
    log_message("tools slot is not converted", message_type = "warning", verbose = verbose)
  }
  if (length(uns_list) > 0) adata$uns <- uns_list
  log_message("Convert Seurat to AnnData object completed", message_type = "success", verbose = verbose)
  adata
}

srt_to_h5ad <- function(
  srt,
  path,
  features = NULL,
  assay_x = "RNA",
  layer_x = "counts",
  assay_y = c("spliced", "unspliced"),
  layer_y = "counts",
  reductions = NULL,
  graphs = NULL,
  neighbors = NULL,
  convert_tools = FALSE,
  convert_misc = FALSE,
  overwrite = FALSE,
  verbose = TRUE
) {
  path <- normalizePath(path.expand(path), mustWork = FALSE, winslash = "/")
  if (file.exists(path) && !isTRUE(overwrite)) {
    log_message(path, " already exists. Set overwrite = TRUE to overwrite", message_type = "error")
  }
  log_message("Converting Seurat to AnnData and writing to ", path, " ...", verbose = verbose)
  adata <- srt_to_adata(
    srt = srt, features = features, assay_x = assay_x, layer_x = layer_x,
    assay_y = assay_y, layer_y = layer_y, reductions = reductions,
    graphs = graphs, neighbors = neighbors, convert_tools = convert_tools,
    convert_misc = convert_misc, verbose = verbose
  )
  adata$write_h5ad(path)
  log_message("Successfully written to ", path, message_type = "success", verbose = verbose)
  invisible(path)
}
