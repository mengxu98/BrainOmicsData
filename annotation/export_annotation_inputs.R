suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

source("functions/data_paths.R")
source("functions/dataset_metadata.R")
source("functions/processed_object.R")
source("functions/integration.R")


integration_dir <- brainomics_data_path("integration_25")
object_file <- file.path(integration_dir, "objects_integrated.rds")
metadata_file <- file.path(integration_dir, "metadata_filtered.rds")
annotation_dir <- file.path(integration_dir, "annotation")
reduction_dir <- file.path(annotation_dir, "reductions")
dir.create(reduction_dir, recursive = TRUE, showWarnings = FALSE)

# Existing fitted spaces can be validated without loading the expression object.
files <- file.path(reduction_dir, paste0(tolower(method_levels), "_latent.rds"))
if (all(file.exists(files))) {
  cells <- as.character(readRDS(metadata_file)$Cells)
  rows <- lapply(seq_along(method_levels), function(i) {
    embedding <- readRDS(files[i])
    stopifnot(
      is.matrix(embedding), identical(dim(embedding), c(2602031L, 50L)),
      identical(rownames(embedding), cells), all(is.finite(embedding))
    )
    data.frame(
      Method = method_levels[i], Cells = length(cells), Dimensions = 50L,
      Cell_Order_SHA256 = sha256_text_vector(cells), Finite = TRUE,
      File = files[i], File_SHA256 = processed_file_sha256(files[i])
    )
  })
  write_tsv(do.call(rbind, rows), file.path(annotation_dir, "annotation_reduction_audit.tsv"))
  message("Validated four existing fitted spaces without refitting or loading expression")
  quit(save = "no", status = 0L)
}

if (!file.exists(object_file) || !file.exists(metadata_file)) {
  stop("Final integration object or metadata is missing")
}

thisutils::log_message("[annotation-input-export] ", "Loading final four-method integration object")
object <- load_processed_object(object_file)
metadata <- readRDS(metadata_file)
if (!identical(colnames(object), as.character(metadata$Cells)) ||
  !identical(rownames(object[[]]), colnames(object))) {
  stop("Object, matrix and metadata cell order differ")
}

method_identity <- integration_method_identity()
method_identity <- method_identity[method_identity$Space == "latent", ]
method_identity <- method_identity[match(
  method_levels,
  method_identity$Display_Name
), , drop = FALSE]

cell_hash <- sha256_text_vector(colnames(object))
audit_rows <- vector("list", nrow(method_identity))
for (index in seq_len(nrow(method_identity))) {
  method <- method_identity$Display_Name[[index]]
  reduction <- method_identity$Reduction[[index]]
  thisutils::log_message("[annotation-input-export] ", paste("Exporting", method, "latent representation"))
  embedding <- Embeddings(object, reduction = reduction)
  if (!identical(rownames(embedding), colnames(object)) ||
    any(!is.finite(embedding))) {
    stop(method, " latent representation failed validation")
  }
  output_file <- file.path(
    reduction_dir,
    paste0(tolower(method), "_latent.rds")
  )
  temporary_file <- paste0(output_file, ".tmp.", Sys.getpid())
  saveRDS(embedding, temporary_file, compress = FALSE)
  if (!file.rename(temporary_file, output_file)) {
    unlink(temporary_file)
    stop("Could not publish ", output_file)
  }
  audit_rows[[index]] <- data.frame(
    Method = method,
    Reduction = reduction,
    Cells = nrow(embedding),
    Dimensions = ncol(embedding),
    Cell_Order_SHA256 = cell_hash,
    Finite = TRUE,
    File = output_file,
    File_SHA256 = processed_file_sha256(output_file),
    stringsAsFactors = FALSE
  )
  rm(embedding)
  gc()
}

features <- rownames(object)
if (length(features) == 0L || anyDuplicated(features)) {
  stop("RNA feature identifiers are empty or duplicated")
}
feature_audit <- data.frame(
  Feature_Order = seq_along(features),
  Feature = features,
  stringsAsFactors = FALSE
)
write_tsv(feature_audit, file.path(annotation_dir, "rna_features.tsv.gz"))
write_tsv(
  do.call(rbind, audit_rows),
  file.path(annotation_dir, "annotation_reduction_audit.tsv")
)

thisutils::log_message("[annotation-input-export] ", paste(
  "Exported",
  nrow(method_identity),
  "latent representations for",
  format(ncol(object), big.mark = ","),
  "cells and audited",
  format(length(features), big.mark = ","),
  "RNA features"
))
