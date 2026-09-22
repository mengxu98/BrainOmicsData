process_allen <- function() {
  matrix_file <- file.path(raw_dir, "matrix.csv")
  metadata_file <- file.path(raw_dir, "metadata.csv")
  shard_dir <- file.path(processed_dir, "source_h5ad_shards")
  manifest <- file.path(shard_dir, "manifest.tsv")
  required_source_files(c(matrix_file, metadata_file))
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required for the complete Allen metadata")
  }
  converter <- file.path(
    brainomics_repo_root(),
    "processing",
    "allen_csv_to_h5ad_shards.py"
  )
  if (overwrite || !file.exists(manifest)) {
    python <- "python3"
    status <- system2(
      python,
      c(
        shQuote(converter),
        "--matrix", shQuote(matrix_file),
        "--metadata", shQuote(metadata_file),
        "--output-dir", shQuote(shard_dir)
      )
    )
    if (status != 0L) stop("Allen source-to-shard conversion failed")
  }
  source_meta <- as.data.frame(
    data.table::fread(
      metadata_file,
      data.table = FALSE,
      check.names = FALSE,
      showProgress = TRUE
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (anyDuplicated(source_meta$sample_name)) {
    stop("Allen source metadata contains duplicated cell IDs")
  }
  shard_manifest <- read.delim(
    manifest,
    sep = "\t",
    quote = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required_columns <- c("donor", "file", "cells", "features", "sha256")
  if (!all(required_columns %in% names(shard_manifest)) ||
    anyDuplicated(shard_manifest$donor)) {
    stop("Allen shard manifest is invalid")
  }
  layers <- list()
  metadata_rows <- list()
  observed_cells <- character()
  for (manifest_index in seq_len(nrow(shard_manifest))) {
    record <- shard_manifest[manifest_index, , drop = FALSE]
    shard_file <- file.path(shard_dir, record$file)
    if (!identical(processed_file_sha256(shard_file), record$sha256)) {
      stop("Allen shard SHA-256 mismatch: ", record$file)
    }
    matrix <- open_h5ad_counts(shard_file, "X")
    counts <- materialize_sparse_counts(
      matrix,
      expected_features = record$features,
      expected_cells = record$cells
    )
    matrix_cells <- colnames(counts)
    index <- match(matrix_cells, source_meta$sample_name)
    if (anyNA(index) || anyDuplicated(matrix_cells)) {
      stop("Allen shard and metadata cell IDs are not one-to-one")
    }
    donor_meta <- source_meta[index, , drop = FALSE]
    donor <- as.character(donor_meta$external_donor_name_label)
    if (any(donor != record$donor)) {
      stop("Allen shard donor differs from source metadata: ", record$donor)
    }
    age <- unname(c(
      "H18.30.001" = "60 years",
      "H18.30.002" = "50 years"
    )[donor])
    if (anyNA(age)) stop("Allen source contains an unexpected donor")
    sex <- c(F = "Female", M = "Male")[as.character(donor_meta$donor_sex_label)]
    metadata <- base_source_metadata(
      matrix_cells,
      donor,
      donor,
      donor,
      donor,
      age,
      sex,
      "Primary motor cortex",
      "Neurologically unaffected",
      donor_meta$cluster_label
    )
    metadata$cluster_label <- donor_meta$cluster_label
    metadata$class_label <- donor_meta$class_label
    metadata$subclass_label <- donor_meta$subclass_label
    metadata$cell_type_accession_label <- donor_meta$cell_type_accession_label
    metadata$cell_type_alias_label <- donor_meta$cell_type_alias_label
    metadata$Source_Outlier_Call <- as.logical(donor_meta$outlier_call)
    metadata$Source_CellType_Annotation_Origin <- "original_study_author"
    metadata$Library_Chemistry <- "10x Chromium 3' v3"
    metadata$Sequencing_Platform <- "Illumina NovaSeq 6000"
    metadata$Analysis_Include <- !metadata$Source_Outlier_Call
    metadata$Analysis_Role <- ifelse(
      metadata$Analysis_Include,
      "reference",
      "reference_excluded"
    )
    metadata$Exclusion_Reason <- ifelse(
      metadata$Analysis_Include,
      NA_character_,
      "Allen source outlier_call"
    )
    metadata <- age_provenance(
      metadata,
      age,
      sub(" years$", "", age),
      "years",
      "postnatal age at death",
      "Allen Brain Map donor documentation and source metadata",
      "source donor age retained in years",
      "source donor documentation",
      FALSE
    )
    layers[[record$donor]] <- counts
    metadata_rows[[record$donor]] <- metadata
    observed_cells <- c(observed_cells, matrix_cells)
  }
  if (anyDuplicated(observed_cells) ||
    !setequal(observed_cells, source_meta$sample_name) ||
    length(observed_cells) != nrow(source_meta)) {
    stop("Allen shards do not cover every source cell exactly once")
  }
  layers <- align_layer_features(layers)
  metadata <- bind_rows_by_name(metadata_rows)
  create_complete_source_object(layers, metadata, dataset)
}

if (sys.nframe() == 0L) {
  PROCESSING_DATASET <- "AllenM1"
  source("processing/process_original_source.R")
}
