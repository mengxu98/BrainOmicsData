write_complete_metadata <- function(metadata, file) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(
      metadata,
      file,
      sep = "\t",
      quote = FALSE,
      na = "NA",
      compress = "gzip"
    )
    return(invisible(file))
  }

  connection <- gzfile(file, "wt")
  tryCatch(
    write.table(
      metadata,
      connection,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      na = "NA"
    ),
    finally = close(connection)
  )
  invisible(file)
}

atomic_write_table <- function(x, file) {
  temporary <- paste0(file, ".tmp.", Sys.getpid())
  on.exit(unlink(temporary), add = TRUE)
  write.table(
    x,
    temporary,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
  if (file.exists(file)) {
    unlink(file)
  }
  if (!file.rename(temporary, file)) {
    stop("failed to install table: ", file)
  }
  invisible(file)
}

process_layered_h5ad <- function(
  dataset,
  config,
  counts_h5,
  metadata,
  feature_meta,
  input_cells,
  input_features,
  processed_dir,
  canonical_metadata_file,
  processing_audit_file,
  object_format_file,
  full_object_file,
  analysis_object_file,
  external_metadata_audit = NULL,
  source_record = NULL,
  source_meta = NULL
) {
  if (!identical(config$matrix_storage, "dgCMatrix_shards")) {
    stop("sharded preprocessing requires matrix_storage=dgCMatrix_shards")
  }
  layer_by <- config$layer_by
  if (is.null(layer_by) || !layer_by %in% names(metadata)) {
    stop("sharded preprocessing column is missing: ", layer_by)
  }
  layer_values <- normalize_missing_metadata(metadata[[layer_by]])
  if (anyNA(layer_values)) {
    stop("sharded preprocessing column contains missing values")
  }
  source_layer_values <- unique(layer_values)
  if (!is.null(config$expected_layers) &&
    length(source_layer_values) != config$expected_layers) {
    stop(
      "shard count differs from the verified source: expected ",
      config$expected_layers,
      ", observed ",
      length(source_layer_values)
    )
  }

  metadata$Source_Cell_Index <- seq_len(nrow(metadata))
  validate_dataset_metadata(
    metadata,
    matrix_cells = input_cells,
    expected_cells = length(input_cells)
  )
  analysis_metadata <- metadata[
    metadata$Analysis_Include, ,
    drop = FALSE
  ]
  analysis_donors <- length(unique(stats::na.omit(
    analysis_metadata$Global_Donor_ID
  )))
  if (nrow(analysis_metadata) != config$expected_analysis_cells ||
    analysis_donors != config$expected_analysis_donors) {
    stop("sharded analysis inclusion counts differ from the config")
  }

  thisutils::log_message(
    paste0(
      "Writing ",
      length(source_layer_values),
      " self-contained dgCMatrix shards because the complete matrix ",
      "exceeds the 32-bit nonzero-index limit of one dgCMatrix"
    )
  )

  shard_dir <- file.path(processed_dir, "count_shards")
  stage_dir <- file.path(
    processed_dir,
    paste0("count_shards.incoming.", Sys.getpid())
  )
  if (dir.exists(stage_dir)) {
    stop("count-shard staging directory already exists: ", stage_dir)
  }
  dir.create(stage_dir, recursive = TRUE)
  stage_installed <- FALSE
  backup_dir <- NA_character_
  manifest_file <- file.path(
    processed_dir,
    "count_shard_manifest.tsv"
  )
  manifest_backup <- NA_character_
  install_complete <- FALSE
  cleanup_install <- function() {
    if (dir.exists(stage_dir)) {
      unlink(stage_dir, recursive = TRUE)
    }
    if (isTRUE(stage_installed) && !isTRUE(install_complete)) {
      if (dir.exists(shard_dir)) {
        unlink(shard_dir, recursive = TRUE)
      }
      if (!is.na(backup_dir) && dir.exists(backup_dir)) {
        file.rename(backup_dir, shard_dir)
      }
      if (file.exists(manifest_file)) {
        unlink(manifest_file)
      }
      if (!is.na(manifest_backup) &&
        file.exists(manifest_backup)) {
        file.rename(manifest_backup, manifest_file)
      }
    }
  }
  on.exit(cleanup_install(), add = TRUE)

  source_matrix_checksum <- BPCells::checksum(counts_h5)
  shard_rows <- vector("list", length(source_layer_values))
  source_total_counts <- 0
  materialized_total_counts <- 0

  for (layer_index in seq_along(source_layer_values)) {
    source_value <- source_layer_values[[layer_index]]
    selected <- which(layer_values == source_value)
    shard_id <- sprintf("batch_%02d", layer_index)
    shard_file_name <- paste0(dataset, "_", shard_id, ".rds")
    stage_file <- file.path(stage_dir, shard_file_name)
    thisutils::log_message(
      paste0(
        "Materializing ",
        shard_id,
        " (source ",
        source_value,
        "): ",
        length(selected),
        " cells"
      )
    )

    source_shard <- counts_h5[, selected]
    source_checksum <- BPCells::checksum(source_shard)
    source_shard_total <- sum(BPCells::colSums(source_shard))
    materialized <- materialize_sparse_counts(
      matrix = source_shard,
      expected_features = length(input_features),
      expected_cells = length(selected),
      expected_feature_names = input_features,
      expected_cell_names = input_cells[selected]
    )
    materialized_iterable <- BPCells::write_matrix_memory(
      BPCells::convert_matrix_type(materialized, "uint32_t")
    )
    materialized_checksum <- BPCells::checksum(
      materialized_iterable
    )
    materialized_shard_total <- sum(
      Matrix::colSums(materialized)
    )
    if (!identical(source_checksum, materialized_checksum) ||
      !isTRUE(all.equal(
        as.numeric(source_shard_total),
        as.numeric(materialized_shard_total),
        tolerance = 0
      ))) {
      stop("H5AD-to-dgCMatrix verification failed for ", shard_id)
    }

    shard_metadata <- metadata[selected, , drop = FALSE]
    shard_metadata$Count_Shard <- shard_id
    shard_object <- create_processed_object(
      counts = materialized,
      metadata = shard_metadata,
      dataset = dataset,
      input_cells = input_cells[selected],
      input_features = input_features
    )
    save_processed_object(shard_object, stage_file)

    shard_rows[[layer_index]] <- data.frame(
      Dataset = dataset,
      Shard_ID = shard_id,
      Layer_By = layer_by,
      Source_Layer_Value = source_value,
      Source_First_Cell_Index = min(selected),
      Source_Last_Cell_Index = max(selected),
      RDS_File = file.path("count_shards", shard_file_name),
      RDS_Bytes = as.numeric(file.info(stage_file)$size),
      RDS_SHA256 = processed_file_sha256(stage_file),
      Cells = length(selected),
      Analysis_Cells = sum(shard_metadata$Analysis_Include),
      Excluded_Cells = sum(!shard_metadata$Analysis_Include),
      Features = nrow(materialized),
      Nonzero_Values = as.numeric(Matrix::nnzero(materialized)),
      Total_Counts = as.numeric(materialized_shard_total),
      Source_Matrix_Checksum = source_checksum,
      Materialized_Matrix_Checksum = materialized_checksum,
      Source_Total_Counts = as.numeric(source_shard_total),
      Materialized_Total_Counts =
        as.numeric(materialized_shard_total),
      Matrix_Class = class(materialized)[[1L]],
      stringsAsFactors = FALSE
    )
    source_total_counts <- source_total_counts + source_shard_total
    materialized_total_counts <-
      materialized_total_counts + materialized_shard_total

    rm(
      source_shard,
      materialized,
      materialized_iterable,
      shard_metadata,
      shard_object
    )
    gc()
  }

  manifest <- do.call(rbind, shard_rows)
  if (sum(manifest$Cells) != length(input_cells) ||
    sum(manifest$Analysis_Cells) !=
      config$expected_analysis_cells ||
    any(manifest$Source_Matrix_Checksum !=
      manifest$Materialized_Matrix_Checksum) ||
    !isTRUE(all.equal(
      as.numeric(source_total_counts),
      as.numeric(materialized_total_counts),
      tolerance = 0
    ))) {
    stop("complete count-shard verification failed")
  }

  install_timestamp <- format(
    Sys.time(),
    "%Y%m%dT%H%M%S",
    tz = "UTC"
  )
  if (dir.exists(shard_dir)) {
    backup_dir <- paste0(
      shard_dir,
      ".previous.",
      install_timestamp
    )
    if (!file.rename(shard_dir, backup_dir)) {
      stop("failed to preserve the previous count-shard directory")
    }
  }
  if (file.exists(manifest_file)) {
    manifest_backup <- paste0(
      manifest_file,
      ".previous.",
      install_timestamp
    )
    if (!file.rename(manifest_file, manifest_backup)) {
      if (!is.na(backup_dir) && dir.exists(backup_dir)) {
        file.rename(backup_dir, shard_dir)
      }
      stop("failed to preserve the previous count-shard manifest")
    }
  }
  if (!file.rename(stage_dir, shard_dir)) {
    if (!is.na(backup_dir) && dir.exists(backup_dir)) {
      file.rename(backup_dir, shard_dir)
    }
    if (!is.na(manifest_backup) && file.exists(manifest_backup)) {
      file.rename(manifest_backup, manifest_file)
    }
    stop("failed to install the new count-shard directory")
  }
  stage_installed <- TRUE
  atomic_write_table(manifest, manifest_file)

  bundle_summary <- validate_processed_shard_bundle(
    processed_dir = processed_dir,
    expected_dataset = dataset,
    verify_sha256 = TRUE,
    verify_canonical_metadata = FALSE
  )
  if (bundle_summary$Full_Cells != config$expected_cells ||
    bundle_summary$Analysis_Cells !=
      config$expected_analysis_cells ||
    bundle_summary$Features != config$expected_features ||
    bundle_summary$Shards != config$expected_layers ||
    !isTRUE(all.equal(
      as.numeric(bundle_summary$Total_Counts),
      as.numeric(materialized_total_counts),
      tolerance = 0
    ))) {
    stop("reloaded count-shard bundle differs from the source")
  }

  write_complete_metadata(metadata, canonical_metadata_file)
  write_canonical_metadata_crosswalks(
    metadata,
    processed_dir,
    source_meta = source_meta
  )
  audit <- dataset_metadata_audit(metadata)
  if (is.null(source_record) || nrow(source_record) != 1L) {
    stop("sharded preprocessing requires one validated source record")
  }
  audit$source_accession <- source_record$source_accession[[1L]]
  audit$source_repository <- source_record$source_repository[[1L]]
  audit$source_publication_doi <- source_record$publication_doi[[1L]]
  audit$source_publication_title <- source_record$verified_title[[1L]]
  audit$source_publication_journal <- source_record$journal[[1L]]
  audit$source_publication_year <- source_record$publication_year[[1L]]
  audit$source_repository_record_url <-
    source_record$repository_record_url[[1L]]
  audit$input_features <- length(input_features)
  audit$full_object_cells <- config$expected_cells
  audit$full_object_features <- config$expected_features
  audit$analysis_object_cells <- config$expected_analysis_cells
  audit$analysis_object_features <- config$expected_features
  audit$analysis_donors <- analysis_donors
  audit$expected_analysis_cells <- config$expected_analysis_cells
  audit$expected_analysis_donors <- config$expected_analysis_donors
  audit$matrix_group <- config$matrix_group
  audit$expected_input_cells <- config$expected_cells
  audit$expected_input_features <- config$expected_features
  audit$all_input_cells_processed <- TRUE
  audit$all_input_features_processed <- TRUE
  audit$source_matrix_checksum <- source_matrix_checksum
  audit$materialized_matrix_checksum <-
    "verified independently for every count shard"
  audit$source_total_counts <- as.numeric(source_total_counts)
  audit$materialized_total_counts <-
    as.numeric(materialized_total_counts)
  audit$analysis_source_checksum <- NA_character_
  audit$analysis_written_checksum <- NA_character_
  audit$analysis_source_total_counts <- NA_real_
  audit$analysis_written_total_counts <- NA_real_
  audit$matrix_storage <- "dgCMatrix_shards"
  audit$count_shards <- nrow(manifest)
  audit$processed_object_self_contained <- TRUE
  audit$preprocessing_reader <- "BPCells H5AD streaming reader"
  audit$count_shard_manifest_sha256 <-
    processed_file_sha256(manifest_file)
  audit$analysis_access <- paste(
    "load shards with analysis_only=TRUE;",
    "no duplicate analysis matrix stored"
  )
  if (!is.null(external_metadata_audit)) {
    for (column in names(external_metadata_audit)) {
      audit[[paste0("external_", column)]] <-
        external_metadata_audit[[column]]
    }
  }
  atomic_write_table(audit, processing_audit_file)
  writeLines(
    capture.output(sessionInfo()),
    file.path(processed_dir, "sessionInfo.txt")
  )

  format_audit <- data.frame(
    Dataset = dataset,
    Format_Type = "sharded_full_object_with_analysis_flag",
    Full_Object = basename(manifest_file),
    Analysis_Object = NA_character_,
    Full_Cells = config$expected_cells,
    Analysis_Cells = config$expected_analysis_cells,
    Canonical_Metadata = basename(canonical_metadata_file),
    Feature_Metadata = "features_raw.tsv.gz",
    Age_Crosswalk = "age_crosswalk.tsv",
    Region_Crosswalk = "region_crosswalk.tsv",
    Source_Cell_Type_Crosswalk = "source_cell_type_crosswalk.tsv",
    Technical_Batch_Crosswalk = "technical_batch_crosswalk.tsv",
    Technical_Batch_Design_Audit =
      "technical_batch_design_audit.tsv",
    Source_Technical_Field_Inventory =
      "source_technical_field_inventory.tsv",
    Biological_Unit_Crosswalk =
      "donor_specimen_library_crosswalk.tsv.gz",
    Full_Object_SHA256 = processed_file_sha256(manifest_file),
    Matrix_Class = "dgCMatrix_shards",
    Count_Layers = nrow(manifest),
    Self_Contained = TRUE,
    Matrix_Recomputed = TRUE,
    Source_Accession = source_record$source_accession[[1L]],
    Source_Repository = source_record$source_repository[[1L]],
    Source_Publication_DOI = source_record$publication_doi[[1L]],
    Source_Publication_Title = source_record$verified_title[[1L]],
    Source_Publication_Journal = source_record$journal[[1L]],
    Source_Publication_Year = source_record$publication_year[[1L]],
    Source_Repository_Record_URL =
      source_record$repository_record_url[[1L]],
    Metadata_Schema_Version = brainomics_metadata_schema_version(),
    Completed_UTC = format(
      Sys.time(),
      "%Y-%m-%dT%H:%M:%SZ",
      tz = "UTC"
    ),
    stringsAsFactors = FALSE
  )
  atomic_write_table(format_audit, object_format_file)

  install_complete <- TRUE
  if (!is.na(backup_dir) && dir.exists(backup_dir)) {
    unlink(backup_dir, recursive = TRUE)
  }
  if (!is.na(manifest_backup) && file.exists(manifest_backup)) {
    unlink(manifest_backup)
  }
  for (obsolete_object in c(full_object_file, analysis_object_file)) {
    if (file.exists(obsolete_object)) {
      unlink(obsolete_object)
    }
  }
  thisutils::log_message(
    paste0(
      dataset,
      " full sharded processing completed: ",
      config$expected_cells,
      " cells, ",
      config$expected_features,
      " features, ",
      nrow(manifest),
      " self-contained dgCMatrix shards; analysis-eligible cells: ",
      config$expected_analysis_cells
    ),
    message_type = "success"
  )
  invisible(manifest_file)
}
