source("functions/data_paths.R")
source("functions/metadata_schema.R")
source("functions/dataset_metadata.R")
source("functions/dataset_config.R")
source("functions/processed_object.R")

required_packages <- c("BPCells", "Matrix", "SeuratObject")
missing_packages <- required_packages[!vapply(
  required_packages,
  requireNamespace,
  logical(1),
  quietly = TRUE
)]
if (length(missing_packages) > 0L) {
  stop(
    "Missing preprocessing packages: ",
    paste(missing_packages, collapse = ", ")
  )
}


args <- commandArgs(trailingOnly = TRUE)
dataset <- if (exists("PROCESSING_DATASET", inherits = FALSE)) {
  get("PROCESSING_DATASET", inherits = FALSE)
} else if (length(args) >= 1L) {
  args[[1L]]
} else {
  stop("Usage: Rscript processing/process_h5ad.R <dataset> [overwrite]")
}
overwrite_arg <- if (exists("PROCESSING_DATASET", inherits = FALSE)) {
  if (length(args) >= 1L) args[[1L]] else "false"
} else {
  if (length(args) >= 2L) args[[2L]] else "false"
}
overwrite <- tolower(overwrite_arg) %in% c("true", "t", "1")
metadata_only <- "--metadata-only" %in% args
config <- dataset_config(dataset)

repo_dir <- brainomics_repo_root()
source_record <- formal_source_access_record(
  dataset,
  file.path(repo_dir, "data", "source_access_summary.tsv")
)
raw_dir <- brainomics_data_path("raw", dataset)
processed_dir <- brainomics_data_path("processed", dataset)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

input_root <- if (is.null(config$input_root)) "raw" else config$input_root
h5ad_file <- brainomics_data_path(
  input_root,
  dataset,
  config$input_file
)
obs_file <- file.path(processed_dir, "metadata_raw.tsv.gz")
var_file <- file.path(processed_dir, "features_raw.tsv.gz")
export_audit_file <- file.path(
  processed_dir,
  "metadata_export_audit.json"
)
canonical_metadata_file <- file.path(
  processed_dir,
  "metadata_canonical.tsv.gz"
)
processing_audit_file <- file.path(
  processed_dir,
  "processing_audit.tsv"
)
metadata_standardization_audit_file <- file.path(
  processed_dir,
  "metadata_standardization_audit.tsv"
)
metadata_standardization_session_file <- file.path(
  processed_dir,
  "metadata_standardization_sessionInfo.txt"
)
object_format_file <- file.path(
  processed_dir,
  "processed_object_format.tsv"
)
metadata_standardization_summary_file <- file.path(
  processed_dir,
  "metadata_standardization_summary.tsv"
)
full_object_file <- file.path(
  processed_dir,
  paste0(dataset, "_full_processed.rds")
)
analysis_object_file <- file.path(
  processed_dir,
  paste0(dataset, "_processed.rds")
)

write_current_dataset_summary <- function(
  objects,
  full_cells,
  analysis_cells,
  features,
  matrix_class,
  status
) {
  summary <- data.frame(
    Dataset = dataset,
    Objects = as.integer(objects),
    Full_Cells = as.numeric(full_cells),
    Analysis_Cells = as.numeric(analysis_cells),
    Features = as.numeric(features),
    Matrix_Class = as.character(matrix_class),
    Status = as.character(status),
    stringsAsFactors = FALSE
  )
  write_metadata_table_atomic(
    summary,
    metadata_standardization_summary_file
  )
}

metadata_sidecars_complete <- all(file.exists(c(
  obs_file,
  var_file,
  export_audit_file
)))
metadata_export_required <- overwrite || !metadata_sidecars_complete
h5ad_required <- !metadata_only || metadata_export_required
if (h5ad_required && !file.exists(h5ad_file)) {
  stop("Missing H5AD input: ", h5ad_file)
}
if (h5ad_required && !is.null(config$expected_input_bytes)) {
  observed_bytes <- file.info(h5ad_file)$size
  if (!identical(
    as.numeric(observed_bytes),
    as.numeric(config$expected_input_bytes)
  )) {
    stop(
      "H5AD byte size differs from the verified source record: expected ",
      config$expected_input_bytes,
      ", observed ",
      observed_bytes
    )
  }
}
input_sha256 <- NA_character_
if (h5ad_required && !is.null(config$expected_input_sha256)) {
  input_sha256 <- processed_file_sha256(h5ad_file)
  if (!identical(input_sha256, config$expected_input_sha256)) {
    stop("H5AD SHA-256 differs from the verified source record")
  }
}

export_script <- file.path(repo_dir, "processing/export_h5ad_metadata.py")
if (metadata_export_required) {
  python <- "python3"
  status <- system2(
    python,
    c(
      shQuote(export_script),
      "--input", h5ad_file,
      "--obs-output", obs_file,
      "--var-output", var_file,
      "--audit-output", export_audit_file,
      "--matrix-group", config$matrix_group
    )
  )
  if (status != 0L) {
    stop("complete H5AD metadata export failed for ", dataset)
  }
}

read_complete_table <- function(file, sep = "\t") {
  if (grepl("[.]xlsx?$", file, ignore.case = TRUE)) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("readxl is required to read source Excel metadata: ", file)
    }
    return(as.data.frame(
      readxl::read_excel(file),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }
  read_brainomics_table(file, sep = sep)
}

raw_meta <- read_complete_table(obs_file)
feature_meta <- read_complete_table(var_file)
external_metadata_audit <- NULL
if (!is.null(config$external_metadata)) {
  external_configs <- if (!is.null(config$external_metadata$file)) {
    list(primary = config$external_metadata)
  } else {
    config$external_metadata
  }
  if (!is.list(external_configs) || length(external_configs) == 0L ||
    any(names(external_configs) == "") || anyNA(names(external_configs))) {
    stop("external metadata configuration must be a named non-empty list")
  }
  external_audits <- vector("list", length(external_configs))
  names(external_audits) <- names(external_configs)
  for (external_name in names(external_configs)) {
    external_config <- external_configs[[external_name]]
    external_file <- brainomics_data_path(
      external_config$root,
      dataset,
      external_config$file
    )
    if (!file.exists(external_file)) {
      stop("Missing external metadata input: ", external_file)
    }
    external_sep <- if (is.null(external_config$sep)) {
      "\t"
    } else {
      external_config$sep
    }
    external_meta <- read_complete_table(
      external_file,
      sep = external_sep
    )
    external_meta <- transform_external_metadata(
      external_meta,
      external_config
    )
    raw_meta <- join_external_metadata(
      matrix_meta = raw_meta,
      external_meta = external_meta,
      key_column = external_config$key_column,
      matrix_key_column = if (is.null(
        external_config$matrix_key_column
      )) {
        "Original_Cell_ID"
      } else {
        external_config$matrix_key_column
      },
      match_column = if (is.null(external_config$match_column)) {
        "External_Metadata_Matched"
      } else {
        external_config$match_column
      },
      require_all_matrix_cells =
        external_config$require_all_matrix_cells,
      require_all_metadata_cells =
        external_config$require_all_metadata_cells
    )
    external_audits[[external_name]] <- attr(
      raw_meta,
      "external_metadata_audit"
    )
  }
  if (length(external_audits) == 1L) {
    external_metadata_audit <- external_audits[[1L]]
  } else {
    external_metadata_audit <- do.call(
      cbind,
      lapply(names(external_audits), function(external_name) {
        audit <- external_audits[[external_name]]
        names(audit) <- paste0(external_name, "_", names(audit))
        audit
      })
    )
  }
}
raw_meta <- transform_raw_metadata(raw_meta, config)
metadata <- build_dataset_metadata(
  raw_meta = raw_meta,
  dataset = dataset,
  field_map = config$field_map,
  constants = config$constants,
  default_analysis_role = config$default_analysis_role,
  default_analysis_include = if (is.null(
    config$default_analysis_include
  )) {
    TRUE
  } else {
    config$default_analysis_include
  },
  default_exclusion_reason = if (is.null(
    config$default_exclusion_reason
  )) {
    NA_character_
  } else {
    config$default_exclusion_reason
  }
)
metadata <- add_source_publication_metadata(metadata, source_record)
metadata <- apply_donor_crosswalk(
  metadata,
  file.path(repo_dir, "data/donor_crosswalk.tsv")
)
metadata <- apply_required_filters(metadata, config)
metadata <- add_age_schema(metadata)
metadata <- add_metadata_schema(metadata)

input_cells <- as.character(raw_meta$Original_Cell_ID)
input_features <- as.character(feature_meta$Original_Feature_ID)
if (length(input_cells) != config$expected_cells) {
  stop(
    "input cell count differs from the verified source record: expected ",
    config$expected_cells,
    ", observed ",
    length(input_cells)
  )
}
if (length(input_features) != config$expected_features) {
  stop(
    "input feature count differs from the verified source record: expected ",
    config$expected_features,
    ", observed ",
    length(input_features)
  )
}
validate_dataset_metadata(
  metadata,
  matrix_cells = input_cells,
  expected_cells = length(input_cells)
)

if (metadata_only) {
  if (!file.exists(canonical_metadata_file)) {
    stop(
      "metadata-only standardization requires the existing canonical ",
      "metadata: ",
      canonical_metadata_file
    )
  }
  if (!file.exists(object_format_file)) {
    stop(
      "metadata-only standardization requires processed_object_format.tsv"
    )
  }
  previous_metadata <- read_complete_table(canonical_metadata_file)
  analysis_donors <- length(unique(stats::na.omit(
    metadata$Global_Donor_ID[metadata$Analysis_Include]
  )))
  if (!is.null(config$expected_analysis_cells) &&
    sum(metadata$Analysis_Include) != config$expected_analysis_cells) {
    stop(
      "metadata-only analysis cell count differs from the verified ",
      "inclusion rules: expected ",
      config$expected_analysis_cells,
      ", observed ",
      sum(metadata$Analysis_Include)
    )
  }
  if (!is.null(config$expected_analysis_donors) &&
    analysis_donors != config$expected_analysis_donors) {
    stop(
      "metadata-only analysis donor count differs from the verified ",
      "inclusion rules: expected ",
      config$expected_analysis_donors,
      ", observed ",
      analysis_donors
    )
  }
  if (!identical(
    as.character(previous_metadata$Cells),
    as.character(metadata$Cells)
  ) || any(as.character(previous_metadata$Dataset) != dataset)) {
    stop("metadata-only standardization changed source cell identity/order")
  }
  metadata <- restore_deterministic_metadata_provenance(
    metadata,
    previous_metadata
  )
  protected_change_audit <- character()
  for (column in c(
    "Analysis_Include", "Global_Donor_ID", "Specimen_ID", "Library_ID"
  )) {
    if (!column %in% names(previous_metadata) ||
      !column %in% names(metadata)) {
      stop(
        "metadata-only standardization lacks protected column: ",
        column
      )
    }
    if (!identical(
      as.character(previous_metadata[[column]]),
      as.character(metadata[[column]])
    )) {
      protected_rules <- config$metadata_only_protected_changes
      rule <- if (is.null(protected_rules)) {
        NULL
      } else {
        protected_rules[[column]]
      }
      if (is.null(rule)) {
        stop(
          "metadata-only standardization changed protected column: ",
          column
        )
      }
      previous_unique <- length(unique(stats::na.omit(
        previous_metadata[[column]]
      )))
      new_unique <- length(unique(stats::na.omit(metadata[[column]])))
      if (!identical(
        as.integer(previous_unique),
        as.integer(rule$expected_previous_unique)
      ) || !identical(
        as.integer(new_unique),
        as.integer(rule$expected_new_unique)
      )) {
        stop(
          "approved metadata-only migration count mismatch for ",
          column
        )
      }
      for (version in c("previous", "new")) {
        source_column <- rule[[paste0(version, "_source_column")]]
        expected_source_unique <- rule[[paste0(
          "expected_", version, "_source_unique"
        )]]
        source_metadata <- if (version == "previous") {
          previous_metadata
        } else {
          metadata
        }
        if (is.null(source_column) ||
          !source_column %in% names(source_metadata) ||
          !identical(
            as.integer(length(unique(stats::na.omit(
              source_metadata[[source_column]]
            )))),
            as.integer(expected_source_unique)
          )) {
          stop(
            "approved metadata-only migration source audit failed for ",
            column,
            " (",
            version,
            ")"
          )
        }
      }
      evidence <- normalize_missing_metadata(rule$evidence)
      if (length(evidence) != 1L || is.na(evidence)) {
        stop(
          "approved metadata-only migration lacks evidence for ",
          column
        )
      }
      protected_change_audit <- c(
        protected_change_audit,
        paste0(
          column,
          ":",
          previous_unique,
          "->",
          new_unique,
          " [",
          evidence,
          "]"
        )
      )
    }
  }

  format_audit <- read.delim(
    object_format_file,
    sep = "\t",
    quote = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (nrow(format_audit) != 1L ||
    as.numeric(format_audit$Full_Cells[[1L]]) != nrow(metadata) ||
    as.numeric(format_audit$Analysis_Cells[[1L]]) !=
      sum(metadata$Analysis_Include)) {
    stop("processed-object manifest and rebuilt metadata counts differ")
  }

  previous_sha256 <- processed_file_sha256(canonical_metadata_file)
  common_columns <- intersect(names(previous_metadata), names(metadata))
  changed_columns <- common_columns[!vapply(
    common_columns,
    function(column) {
      identical(
        as.character(previous_metadata[[column]]),
        as.character(metadata[[column]])
      )
    },
    logical(1)
  )]
  added_columns <- setdiff(names(metadata), names(previous_metadata))
  removed_columns <- setdiff(names(previous_metadata), names(metadata))
  if (length(removed_columns) > 0L) {
    stop(
      "metadata-only standardization would remove columns: ",
      paste(removed_columns, collapse = ", ")
    )
  }

  write_metadata_table_atomic(metadata, canonical_metadata_file)
  write_canonical_metadata_crosswalks(
    metadata,
    processed_dir,
    source_meta = raw_meta
  )
  metadata_audit <- dataset_metadata_audit(metadata)
  metadata_audit$source_accession <- source_record$source_accession[[1L]]
  metadata_audit$source_repository <- source_record$source_repository[[1L]]
  metadata_audit$source_publication_doi <-
    source_record$publication_doi[[1L]]
  metadata_audit$source_publication_title <-
    source_record$verified_title[[1L]]
  metadata_audit$source_repository_record_url <-
    source_record$repository_record_url[[1L]]
  metadata_audit$mode <- "metadata_only"
  metadata_audit$previous_metadata_sha256 <- previous_sha256
  metadata_audit$canonical_metadata_sha256 <-
    processed_file_sha256(canonical_metadata_file)
  metadata_audit$added_columns <- paste(added_columns, collapse = ";")
  metadata_audit$changed_columns <- paste(
    changed_columns,
    collapse = ";"
  )
  metadata_audit$removed_columns <- ""
  metadata_audit$approved_protected_column_changes <- paste(
    protected_change_audit,
    collapse = ";"
  )
  metadata_audit$matrix_recomputed <- FALSE
  write_metadata_table_atomic(
    metadata_audit,
    metadata_standardization_audit_file
  )
  if (file.exists(processing_audit_file)) {
    processing_audit <- read.delim(
      processing_audit_file,
      sep = "\t",
      quote = "",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    if (nrow(processing_audit) != 1L ||
      as.character(processing_audit$Dataset[[1L]]) != dataset) {
      stop("processing audit differs from metadata-only dataset")
    }
    shared_audit_columns <- intersect(
      names(processing_audit),
      names(metadata_audit)
    )
    for (column in shared_audit_columns) {
      processing_audit[[column]] <- metadata_audit[[column]]
    }
    processing_audit$metadata_standardized_utc <- format(
      Sys.time(),
      "%Y-%m-%dT%H:%M:%SZ",
      tz = "UTC"
    )
    processing_audit$metadata_standardization_matrix_recomputed <- FALSE
    processing_audit$approved_protected_column_changes <-
      metadata_audit$approved_protected_column_changes
    write_metadata_table_atomic(
      processing_audit,
      processing_audit_file
    )
  }

  set_format_value <- function(name, value) {
    format_audit[[name]] <<- value
  }
  set_format_value(
    "Canonical_Metadata",
    basename(canonical_metadata_file)
  )
  set_format_value("Feature_Metadata", basename(var_file))
  set_format_value("Age_Crosswalk", "age_crosswalk.tsv")
  set_format_value("Region_Crosswalk", "region_crosswalk.tsv")
  set_format_value(
    "Source_Cell_Type_Crosswalk",
    "source_cell_type_crosswalk.tsv"
  )
  set_format_value(
    "Technical_Batch_Crosswalk",
    "technical_batch_crosswalk.tsv"
  )
  set_format_value(
    "Technical_Batch_Design_Audit",
    "technical_batch_design_audit.tsv"
  )
  set_format_value(
    "Source_Technical_Field_Inventory",
    "source_technical_field_inventory.tsv"
  )
  set_format_value(
    "Biological_Unit_Crosswalk",
    "donor_specimen_library_crosswalk.tsv.gz"
  )
  set_format_value(
    "Metadata_Standardization_Audit",
    basename(metadata_standardization_audit_file)
  )
  analysis_object_name <- normalize_missing_metadata(
    format_audit$Analysis_Object[[1L]]
  )
  format_type <- if (
    identical(
      as.character(format_audit$Matrix_Class[[1L]]),
      "dgCMatrix_shards"
    )
  ) {
    "sharded_full_object_with_analysis_flag"
  } else if (is.na(analysis_object_name) &&
    as.numeric(format_audit$Analysis_Cells[[1L]]) == 0) {
    "single_full_object_without_analysis"
  } else {
    "paired_full_analysis_objects"
  }
  set_format_value("Format_Type", format_type)
  full_object_name <- normalize_missing_metadata(
    format_audit$Full_Object[[1L]]
  )
  if (is.na(full_object_name) ||
    grepl("^/|(^|/)[.][.](/|$)", full_object_name)) {
    stop("processed-object manifest has an unsafe full object")
  }
  full_object_path <- file.path(processed_dir, full_object_name)
  if (!file.exists(full_object_path)) {
    stop("processed full object is missing: ", full_object_path)
  }
  full_object_sha256 <- processed_file_sha256(full_object_path)
  if ("Full_Object_SHA256" %in% names(format_audit)) {
    recorded_sha256 <- normalize_missing_metadata(
      format_audit$Full_Object_SHA256[[1L]]
    )
    if (!is.na(recorded_sha256) &&
      !identical(recorded_sha256, full_object_sha256)) {
      stop("processed full object SHA-256 differs from its manifest")
    }
  }
  object_refresh_files <- character()
  if (identical(format_type, "sharded_full_object_with_analysis_flag") &&
    isTRUE(config$refresh_shard_metadata_on_metadata_only)) {
    shard_refresh <- refresh_processed_shard_metadata(
      processed_dir = processed_dir,
      metadata = metadata,
      expected_dataset = dataset,
      verify_sha256 = TRUE
    )
    object_refresh_files <- shard_refresh$Changed_Files
    full_object_sha256 <- processed_file_sha256(full_object_path)
    set_format_value(
      "Shard_Metadata_Refreshed_UTC",
      shard_refresh$Refreshed_UTC
    )
    set_format_value("Shard_Matrix_Recomputed", FALSE)
  }
  if ("Matrix_Recomputed" %in% names(format_audit)) {
    recorded_recomputed <- as.logical(
      format_audit$Matrix_Recomputed[[1L]]
    )
    if (!is.na(recorded_recomputed) && !recorded_recomputed) {
      stop(
        "processed-object manifest contradicts the complete matrix ",
        "materialization workflow"
      )
    }
  }
  set_format_value("Full_Object_SHA256", full_object_sha256)
  set_format_value("Matrix_Recomputed", TRUE)
  set_format_value(
    "Metadata_Schema_Version",
    brainomics_metadata_schema_version()
  )
  set_format_value("Source_Accession", source_record$source_accession[[1L]])
  set_format_value("Source_Repository", source_record$source_repository[[1L]])
  set_format_value(
    "Source_Publication_DOI",
    source_record$publication_doi[[1L]]
  )
  set_format_value(
    "Source_Publication_Title",
    source_record$verified_title[[1L]]
  )
  set_format_value(
    "Source_Publication_Journal",
    source_record$journal[[1L]]
  )
  set_format_value(
    "Source_Publication_Year",
    source_record$publication_year[[1L]]
  )
  set_format_value(
    "Source_Repository_Record_URL",
    source_record$repository_record_url[[1L]]
  )
  set_format_value("Metadata_Standardization_Matrix_Recomputed", FALSE)
  set_format_value(
    "Metadata_Standardized_UTC",
    format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  )
  write_metadata_table_atomic(format_audit, object_format_file)
  write_current_dataset_summary(
    objects = if (identical(format_type, "paired_full_analysis_objects")) {
      2L
    } else {
      1L
    },
    full_cells = format_audit$Full_Cells[[1L]],
    analysis_cells = format_audit$Analysis_Cells[[1L]],
    features = length(input_features),
    matrix_class = format_audit$Matrix_Class[[1L]],
    status = "metadata sidecars regenerated and validated"
  )
  writeLines(
    capture.output(sessionInfo()),
    metadata_standardization_session_file
  )
  refresh_processed_bundle_sha256(
    processed_dir,
    c(
      basename(canonical_metadata_file),
      basename(obs_file),
      basename(var_file),
      basename(export_audit_file),
      "age_crosswalk.tsv",
      "region_crosswalk.tsv",
      "source_cell_type_crosswalk.tsv",
      "technical_batch_crosswalk.tsv",
      "technical_batch_design_audit.tsv",
      "source_technical_field_inventory.tsv",
      "donor_specimen_library_crosswalk.tsv.gz",
      basename(metadata_standardization_audit_file),
      basename(object_format_file),
      basename(metadata_standardization_summary_file),
      basename(metadata_standardization_session_file),
      basename(processing_audit_file),
      object_refresh_files
    )
  )
  thisutils::log_message("", 
    paste0(
      dataset,
      " metadata standardized for all ",
      nrow(metadata),
      " cells; matrices unchanged"
    ),
    message_type = "success"
  )
  quit(save = "no", status = 0L)
}

counts_h5 <- open_h5ad_counts(
  h5ad_file = h5ad_file,
  matrix_group = config$matrix_group
)
if (!identical(colnames(counts_h5), input_cells)) {
  stop("H5AD reader changed or reordered input cells")
}
if (!identical(rownames(counts_h5), input_features)) {
  stop("H5AD reader changed or reordered input features")
}

matrix_storage <- if (is.null(config$matrix_storage)) {
  "dgCMatrix"
} else {
  config$matrix_storage
}
if (identical(matrix_storage, "dgCMatrix_shards")) {
  source("processing/process_layered_h5ad.R")
  process_layered_h5ad(
    dataset = dataset,
    config = config,
    counts_h5 = counts_h5,
    metadata = metadata,
    feature_meta = feature_meta,
    input_cells = input_cells,
    input_features = input_features,
    processed_dir = processed_dir,
    canonical_metadata_file = canonical_metadata_file,
    processing_audit_file = processing_audit_file,
    object_format_file = object_format_file,
    full_object_file = full_object_file,
    analysis_object_file = analysis_object_file,
    external_metadata_audit = external_metadata_audit,
    source_record = source_record,
    source_meta = raw_meta
  )
  quit(save = "no", status = 0L)
}
if (!identical(matrix_storage, "dgCMatrix")) {
  stop("unsupported processed matrix storage: ", matrix_storage)
}

source_matrix_checksum <- BPCells::checksum(counts_h5)
source_column_sums <- as.numeric(BPCells::colSums(counts_h5))
source_row_sums <- as.numeric(BPCells::rowSums(counts_h5))
counts <- materialize_sparse_counts(
  matrix = counts_h5,
  expected_features = length(input_features),
  expected_cells = length(input_cells),
  expected_feature_names = input_features,
  expected_cell_names = input_cells
)
matrix_content_verification <- verify_materialized_matrix_content(
  source_matrix = counts_h5,
  materialized_matrix = counts,
  block_cells = 5000L
)
materialized_iterable <- BPCells::write_matrix_memory(
  BPCells::convert_matrix_type(counts, "uint32_t")
)
materialized_matrix_checksum <- BPCells::checksum(
  materialized_iterable
)
materialized_column_sums <- as.numeric(Matrix::colSums(counts))
materialized_row_sums <- as.numeric(Matrix::rowSums(counts))
source_total_counts <- sum(source_column_sums)
materialized_total_counts <- sum(materialized_column_sums)
materialized_nonzero_values <- as.numeric(Matrix::nnzero(counts))
if (!is.null(config$expected_nonzero_values) &&
  materialized_nonzero_values != config$expected_nonzero_values) {
  stop(
    "input nonzero-value count differs from the verified source record: ",
    "expected ", config$expected_nonzero_values, ", observed ",
    materialized_nonzero_values
  )
}
if (!is.null(config$expected_total_counts) &&
  materialized_total_counts != config$expected_total_counts) {
  stop(
    "input total count differs from the verified source record: expected ",
    config$expected_total_counts, ", observed ",
    materialized_total_counts
  )
}
if (!isTRUE(matrix_content_verification$verified) || !identical(
  source_column_sums,
  materialized_column_sums
) || !identical(
  source_row_sums,
  materialized_row_sums
) || !isTRUE(all.equal(
  source_total_counts,
  materialized_total_counts,
  tolerance = 0
))) {
  stop("complete H5AD-to-dgCMatrix content verification failed")
}
rm(
  counts_h5,
  materialized_iterable,
  source_column_sums,
  source_row_sums,
  materialized_column_sums,
  materialized_row_sums
)
gc()
validate_dataset_metadata(
  metadata,
  matrix_cells = colnames(counts),
  expected_cells = ncol(counts)
)

full_object <- create_processed_object(
  counts = counts,
  metadata = metadata,
  dataset = dataset,
  input_cells = input_cells,
  input_features = input_features
)
save_processed_object(full_object, full_object_file)
full_object_cells <- ncol(full_object)
full_object_features <- nrow(full_object)
rm(full_object)
gc()

included <- metadata$Analysis_Include
analysis_metadata <- metadata[included, , drop = FALSE]
analysis_donors <- length(unique(stats::na.omit(
  analysis_metadata$Global_Donor_ID
)))
if (!is.null(config$expected_analysis_cells) &&
  nrow(analysis_metadata) != config$expected_analysis_cells) {
  stop(
    "analysis cell count differs from the verified inclusion rules: ",
    "expected ",
    config$expected_analysis_cells,
    ", observed ",
    nrow(analysis_metadata)
  )
}
if (!is.null(config$expected_analysis_donors) &&
  analysis_donors != config$expected_analysis_donors) {
  stop(
    "analysis donor count differs from the verified inclusion rules: ",
    "expected ",
    config$expected_analysis_donors,
    ", observed ",
    analysis_donors
  )
}
analysis_object <- NULL
analysis_object_cells <- 0L
analysis_object_features <- 0L
analysis_source_total_counts <- 0
analysis_written_total_counts <- 0
analysis_source_checksum <- NA_character_
analysis_written_checksum <- NA_character_
if (sum(included) > 0L) {
  analysis_source_counts <- counts[, included, drop = FALSE]
  if (!inherits(analysis_source_counts, "dgCMatrix")) {
    analysis_source_counts <- methods::as(
      analysis_source_counts,
      "dgCMatrix"
    )
  }
  analysis_source_total_counts <- sum(
    Matrix::colSums(analysis_source_counts)
  )
  analysis_iterable <- BPCells::write_matrix_memory(
    BPCells::convert_matrix_type(
      analysis_source_counts,
      "uint32_t"
    )
  )
  analysis_source_checksum <- BPCells::checksum(analysis_iterable)
  analysis_counts <- analysis_source_counts
  if (!identical(
    as.numeric(dim(analysis_counts)),
    c(as.numeric(length(input_features)), as.numeric(sum(included)))
  ) || !identical(
    colnames(analysis_counts),
    analysis_metadata$Cells
  )) {
    stop("deterministic analysis filtering produced a cell mismatch")
  }
  analysis_written_total_counts <- sum(Matrix::colSums(analysis_counts))
  analysis_written_checksum <- analysis_source_checksum
  if (!isTRUE(all.equal(
    as.numeric(analysis_source_total_counts),
    as.numeric(analysis_written_total_counts),
    tolerance = 0
  ))) {
    stop("analysis dgCMatrix count-sum verification failed")
  }
  analysis_object <- create_processed_object(
    counts = analysis_counts,
    metadata = analysis_metadata,
    dataset = dataset,
    input_cells = analysis_metadata$Cells,
    input_features = input_features
  )
  save_processed_object(analysis_object, analysis_object_file)
  analysis_object_cells <- ncol(analysis_object)
  analysis_object_features <- nrow(analysis_object)
  rm(analysis_object, analysis_counts, analysis_source_counts)
  rm(analysis_iterable)
  gc()
}

audit <- dataset_metadata_audit(metadata)
audit$source_accession <- source_record$source_accession[[1L]]
audit$source_repository <- source_record$source_repository[[1L]]
audit$source_publication_doi <- source_record$publication_doi[[1L]]
audit$source_publication_title <- source_record$verified_title[[1L]]
audit$source_publication_journal <- source_record$journal[[1L]]
audit$source_publication_year <- source_record$publication_year[[1L]]
audit$source_repository_record_url <-
  source_record$repository_record_url[[1L]]
audit$input_features <- length(input_features)
audit$full_object_cells <- full_object_cells
audit$full_object_features <- full_object_features
audit$analysis_object_cells <- analysis_object_cells
audit$analysis_object_features <- analysis_object_features
audit$analysis_donors <- analysis_donors
audit$expected_analysis_cells <- config$expected_analysis_cells
audit$expected_analysis_donors <- config$expected_analysis_donors
audit$matrix_group <- config$matrix_group
audit$expected_input_cells <- config$expected_cells
audit$expected_input_features <- config$expected_features
audit$input_file <- basename(h5ad_file)
audit$input_bytes <- as.numeric(file.info(h5ad_file)$size)
audit$input_sha256 <- if (is.na(input_sha256)) {
  processed_file_sha256(h5ad_file)
} else {
  input_sha256
}
audit$all_input_cells_processed <-
  full_object_cells == length(input_cells)
audit$all_input_features_processed <-
  full_object_features == length(input_features)
audit$source_matrix_checksum <- source_matrix_checksum
audit$materialized_matrix_checksum <- materialized_matrix_checksum
audit$matrix_content_verification <-
  "exact dgCMatrix slot comparison in complete cell-column blocks"
audit$matrix_content_verification_block_cells <-
  matrix_content_verification$block_cells
audit$matrix_content_verification_blocks <-
  matrix_content_verification$blocks
audit$matrix_backend_checksums_equal <- identical(
  source_matrix_checksum,
  materialized_matrix_checksum
)
audit$source_total_counts <- as.numeric(source_total_counts)
audit$materialized_total_counts <- as.numeric(
  materialized_total_counts
)
audit$materialized_nonzero_values <- materialized_nonzero_values
audit$analysis_source_checksum <- analysis_source_checksum
audit$analysis_written_checksum <- analysis_written_checksum
audit$analysis_source_total_counts <-
  as.numeric(analysis_source_total_counts)
audit$analysis_written_total_counts <-
  as.numeric(analysis_written_total_counts)
audit$matrix_storage <- "dgCMatrix"
audit$processed_object_self_contained <- TRUE
audit$preprocessing_reader <- "BPCells H5AD streaming reader"
if (!is.null(external_metadata_audit)) {
  for (column in names(external_metadata_audit)) {
    audit[[paste0("external_", column)]] <-
      external_metadata_audit[[column]]
  }
}
write.table(
  audit,
  processing_audit_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
if (requireNamespace("data.table", quietly = TRUE)) {
  data.table::fwrite(
    metadata,
    canonical_metadata_file,
    sep = "\t",
    quote = FALSE,
    na = "NA",
    compress = "gzip"
  )
} else {
  connection <- gzfile(canonical_metadata_file, "wt")
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
}
write_canonical_metadata_crosswalks(
  metadata,
  processed_dir,
  source_meta = raw_meta
)
metadata_standardization_audit <- dataset_metadata_audit(metadata)
metadata_standardization_audit$source_accession <-
  source_record$source_accession[[1L]]
metadata_standardization_audit$source_repository <-
  source_record$source_repository[[1L]]
metadata_standardization_audit$source_publication_doi <-
  source_record$publication_doi[[1L]]
metadata_standardization_audit$source_publication_title <-
  source_record$verified_title[[1L]]
metadata_standardization_audit$source_repository_record_url <-
  source_record$repository_record_url[[1L]]
metadata_standardization_audit$mode <- "full_processing"
metadata_standardization_audit$previous_metadata_sha256 <- NA_character_
metadata_standardization_audit$canonical_metadata_sha256 <-
  processed_file_sha256(canonical_metadata_file)
metadata_standardization_audit$added_columns <- NA_character_
metadata_standardization_audit$changed_columns <- NA_character_
metadata_standardization_audit$removed_columns <- NA_character_
metadata_standardization_audit$matrix_recomputed <- TRUE
write_metadata_table_atomic(
  metadata_standardization_audit,
  metadata_standardization_audit_file
)
writeLines(
  capture.output(sessionInfo()),
  file.path(processed_dir, "sessionInfo.txt")
)
writeLines(
  capture.output(sessionInfo()),
  metadata_standardization_session_file
)

format_audit <- data.frame(
  Dataset = dataset,
  Format_Type = if (analysis_object_cells > 0L) {
    "paired_full_analysis_objects"
  } else {
    "single_full_object_without_analysis"
  },
  Full_Object = basename(full_object_file),
  Analysis_Object = if (analysis_object_cells > 0L) {
    basename(analysis_object_file)
  } else {
    NA_character_
  },
  Full_Cells = full_object_cells,
  Analysis_Cells = analysis_object_cells,
  Canonical_Metadata = basename(canonical_metadata_file),
  Feature_Metadata = basename(var_file),
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
  Metadata_Standardization_Audit =
    basename(metadata_standardization_audit_file),
  Full_Object_SHA256 = processed_file_sha256(full_object_file),
  Matrix_Class = "dgCMatrix",
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
  Metadata_Standardization_Matrix_Recomputed = TRUE,
  Metadata_Standardized_UTC = format(
    Sys.time(),
    "%Y-%m-%dT%H:%M:%SZ",
    tz = "UTC"
  ),
  Completed_UTC = format(
    Sys.time(),
    "%Y-%m-%dT%H:%M:%SZ",
    tz = "UTC"
  ),
  stringsAsFactors = FALSE
)
temporary_format_file <- paste0(
  object_format_file,
  ".tmp.",
  Sys.getpid()
)
write.table(
  format_audit,
  temporary_format_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)
if (file.exists(object_format_file)) {
  unlink(object_format_file)
}
if (!file.rename(temporary_format_file, object_format_file)) {
  stop("failed to install processed-object format audit")
}

if (!all(audit$all_input_cells_processed) ||
  !all(audit$all_input_features_processed)) {
  stop("full-data processing guard failed")
}

write_current_dataset_summary(
  objects = 1L + as.integer(analysis_object_cells > 0L),
  full_cells = full_object_cells,
  analysis_cells = analysis_object_cells,
  features = full_object_features,
  matrix_class = "dgCMatrix",
  status = "complete H5AD reconstruction and metadata validation passed"
)

bundle_files <- c(
  basename(full_object_file),
  if (analysis_object_cells > 0L) basename(analysis_object_file),
  basename(canonical_metadata_file),
  basename(var_file),
  basename(obs_file),
  basename(export_audit_file),
  basename(processing_audit_file),
  "age_crosswalk.tsv",
  "region_crosswalk.tsv",
  "source_cell_type_crosswalk.tsv",
  "technical_batch_crosswalk.tsv",
  "technical_batch_design_audit.tsv",
  "source_technical_field_inventory.tsv",
  "donor_specimen_library_crosswalk.tsv.gz",
  basename(metadata_standardization_audit_file),
  basename(metadata_standardization_session_file),
  basename(metadata_standardization_summary_file),
  "sessionInfo.txt",
  basename(object_format_file),
  if (file.exists(file.path(processed_dir, "counts_conversion_audit.json"))) {
    "counts_conversion_audit.json"
  }
)
bundle_files <- unique(bundle_files[!is.na(bundle_files)])
refresh_processed_bundle_sha256(processed_dir, bundle_files)

thisutils::log_message("", 
  paste0(
    dataset,
    " full processing completed: ",
    full_object_cells,
    " cells, ",
    full_object_features,
    " features; analysis object: ",
    analysis_object_cells,
    " cells"
  ),
  message_type = "success"
)
