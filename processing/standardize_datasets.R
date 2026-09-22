#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

source("functions/data_paths.R")
source("functions/metadata_schema.R")
source("functions/dataset_metadata.R")
source("functions/processed_object.R")
source("functions/integration.R")

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = NULL) {
  index <- match(flag, args)
  if (is.na(index) || index == length(args)) {
    return(default)
  }
  args[[index + 1L]]
}
has_flag <- function(flag) flag %in% args

worker_value <- value_after(
  "--workers",
  Sys.getenv("BRAINOMICS_DATASET_WORKERS", unset = "1")
)
dataset_workers <- suppressWarnings(as.integer(worker_value))
if (is.na(dataset_workers) || dataset_workers < 1L) {
  stop("--workers must be a positive integer")
}

selected_dataset <- value_after("--dataset")
overwrite <- has_flag("--overwrite")
combine_only <- has_flag("--combine-only")
processed_root <- value_after(
  "--processed-root",
  brainomics_data_path("processed")
)
crosswalk_file <- value_after(
  "--donor-crosswalk",
  file.path(
    brainomics_repo_root(),
    "data",
    "donor_crosswalk.tsv"
  )
)
summary_override <- value_after("--summary")
summary_file <- if (!is.null(summary_override)) {
  summary_override
} else if (!is.null(selected_dataset)) {
  file.path(
    processed_root,
    selected_dataset,
    "metadata_standardization_summary.tsv"
  )
} else {
  file.path(processed_root, "dataset_standardization_summary.tsv")
}
metadata_schema_version <- brainomics_metadata_schema_version()

datasets <- formal_datasets()
if (!is.null(selected_dataset)) {
  if (!selected_dataset %in% datasets) {
    stop("unknown formal dataset: ", selected_dataset)
  }
  datasets <- selected_dataset
}
if (!file.exists(crosswalk_file)) {
  stop("donor crosswalk is missing: ", crosswalk_file)
}
source_access_file <- file.path(
  brainomics_repo_root(),
  "data",
  "source_access_summary.tsv"
)
source_records <- stats::setNames(
  lapply(
    datasets,
    formal_source_access_record,
    source_access_file = source_access_file
  ),
  datasets
)
redistribution_evidence_file <- file.path(
  brainomics_repo_root(), "data", "source_redistribution_evidence.tsv"
)
redistribution_evidence <- read.delim(
  redistribution_evidence_file, sep = "\t", quote = "",
  stringsAsFactors = FALSE, check.names = FALSE
)
if (anyDuplicated(redistribution_evidence$Dataset) ||
  !all(c("Dataset", "Data_License", "Derived_Matrix_Scope") %in%
    names(redistribution_evidence))) {
  stop("source redistribution evidence has invalid columns or duplicate datasets")
}

write_large_tsv <- function(x, file) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    return(write_tsv(x, file))
  }
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(file, ".tmp.", Sys.getpid())
  on.exit(unlink(temporary), add = TRUE)
  data.table::fwrite(
    x,
    temporary,
    sep = "\t",
    quote = FALSE,
    na = "NA",
    compress = if (grepl("[.]gz$", file, ignore.case = TRUE)) {
      "gzip"
    } else {
      "none"
    }
  )
  if (file.exists(file)) {
    unlink(file)
  }
  if (!file.rename(temporary, file)) {
    stop("failed to atomically write table: ", file)
  }
  invisible(file)
}

read_standardized_table <- function(file) {
  read_brainomics_table(file)
}

standardize_retained_dataset <- function(dataset) {
  source_record <- source_records[[dataset]]
  dataset_dir <- file.path(processed_root, dataset)
  object_file <- file.path(
    dataset_dir,
    paste0(dataset, "_processed.rds")
  )
  raw_metadata_file <- file.path(
    dataset_dir,
    "metadata_raw.tsv.gz"
  )
  canonical_metadata_file <- file.path(
    dataset_dir,
    "metadata_canonical.tsv.gz"
  )
  feature_metadata_file <- file.path(
    dataset_dir,
    "features_raw.tsv.gz"
  )
  age_crosswalk_file <- file.path(
    dataset_dir,
    "age_crosswalk.tsv"
  )
  region_crosswalk_file <- file.path(
    dataset_dir,
    "region_crosswalk.tsv"
  )
  source_cell_type_crosswalk_file <- file.path(
    dataset_dir,
    "source_cell_type_crosswalk.tsv"
  )
  technical_batch_crosswalk_file <- file.path(
    dataset_dir,
    "technical_batch_crosswalk.tsv"
  )
  technical_batch_design_audit_file <- file.path(
    dataset_dir,
    "technical_batch_design_audit.tsv"
  )
  source_technical_field_inventory_file <- file.path(
    dataset_dir,
    "source_technical_field_inventory.tsv"
  )
  biological_unit_crosswalk_file <- file.path(
    dataset_dir,
    "donor_specimen_library_crosswalk.tsv.gz"
  )
  processing_audit_file <- file.path(
    dataset_dir,
    "processing_audit.tsv"
  )
  format_file <- file.path(
    dataset_dir,
    "processed_object_format.tsv"
  )
  session_file <- file.path(
    dataset_dir,
    "metadata_standardization_sessionInfo.txt"
  )

  shard_manifest_file <- file.path(
    dataset_dir,
    "count_shard_manifest.tsv"
  )
  if (file.exists(shard_manifest_file) && file.exists(format_file)) {
    current_format <- read_processed_object_format(dataset_dir)
    if (identical(
      current_format$Format_Type[[1L]],
      "sharded_full_object_with_analysis_flag"
    )) {
      result <- validate_processed_object_bundle(
        dataset_dir,
        expected_dataset = dataset
      )
      result$Status <- "existing complete count-shard bundle validated"
      return(result)
    }
  }

  if (!file.exists(object_file)) {
    stop("retained processed object is missing: ", object_file)
  }
  if (file.exists(format_file) && !overwrite) {
    current_format <- read.delim(
      format_file,
      sep = "\t",
      quote = "",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    current_version <- if (
      "Metadata_Schema_Version" %in% names(current_format)
    ) {
      as.character(current_format$Metadata_Schema_Version[[1L]])
    } else {
      NA_character_
    }
    schema_current <- identical(
      current_version,
      metadata_schema_version
    ) &&
      "Source_Cell_Type_Crosswalk" %in% names(current_format) &&
      file.exists(source_cell_type_crosswalk_file) &&
      "Technical_Batch_Crosswalk" %in% names(current_format) &&
      file.exists(technical_batch_crosswalk_file) &&
      "Technical_Batch_Design_Audit" %in% names(current_format) &&
      file.exists(technical_batch_design_audit_file) &&
      "Source_Technical_Field_Inventory" %in% names(current_format) &&
      file.exists(source_technical_field_inventory_file)
    if (schema_current) {
      result <- validate_processed_object_bundle(
        dataset_dir,
        expected_dataset = dataset
      )
      result$Status <- "existing sidecars validated"
      return(result)
    }
    message(
      dataset,
      " metadata sidecars use an older schema; regenerating sidecars only"
    )
  }

  object_sha256_before <- processed_file_sha256(object_file)
  crosswalk_sha256 <- processed_file_sha256(crosswalk_file)
  object <- load_processed_object(object_file)
  original_meta <- object[[]]
  if (!identical(rownames(original_meta), colnames(object))) {
    stop("source object metadata and cell order differ for ", dataset)
  }
  dataset_values <- unique(as.character(original_meta$Dataset))
  if (length(dataset_values) != 1L ||
    !identical(dataset_values, dataset)) {
    stop("source object dataset label differs for ", dataset)
  }

  raw_export <- data.frame(
    Source_Object_Cell_ID = colnames(object),
    original_meta,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  canonical <- canonicalize_existing_reference_metadata(
    object,
    dataset,
    donor_crosswalk_file = crosswalk_file,
    source_record = source_record
  )
  if (!identical(canonical$Cells, colnames(object)) ||
    !identical(rownames(canonical), colnames(object))) {
    stop("canonical metadata changed the source cell order for ", dataset)
  }
  feature_metadata <- data.frame(
    Feature_Index = seq_len(nrow(object)),
    Original_Feature_ID = rownames(object),
    stringsAsFactors = FALSE
  )
  crosswalks <- canonical_metadata_crosswalks(canonical)
  age_crosswalk <- crosswalks$age
  region_crosswalk <- crosswalks$region
  source_cell_type_crosswalk <- crosswalks$source_cell_type
  technical_batch_crosswalk <- crosswalks$technical_batch
  biological_unit_crosswalk <- crosswalks$biological_unit
  metadata_audits <- canonical_metadata_audits(
    canonical,
    source_meta = raw_export
  )

  known_donor <- !is.na(normalize_missing_metadata(
    canonical$Global_Donor_ID
  ))
  known_specimen <- !is.na(normalize_missing_metadata(
    canonical$Specimen_ID
  ))
  known_library <- !is.na(normalize_missing_metadata(
    canonical$Library_ID
  ))
  known_source_record <- !is.na(normalize_missing_metadata(
    canonical$Source_Record_ID
  ))
  known_technical_batch <- !is.na(normalize_missing_metadata(
    canonical$Technical_Batch_ID
  ))
  age_reported <- !is.na(normalize_missing_metadata(
    canonical$Age
  ))
  age_parsed <- !is.na(canonical$Age_num)
  age_interval_assigned <- !is.na(normalize_missing_metadata(
    canonical$Stage
  ))
  source_cell_type_available <- as.logical(
    canonical$source_cell_type_available
  )
  source_cell_type_original_available <- !is.na(
    normalize_missing_metadata(canonical$source_cell_type_original_label)
  )
  source_cell_type_level_1_available <- !is.na(
    normalize_missing_metadata(canonical$source_cell_type_level_1)
  )
  source_cell_type_level_2_available <- !is.na(
    normalize_missing_metadata(canonical$source_cell_type_level_2)
  )
  source_cell_type_level_3_available <- !is.na(
    normalize_missing_metadata(canonical$source_cell_type_level_3)
  )
  source_cluster_available <- !is.na(
    normalize_missing_metadata(canonical$source_cluster_id)
  )
  region_ontology_mapped <- !is.na(normalize_missing_metadata(
    canonical$brain_region_ontology_id
  ))
  included <- as.logical(canonical$Analysis_Include)
  technical_batch_eligible <- grepl(
    "^eligible([ :]|$)",
    tolower(as.character(canonical$Technical_Batch_Use_Status))
  )
  audited_single_value <- function(column, fallback) {
    if (!column %in% names(original_meta)) {
      return(fallback)
    }
    values <- unique(stats::na.omit(normalize_missing_metadata(
      original_meta[[column]]
    )))
    if (length(values) != 1L) {
      stop(dataset, " must have exactly one retained ", column, " value")
    }
    values[[1L]]
  }
  retained_object_scope <- audited_single_value(
    "Source_Preprocessing_Retained_Object_Scope",
    paste(
      "previously processed retained object;",
      "not claimed to represent raw/pre-QC cells"
    )
  )
  raw_to_retained_attrition_status <- audited_single_value(
    "Source_Preprocessing_Attrition_Status",
    paste(
      "not reconstructible from this retained object alone;",
      "requires source-specific raw preprocessing audit"
    )
  )
  diagnosis_audit_status <- audited_single_value(
    "Source_Preprocessing_Diagnosis_Audit_Status",
    paste(
      "diagnosis not retained in the object;",
      "source-level inclusion ledger required"
    )
  )

  audit <- data.frame(
    Dataset = dataset,
    Source_Object = basename(object_file),
    Source_Object_Bytes = as.numeric(file.info(object_file)$size),
    Source_Object_SHA256 = object_sha256_before,
    Donor_Crosswalk_SHA256 = crosswalk_sha256,
    Source_Accession = source_record$source_accession[[1L]],
    Source_Repository = source_record$source_repository[[1L]],
    Source_Publication_DOI = source_record$publication_doi[[1L]],
    Source_Publication_Title = source_record$verified_title[[1L]],
    Source_Publication_Journal = source_record$journal[[1L]],
    Source_Publication_Year = source_record$publication_year[[1L]],
    Source_Repository_Record_URL =
      source_record$repository_record_url[[1L]],
    Retained_Input_Cells = ncol(object),
    Retained_Input_Features = nrow(object),
    Analysis_Cells = sum(canonical$Analysis_Include),
    Analysis_Excluded_Cells = sum(!canonical$Analysis_Include),
    Known_Donor_Cells = sum(known_donor),
    Unknown_Donor_Cells = sum(!known_donor),
    Known_Donors = length(unique(canonical$Global_Donor_ID[
      known_donor
    ])),
    Known_Specimens = length(unique(canonical$Specimen_ID[
      known_specimen
    ])),
    Known_Libraries = length(unique(canonical$Library_ID[
      known_library
    ])),
    Source_Records = length(unique(canonical$Source_Record_ID[
      known_source_record
    ])),
    Analysis_Donors = length(unique(canonical$Global_Donor_ID[
      included & known_donor
    ])),
    Analysis_Specimens = length(unique(canonical$Specimen_ID[
      included & known_specimen
    ])),
    Analysis_Libraries = length(unique(canonical$Library_ID[
      included & known_library
    ])),
    Analysis_Source_Records = length(unique(
      canonical$Source_Record_ID[
        included & known_source_record
      ]
    )),
    Verified_Technical_Batches = length(unique(
      canonical$Technical_Batch_ID[known_technical_batch]
    )),
    Analysis_Verified_Technical_Batches = length(unique(
      canonical$Technical_Batch_ID[
        included & known_technical_batch
      ]
    )),
    Technical_Batch_Correction_Eligible_Cells = sum(
      technical_batch_eligible
    ),
    Age_Reported_Cells = sum(age_reported),
    Age_Parsed_Cells = sum(age_parsed),
    Age_Interval_Assigned_Cells = sum(age_interval_assigned),
    Source_Cell_Type_Available_Cells = sum(source_cell_type_available),
    Source_Cell_Type_Unavailable_Cells = sum(!source_cell_type_available),
    Source_Cell_Type_Original_Label_Cells = sum(
      source_cell_type_original_available
    ),
    Source_Cell_Type_Level_1_Cells = sum(
      source_cell_type_level_1_available
    ),
    Source_Cell_Type_Level_2_Cells = sum(
      source_cell_type_level_2_available
    ),
    Source_Cell_Type_Level_3_Cells = sum(
      source_cell_type_level_3_available
    ),
    Source_Cluster_ID_Cells = sum(source_cluster_available),
    Brain_Region_Ontology_Mapped_Cells = sum(region_ontology_mapped),
    Brain_Region_Ontology_Pending_Cells = sum(!region_ontology_mapped),
    Original_Metadata_Columns = ncol(original_meta),
    Canonical_Metadata_Columns = ncol(canonical),
    Matrix_Class = processed_object_matrix_class(object),
    Matrix_Recomputed = FALSE,
    Retained_Object_Scope = retained_object_scope,
    Raw_To_Retained_Attrition_Status = raw_to_retained_attrition_status,
    Diagnosis_Audit_Status = diagnosis_audit_status,
    Metadata_Schema_Version = metadata_schema_version,
    Completed_UTC = format(
      Sys.time(),
      "%Y-%m-%dT%H:%M:%SZ",
      tz = "UTC"
    ),
    stringsAsFactors = FALSE
  )

  write_large_tsv(raw_export, raw_metadata_file)
  write_large_tsv(canonical, canonical_metadata_file)
  write_large_tsv(feature_metadata, feature_metadata_file)
  write_large_tsv(age_crosswalk, age_crosswalk_file)
  write_large_tsv(region_crosswalk, region_crosswalk_file)
  write_large_tsv(
    source_cell_type_crosswalk,
    source_cell_type_crosswalk_file
  )
  write_large_tsv(
    technical_batch_crosswalk,
    technical_batch_crosswalk_file
  )
  write_large_tsv(
    biological_unit_crosswalk,
    biological_unit_crosswalk_file
  )
  write_large_tsv(
    metadata_audits$technical_batch_design,
    technical_batch_design_audit_file
  )
  write_large_tsv(
    metadata_audits$source_technical_field_inventory,
    source_technical_field_inventory_file
  )
  write_tsv(audit, processing_audit_file)
  writeLines(capture.output(sessionInfo()), session_file)

  object_sha256_after <- processed_file_sha256(object_file)
  if (!identical(object_sha256_before, object_sha256_after)) {
    stop("retained expression object changed for ", dataset)
  }
  format_audit <- data.frame(
    Dataset = dataset,
    Format_Type =
      "single_retained_object_with_canonical_metadata",
    Full_Object = basename(object_file),
    Analysis_Object = NA_character_,
    Canonical_Metadata = basename(canonical_metadata_file),
    Feature_Metadata = basename(feature_metadata_file),
    Age_Crosswalk = basename(age_crosswalk_file),
    Region_Crosswalk = basename(region_crosswalk_file),
    Source_Cell_Type_Crosswalk =
      basename(source_cell_type_crosswalk_file),
    Technical_Batch_Crosswalk =
      basename(technical_batch_crosswalk_file),
    Technical_Batch_Design_Audit =
      basename(technical_batch_design_audit_file),
    Source_Technical_Field_Inventory =
      basename(source_technical_field_inventory_file),
    Biological_Unit_Crosswalk =
      basename(biological_unit_crosswalk_file),
    Raw_Metadata = basename(raw_metadata_file),
    Processing_Audit = basename(processing_audit_file),
    Full_Object_SHA256 = object_sha256_after,
    Donor_Crosswalk_SHA256 = crosswalk_sha256,
    Source_Accession = source_record$source_accession[[1L]],
    Source_Repository = source_record$source_repository[[1L]],
    Source_Publication_DOI = source_record$publication_doi[[1L]],
    Source_Publication_Title = source_record$verified_title[[1L]],
    Source_Publication_Journal = source_record$journal[[1L]],
    Source_Publication_Year = source_record$publication_year[[1L]],
    Source_Repository_Record_URL =
      source_record$repository_record_url[[1L]],
    Full_Cells = ncol(object),
    Analysis_Cells = sum(canonical$Analysis_Include),
    Matrix_Class = processed_object_matrix_class(object),
    Self_Contained = TRUE,
    Matrix_Recomputed = FALSE,
    Retained_Input_Scope = retained_object_scope,
    Metadata_Schema_Version = metadata_schema_version,
    Completed_UTC = format(
      Sys.time(),
      "%Y-%m-%dT%H:%M:%SZ",
      tz = "UTC"
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(format_audit, format_file)

  rm(
    object,
    original_meta,
    raw_export,
    canonical,
    feature_metadata,
    crosswalks,
    age_crosswalk,
    region_crosswalk,
    source_cell_type_crosswalk,
    technical_batch_crosswalk,
    biological_unit_crosswalk,
    metadata_audits
  )
  gc()
  result <- validate_processed_object_bundle(
    dataset_dir,
    expected_dataset = dataset
  )
  result$Status <- "sidecars generated and validated"
  result
}

standardize_configured_dataset <- function(dataset) {
  dataset_dir <- file.path(processed_root, dataset)
  format_file <- file.path(dataset_dir, "processed_object_format.tsv")
  required_sidecars <- c(
    "metadata_canonical.tsv.gz",
    "features_raw.tsv.gz",
    "age_crosswalk.tsv",
    "region_crosswalk.tsv",
    "source_cell_type_crosswalk.tsv",
    "technical_batch_crosswalk.tsv",
    "technical_batch_design_audit.tsv",
    "source_technical_field_inventory.tsv",
    "donor_specimen_library_crosswalk.tsv.gz",
    "processing_audit.tsv"
  )
  if (file.exists(format_file) && !overwrite) {
    current_format <- read.delim(
      format_file,
      sep = "\t",
      quote = "",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    current_version <- if (
      "Metadata_Schema_Version" %in% names(current_format)
    ) {
      as.character(current_format$Metadata_Schema_Version[[1L]])
    } else {
      NA_character_
    }
    sidecars_present <- all(file.exists(file.path(
      dataset_dir,
      required_sidecars
    )))
    if (identical(current_version, metadata_schema_version) &&
      sidecars_present) {
      result <- validate_processed_object_bundle(
        dataset_dir,
        expected_dataset = dataset
      )
      result$Status <- "existing configured bundle validated"
      return(result)
    }
    message(
      dataset,
      " configured metadata sidecars are incomplete or use an older schema; ",
      "regenerating sidecars"
    )
  }
  rscript <- file.path(R.home("bin"), "Rscript")
  script <- if (identical(dataset, "GSE294786")) {
    file.path(brainomics_repo_root(), "processing", "GSE294786.R")
  } else {
    file.path(brainomics_repo_root(), "processing", "process_h5ad.R")
  }
  script_args <- if (identical(dataset, "GSE294786")) {
    if (overwrite) c("true", "--metadata-only") else "--metadata-only"
  } else {
    c(
      dataset,
      if (overwrite) "true" else "false",
      "--metadata-only"
    )
  }
  status <- system2(
    rscript,
    c(shQuote(script), script_args)
  )
  if (status != 0L) {
    stop("configured metadata standardization failed for ", dataset)
  }
  result <- validate_processed_object_bundle(
    dataset_dir,
    expected_dataset = dataset
  )
  result$Status <- "configured sidecars regenerated and validated"
  result
}

refresh_standardized_bundle_manifest <- function(dataset_dir) {
  manifest_file <- file.path(dataset_dir, ".processed_bundle.sha256")
  if (!file.exists(manifest_file)) {
    return(invisible(FALSE))
  }

  format_audit <- read_processed_object_format(dataset_dir)
  file_columns <- intersect(
    c(
      "Full_Object", "Analysis_Object", "Canonical_Metadata",
      "Feature_Metadata", "Age_Crosswalk", "Region_Crosswalk",
      "Source_Cell_Type_Crosswalk", "Technical_Batch_Crosswalk",
      "Technical_Batch_Design_Audit",
      "Source_Technical_Field_Inventory",
      "Biological_Unit_Crosswalk", "Raw_Metadata",
      "Processing_Audit"
    ),
    names(format_audit)
  )
  changed_files <- unlist(
    format_audit[1L, file_columns, drop = FALSE],
    use.names = FALSE
  )
  changed_files <- normalize_missing_metadata(changed_files)
  changed_files <- changed_files[!is.na(changed_files)]
  changed_files <- c(
    changed_files,
    "processed_object_format.tsv",
    "metadata_standardization_sessionInfo.txt"
  )

  shard_manifest <- file.path(dataset_dir, "count_shard_manifest.tsv")
  if (file.exists(shard_manifest)) {
    shards <- read_processed_shard_manifest(dataset_dir)
    changed_files <- c(changed_files, as.character(shards$RDS_File))
  }
  changed_files <- unique(changed_files)
  missing_files <- changed_files[
    !file.exists(file.path(dataset_dir, changed_files))
  ]
  if (length(missing_files) > 0L) {
    stop(
      "standardized bundle references missing files: ",
      paste(missing_files, collapse = ", ")
    )
  }
  refresh_processed_bundle_sha256(dataset_dir, changed_files)
}

standardize_one_dataset <- function(dataset) {
  result <- if (dataset %in% existing_reference_datasets()) {
    standardize_retained_dataset(dataset)
  } else {
    standardize_configured_dataset(dataset)
  }
  refresh_standardized_bundle_manifest(file.path(processed_root, dataset))
  result
}

if (!combine_only) {
  standardize_index <- function(index) {
    dataset <- datasets[[index]]
    message(
      "[",
      index,
      "/",
      length(datasets),
      "] standardizing ",
      dataset
    )
    result <- standardize_one_dataset(dataset)
    gc()
    result
  }
  workers_used <- min(dataset_workers, length(datasets))
  summary_rows <- if (workers_used == 1L) {
    lapply(seq_along(datasets), standardize_index)
  } else {
    parallel::mclapply(
      seq_along(datasets),
      standardize_index,
      mc.cores = workers_used,
      mc.preschedule = FALSE
    )
  }
  failed <- vapply(summary_rows, inherits, logical(1), what = "try-error")
  if (any(failed)) {
    failures <- vapply(which(failed), function(index) {
      detail <- gsub(
        "[\r\n]+",
        " ",
        as.character(summary_rows[[index]])
      )
      paste0(datasets[[index]], ": ", detail)
    }, character(1))
    stop(
      "metadata standardization failed: ",
      paste(failures, collapse = "; ")
    )
  }
  summary <- bind_rows_by_name(summary_rows)
  write_tsv(summary, summary_file)
} else {
  if (!file.exists(summary_file)) {
    stop("standardization summary is missing: ", summary_file)
  }
  existing_summary <- read_standardized_table(summary_file)
  if (anyDuplicated(existing_summary$Dataset)) {
    stop("standardization summary contains duplicated datasets")
  }
  current_summary_row <- function(dataset) {
    dataset_summary_file <- file.path(
      processed_root,
      dataset,
      "metadata_standardization_summary.tsv"
    )
    row <- if (file.exists(dataset_summary_file)) {
      candidate <- read_standardized_table(dataset_summary_file)
      if (nrow(candidate) != 1L || candidate$Dataset[[1L]] != dataset) {
        stop("dataset standardization summary is invalid: ", dataset)
      }
      candidate
    } else {
      candidate <- existing_summary[
        existing_summary$Dataset == dataset, ,
        drop = FALSE
      ]
      if (nrow(candidate) != 1L) {
        stop("standardization summary lacks current dataset: ", dataset)
      }
      candidate
    }
    format_audit <- read_processed_object_format(file.path(
      processed_root,
      dataset
    ))
    validated_bundle <- validate_processed_object_bundle(
      processed_dir = file.path(processed_root, dataset),
      expected_dataset = dataset
    )
    required <- c(
      "Dataset", "Full_Cells", "Analysis_Cells", "Matrix_Class"
    )
    if (any(!required %in% names(row)) ||
      metadata_count_values(
        row$Full_Cells,
        paste0(dataset, " summary Full_Cells")
      ) != metadata_count_values(
        format_audit$Full_Cells,
        paste0(dataset, " format Full_Cells")
      ) ||
      metadata_count_values(
        row$Analysis_Cells,
        paste0(dataset, " summary Analysis_Cells")
      ) != metadata_count_values(
        format_audit$Analysis_Cells,
        paste0(dataset, " format Analysis_Cells")
      ) ||
      as.character(row$Matrix_Class[[1L]]) !=
        as.character(format_audit$Matrix_Class[[1L]]) ||
      metadata_count_values(
        validated_bundle$Full_Cells,
        paste0(dataset, " validated Full_Cells")
      ) != metadata_count_values(
        format_audit$Full_Cells,
        paste0(dataset, " format Full_Cells")
      ) ||
      metadata_count_values(
        validated_bundle$Analysis_Cells,
        paste0(dataset, " validated Analysis_Cells")
      ) != metadata_count_values(
        format_audit$Analysis_Cells,
        paste0(dataset, " format Analysis_Cells")
      )) {
      stop(
        "dataset summary differs from the current processed-object format: ",
        dataset
      )
    }
    if (identical(
      format_audit$Format_Type[[1L]],
      "sharded_full_object_with_analysis_flag"
    )) {
      manifest <- read_standardized_table(file.path(
        processed_root,
        dataset,
        format_audit$Full_Object[[1L]]
      ))
      manifest_features <- unique(metadata_count_values(
        manifest$Features,
        paste0(dataset, " shard Features")
      ))
      if (length(manifest_features) != 1L) {
        stop("count shards do not share one feature count: ", dataset)
      }
      row$Objects <- NA_real_
      row$Shards <- nrow(manifest)
      row$Features <- manifest_features[[1L]]
      row$Nonzero_Values <- sum(metadata_count_values(
        manifest$Nonzero_Values,
        paste0(dataset, " shard Nonzero_Values")
      ))
      row$Total_Counts <- sum(metadata_count_values(
        manifest$Total_Counts,
        paste0(dataset, " shard Total_Counts")
      ))
    } else {
      row$Objects <- if (identical(
        format_audit$Format_Type[[1L]],
        "paired_full_analysis_objects"
      )) {
        2
      } else {
        1
      }
      row$Shards <- NA_real_
    }
    count_columns <- intersect(
      c(
        "Objects", "Shards", "Full_Cells", "Analysis_Cells",
        "Features", "Nonzero_Values", "Total_Counts"
      ),
      names(row)
    )
    for (column in count_columns) {
      text <- normalize_missing_metadata(row[[column]])
      values <- rep(NA_real_, length(text))
      retained <- !is.na(text)
      if (any(retained)) {
        values[retained] <- metadata_count_values(
          text[retained],
          paste0(dataset, " summary ", column)
        )
      }
      row[[column]] <- values
    }
    row
  }
  summary <- bind_rows_by_name(lapply(datasets, current_summary_row))
  if (!identical(as.character(summary$Dataset), datasets)) {
    stop("standardization summary does not cover selected datasets")
  }
  write_tsv(summary, summary_file)
}

if (is.null(selected_dataset)) {
  combined_sidecars <- list(
    datasets_age_crosswalk.tsv = "age_crosswalk.tsv",
    datasets_region_crosswalk.tsv = "region_crosswalk.tsv",
    datasets_source_cell_type_crosswalk.tsv =
      "source_cell_type_crosswalk.tsv",
    datasets_technical_batch_crosswalk.tsv =
      "technical_batch_crosswalk.tsv",
    datasets_donor_specimen_library_crosswalk.tsv.gz =
      "donor_specimen_library_crosswalk.tsv.gz"
  )
  for (output_name in names(combined_sidecars)) {
    input_name <- combined_sidecars[[output_name]]
    combined <- bind_rows_by_name(
      lapply(datasets, function(dataset) {
        read_standardized_table(file.path(
          processed_root,
          dataset,
          input_name
        ))
      })
    )
    if (!"Cells" %in% names(combined)) {
      stop("combined sidecar lacks a Cells count column: ", output_name)
    }
    combined$Cells <- metadata_count_values(
      combined$Cells,
      paste0(output_name, " Cells")
    )
    expected_full_cells <- sum(metadata_count_values(
      summary$Full_Cells,
      "dataset standardization Full_Cells"
    ))
    if (sum(combined$Cells) != expected_full_cells) {
      stop("combined sidecar cell coverage failed: ", output_name)
    }
    if (identical(output_name, "datasets_region_crosswalk.tsv")) {
      canonical_region_columns <- intersect(
        c("brain_region_standardized", "brain_region_harmonized"),
        names(combined)
      )
      for (column in canonical_region_columns) {
        region <- normalize_missing_metadata(combined[[column]])
        known <- !is.na(region)
        case_variants <- split(region[known], tolower(region[known]))
        conflicting <- names(Filter(
          function(value) length(unique(value)) > 1L,
          case_variants
        ))
        if (length(conflicting) > 0L) {
          stop(
            "combined region crosswalk retains case-only canonical labels in ",
            column, ": ", paste(conflicting, collapse = ", ")
          )
        }
      }
    }
    write_large_tsv(
      combined,
      file.path(processed_root, output_name)
    )
  }
  combined_audits <- c(
    datasets_technical_batch_design_audit.tsv =
      "technical_batch_design_audit.tsv",
    datasets_source_technical_field_inventory.tsv =
      "source_technical_field_inventory.tsv"
  )
  for (output_name in names(combined_audits)) {
    input_name <- combined_audits[[output_name]]
    combined <- bind_rows_by_name(lapply(datasets, function(dataset) {
      read_standardized_table(file.path(
        processed_root,
        dataset,
        input_name
      ))
    }))
    summaries <- combined[
      combined$Record_Type == "dataset_summary", ,
      drop = FALSE
    ]
    if (!identical(as.character(summaries$Dataset), datasets)) {
      stop("combined metadata audit dataset coverage failed: ", output_name)
    }
    write_large_tsv(combined, file.path(processed_root, output_name))
  }

  # Source-by-source biological inclusion ledger. Diagnosis/status values are
  # preserved exactly as retained in canonical metadata; missing source values
  # remain explicit rather than being treated as controls. Access to raw data,
  # access to the processed workflow input, and redistribution of the derived
  # release are deliberately recorded as separate concepts.
  access_scope <- function(access_status, layer) {
    status <- tolower(normalize_missing_metadata(access_status))
    if (is.na(status)) {
      return("not reported")
    }
    if (identical(layer, "raw")) {
      if (grepl("restrict|control", status)) {
        return("controlled or restricted")
      }
      if (grepl("mixed", status)) {
        return("mixed public and controlled")
      }
      if (grepl("public|open|accessible", status)) {
        return("publicly accessible")
      }
      return(access_status)
    }
    if (grepl("restrict|control", status)) {
      "processed workflow input available under controlled or source-specific access"
    } else {
      "processed workflow input available through the recorded public source route"
    }
  }
  known_unique_count <- function(values) {
    values <- normalize_missing_metadata(values)
    length(unique(values[!is.na(values)]))
  }
  inclusion_ledger <- bind_rows_by_name(lapply(datasets, function(dataset) {
    canonical <- read_standardized_table(file.path(
      processed_root, dataset, "metadata_canonical.tsv.gz"
    ))
    required <- c(
      "Diagnosis_raw", "Analysis_Include", "Analysis_Role",
      "Exclusion_Reason", "Global_Donor_ID", "Specimen_ID", "Library_ID"
    )
    if (!all(required %in% names(canonical))) {
      stop("canonical metadata lacks inclusion-ledger fields: ", dataset)
    }
    diagnosis <- normalize_missing_metadata(canonical$Diagnosis_raw)
    diagnosis[is.na(diagnosis)] <- "Not reported by source"
    role <- normalize_missing_metadata(canonical$Analysis_Role)
    role[is.na(role)] <- "Not reported"
    reason <- normalize_missing_metadata(canonical$Exclusion_Reason)
    reason[is.na(reason)] <- "None recorded"
    key <- data.frame(
      Diagnosis_Source_Value = diagnosis,
      Analysis_Include = as.logical(canonical$Analysis_Include),
      Analysis_Role = role,
      Exclusion_Reason = reason,
      stringsAsFactors = FALSE
    )
    groups <- split(seq_len(nrow(canonical)), interaction(
      key$Diagnosis_Source_Value, key$Analysis_Include,
      key$Analysis_Role, key$Exclusion_Reason,
      drop = TRUE, lex.order = TRUE
    ))
    source_record <- source_records[[dataset]]
    evidence_row <- redistribution_evidence[
      redistribution_evidence$Dataset == dataset, , drop = FALSE
    ]
    access_override <- function(column, layer) {
      value <- normalize_missing_metadata(source_record[[column]])
      if (length(value) == 1L && !is.na(value)) {
        value
      } else {
        access_scope(source_record$access_status[[1L]], layer)
      }
    }
    bind_rows_by_name(lapply(groups, function(index) {
      data.frame(
        Dataset = dataset,
        Source_Accession = source_record$source_accession[[1L]],
        Source_Publication_DOI = source_record$publication_doi[[1L]],
        Diagnosis_Source_Value = key$Diagnosis_Source_Value[index[[1L]]],
        Analysis_Include = key$Analysis_Include[index[[1L]]],
        Analysis_Role = key$Analysis_Role[index[[1L]]],
        Exclusion_Reason = key$Exclusion_Reason[index[[1L]]],
        Cells_or_Nuclei = length(index),
        Known_Donors = known_unique_count(canonical$Global_Donor_ID[index]),
        Known_Specimens = known_unique_count(canonical$Specimen_ID[index]),
        Known_Libraries = known_unique_count(canonical$Library_ID[index]),
        Raw_Data_Access_Status = access_override(
          "raw_data_access_status", "raw"
        ),
        Processed_Input_Access_Status = access_override(
          "processed_input_access_status", "processed"
        ),
        Derived_Matrix_License = if (nrow(evidence_row) == 1L) {
          evidence_row$Data_License[[1L]]
        } else {
          NA_character_
        },
        Derived_Matrix_Scope = if (nrow(evidence_row) == 1L) {
          evidence_row$Derived_Matrix_Scope[[1L]]
        } else {
          NA_character_
        },
        stringsAsFactors = FALSE
      )
    }))
  }))
  expected_full_cells <- sum(metadata_count_values(
    summary$Full_Cells,
    "dataset standardization Full_Cells"
  ))
  if (sum(metadata_count_values(
    inclusion_ledger$Cells_or_Nuclei,
    "datasets inclusion ledger Cells_or_Nuclei"
  )) != expected_full_cells ||
    !identical(unique(as.character(inclusion_ledger$Dataset)), datasets)) {
    stop("combined inclusion ledger failed dataset or cell coverage")
  }
  write_large_tsv(
    inclusion_ledger,
    file.path(processed_root, "datasets_inclusion_ledger.tsv")
  )
}
message(
  if (combine_only) {
    "Combined crosswalk validation completed for "
  } else {
    "Metadata standardization completed for "
  },
  nrow(summary),
  " datasets; retained expression matrices were not modified"
)
