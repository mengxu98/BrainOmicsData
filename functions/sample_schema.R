normalize_sample_id_value <- function(dataset, value) {
  value <- trimws(as.character(value))
  value[is.na(value) | value == "" | value %in% c("NA", "NULL", "None")] <- "Unknown/not reported"
  paste(dataset, value, sep = ":")
}

add_sample_schema <- function(meta) {
  required <- c("Dataset", "Sample", "Sample_ID")
  missing <- setdiff(required, names(meta))
  if (length(missing) > 0) {
    stop("metadata is missing required sample columns: ", paste(missing, collapse = ", "))
  }

  dataset <- as.character(meta$Dataset)
  original_sample <- if ("Original_Sample" %in% names(meta)) {
    as.character(meta$Original_Sample)
  } else {
    as.character(meta$Sample)
  }
  original_sample_id <- if ("Original_Sample_ID" %in% names(meta)) {
    as.character(meta$Original_Sample_ID)
  } else {
    as.character(meta$Sample_ID)
  }

  donor_id <- normalize_sample_id_value(dataset, original_sample)
  biological_sample_id <- donor_id
  library_id <- normalize_sample_id_value(dataset, original_sample_id)
  sample_schema_rule <- rep("sample_field_as_donor_and_biological_sample", nrow(meta))

  use_original_sample_id_as_sample <- meta$Dataset %in% c(
    "HYPOMAP", "GSE199762", "GSE103723", "Li_et_al_2018",
    "GSE207334", "GSE217511", "GSE261983", "GSE97942"
  )
  biological_sample_id[use_original_sample_id_as_sample] <- normalize_sample_id_value(
    dataset[use_original_sample_id_as_sample],
    original_sample_id[use_original_sample_id_as_sample]
  )
  sample_schema_rule[use_original_sample_id_as_sample] <- "original_sample_id_as_biological_sample"

  donor_sample_only <- meta$Dataset %in% c(
    "GSE104276", "Nowakowski_et_al_2017", "GSE186538"
  )
  biological_sample_id[donor_sample_only] <- normalize_sample_id_value(
    dataset[donor_sample_only],
    original_sample[donor_sample_only]
  )
  library_id[donor_sample_only] <- normalize_sample_id_value(
    dataset[donor_sample_only],
    original_sample_id[donor_sample_only]
  )
  sample_schema_rule[donor_sample_only] <- "original_sample_id_is_cell_or_project_level"

  cell_level_without_sample <- meta$Dataset %in% c("GSE81475", "GSE67835")
  donor_id[cell_level_without_sample] <- normalize_sample_id_value(
    dataset[cell_level_without_sample],
    "Unknown/not reported"
  )
  biological_sample_id[cell_level_without_sample] <- normalize_sample_id_value(
    dataset[cell_level_without_sample],
    "Unknown/not reported"
  )
  library_id[cell_level_without_sample] <- normalize_sample_id_value(
    dataset[cell_level_without_sample],
    original_sample_id[cell_level_without_sample]
  )
  sample_schema_rule[cell_level_without_sample] <- "cell_level_identifier_no_reliable_sample_or_donor"

  is_gse296073 <- meta$Dataset == "GSE296073"
  if (any(is_gse296073)) {
    region <- if ("BrainRegion" %in% names(meta)) as.character(meta$BrainRegion) else "UnknownBrainRegion"
    biological_sample_id[is_gse296073] <- normalize_sample_id_value(
      dataset[is_gse296073],
      paste(original_sample[is_gse296073], region[is_gse296073], sep = "_")
    )
    library_id[is_gse296073] <- normalize_sample_id_value(
      dataset[is_gse296073],
      original_sample_id[is_gse296073]
    )
    sample_schema_rule[is_gse296073] <- "subject_brain_region_as_biological_sample_library_id_preserved"
  }

  meta$Original_Sample <- original_sample
  meta$Original_Sample_ID <- original_sample_id
  meta$Sample <- biological_sample_id
  meta$Donor_ID <- donor_id
  meta$Sample_ID <- biological_sample_id
  meta$Library_ID <- library_id
  meta$sample_schema_rule <- sample_schema_rule
  meta
}

sample_schema_audit <- function(meta) {
  required <- c("Dataset", "Original_Sample", "Original_Sample_ID", "Donor_ID", "Sample_ID", "Library_ID", "sample_schema_rule")
  missing <- setdiff(required, names(meta))
  if (length(missing) > 0) {
    stop("metadata is missing harmonized sample schema columns: ", paste(missing, collapse = ", "))
  }
  split_meta <- split(meta, meta$Dataset)
  do.call(rbind, lapply(split_meta, function(x) {
    data.frame(
      Dataset = unique(x$Dataset),
      cells = nrow(x),
      original_sample_count = length(unique(x$Original_Sample)),
      original_sample_id_count = length(unique(x$Original_Sample_ID)),
      donor_count = length(unique(x$Donor_ID)),
      biological_sample_count = length(unique(x$Sample_ID)),
      library_count = length(unique(x$Library_ID)),
      sample_schema_rule = paste(sort(unique(x$sample_schema_rule)), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }))
}
