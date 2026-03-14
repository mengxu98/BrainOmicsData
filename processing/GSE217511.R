source("functions/prepare_env.R")

data_dir <- "../../data/BrainOmicsData/raw/GSE217511/GSE217511"
res_dir <- check_dir("../../data/BrainOmicsData/processed/GSE217511/")

log_message("Start loading data...")

metadata_dir <- "../../data/BrainOmicsData/raw/GSE217511"
metadata_files <- list.files(
  metadata_dir,
  pattern = ".*_Seuratmetadata\\.csv$",
  full.names = TRUE
)

all_metadata_list <- list()
for (meta_file in metadata_files) {
  meta_data <- read.csv(
    meta_file,
    header = TRUE, stringsAsFactors = FALSE, row.names = 1
  )
  all_metadata_list[[basename(meta_file)]] <- meta_data
}

gsm_dirs <- list.dirs(data_dir, full.names = FALSE, recursive = FALSE)
gsm_dirs <- gsm_dirs[grepl("^GSM", gsm_dirs)]

gsm_sample_mapping_prelim <- list()
for (gsm_id in gsm_dirs) {
  gsm_dir <- file.path(data_dir, gsm_id)
  files <- list.files(
    gsm_dir,
    pattern = paste0("^", gsm_id, "_.*\\.(tsv|mtx)\\.gz$")
  )

  for (file in files) {
    match_result <- regmatches(
      file, regexec(paste0("^", gsm_id, "_([^_]+)_"), file)
    )
    if (length(match_result[[1]]) > 1) {
      sample_id <- match_result[[1]][2]
      combined_key <- paste0(gsm_id, "_", sample_id)
      gsm_sample_mapping_prelim[[combined_key]] <- list(
        gsm_id = gsm_id, sample_id = sample_id
      )
    }
  }
}

sample_to_gsm_prelim <- data.frame(
  sample_id = sapply(gsm_sample_mapping_prelim, function(x) x$sample_id),
  gsm_id = sapply(gsm_sample_mapping_prelim, function(x) x$gsm_id),
  stringsAsFactors = FALSE
)

common_cols <- Reduce(
  intersect, lapply(all_metadata_list, colnames)
)

combined_metadata_list <- list()
for (meta_file in names(all_metadata_list)) {
  meta_data <- all_metadata_list[[meta_file]]
  meta_data <- meta_data[, common_cols, drop = FALSE]

  old_rownames <- rownames(meta_data)
  cell_barcodes <- gsub("^GSE217511_[^.]+\\.csv\\.", "", old_rownames)
  cell_barcodes <- gsub("^GSM[0-9]+_", "", cell_barcodes)

  if ("sample" %in% colnames(meta_data)) {
    sample_matches <- match(meta_data$sample, sample_to_gsm_prelim$sample_id)
    matched_indices <- which(!is.na(sample_matches))

    if (length(matched_indices) > 0) {
      matched_gsm_ids <- sample_to_gsm_prelim$gsm_id[sample_matches[matched_indices]]

      new_rownames <- old_rownames
      new_rownames[matched_indices] <- paste0(
        matched_gsm_ids, "_", cell_barcodes[matched_indices]
      )

      unmatched_indices <- which(is.na(sample_matches))
      if (length(unmatched_indices) > 0) {
        file_prefix <- gsub(
          "^GSE217511_([^_]+)_Seuratmetadata\\.csv$", "\\1", meta_file
        )
        new_rownames[unmatched_indices] <- paste0(
          file_prefix, "_", cell_barcodes[unmatched_indices]
        )
      }

      rownames(meta_data) <- new_rownames
    } else {
      file_prefix <- gsub(
        "^GSE217511_([^_]+)_Seuratmetadata\\.csv$", "\\1", meta_file
      )
      rownames(meta_data) <- paste0(file_prefix, "_", cell_barcodes)
    }
  } else {
    file_prefix <- gsub(
      "^GSE217511_([^_]+)_Seuratmetadata\\.csv$", "\\1", meta_file
    )
    rownames(meta_data) <- paste0(file_prefix, "_", cell_barcodes)
  }

  combined_metadata_list[[meta_file]] <- meta_data
}

combined_metadata <- do.call(rbind, combined_metadata_list)

if (any(duplicated(rownames(combined_metadata)))) {
  log_message("Found duplicate rownames after merging, making them unique...")
  rownames(combined_metadata) <- make.unique(
    rownames(combined_metadata),
    sep = "_"
  )
}

series_matrix_file <- file.path(metadata_dir, "GSE217511_series_matrix.txt.gz")

con <- gzfile(series_matrix_file, "r")
lines <- readLines(con)
close(con)

sample_meta_lines <- lines[grepl("^!Sample_", lines)]
sample_meta_list <- list()
for (line in sample_meta_lines) {
  parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
  if (length(parts) > 1) {
    key <- gsub("^!Sample_", "", parts[1])
    values <- parts[-1]
    values <- gsub('^"|"$', "", values)

    if (key == "characteristics_ch1") {
      first_value <- values[1]
      if (grepl("^tissue:", first_value)) {
        tissue_values <- gsub("^tissue: ", "", values)
        if (!"tissue" %in% names(sample_meta_list)) {
          sample_meta_list[["tissue"]] <- tissue_values
        }
      } else if (grepl("^developmental stage:", first_value)) {
        age_values <- gsub("^developmental stage: ", "", values)
        if (!"developmental_stage" %in% names(sample_meta_list)) {
          sample_meta_list[["developmental_stage"]] <- age_values
        }
      } else if (grepl("^Sex:", first_value)) {
        sex_values <- gsub("^Sex: ", "", values)
        if (!"sex" %in% names(sample_meta_list)) {
          sample_meta_list[["sex"]] <- sex_values
        }
      } else if (grepl("^treatment:", first_value)) {
        treatment_values <- gsub("^treatment: ", "", values)
        if (!"treatment" %in% names(sample_meta_list)) {
          sample_meta_list[["treatment"]] <- treatment_values
        }
      } else {
        if (!"characteristics_ch1" %in% names(sample_meta_list)) {
          sample_meta_list[["characteristics_ch1"]] <- values
        }
      }
    } else {
      sample_meta_list[[key]] <- values
    }
  }
}

gsm_ids <- NULL
gsm_ids <- sample_meta_list[["geo_accession"]]

log_message("Found {.val {length(gsm_ids)}} GSM IDs")

series_meta_df <- data.frame(
  row.names = gsm_ids,
  stringsAsFactors = FALSE
)

for (key in names(sample_meta_list)) {
  series_meta_df[[key]] <- sample_meta_list[[key]]
}
columns_to_keep <- c(
  "geo_accession", "developmental_stage", "sex",
  "source_name_ch1", "tissue"
)
available_columns <- intersect(columns_to_keep, colnames(series_meta_df))
series_meta_df <- series_meta_df[, available_columns, drop = FALSE]
colnames(series_meta_df) <- c(
  "Sample_ID", "Age", "Sex", "Brain_Region", "Region"
)

gsm_sample_mapping <- list()
for (gsm_id in gsm_dirs) {
  gsm_dir <- file.path(data_dir, gsm_id)
  files <- list.files(
    gsm_dir,
    pattern = paste0("^", gsm_id, "_.*\\.(tsv|mtx)\\.gz$")
  )

  for (file in files) {
    match_result <- regmatches(file, regexec(paste0("^", gsm_id, "_([^_]+)_"), file))
    if (length(match_result[[1]]) > 1) {
      sample_id <- match_result[[1]][2]
      combined_key <- paste0(gsm_id, "_", sample_id)
      gsm_sample_mapping[[combined_key]] <- list(gsm_id = gsm_id, sample_id = sample_id)
    }
  }
}

sample_to_gsm <- data.frame(
  sample_id = sapply(gsm_sample_mapping, function(x) x$sample_id),
  gsm_id = sapply(gsm_sample_mapping, function(x) x$gsm_id),
  stringsAsFactors = FALSE
)
sample_to_gsm <- sample_to_gsm[!duplicated(sample_to_gsm$sample_id), ]

sample_matches <- match(combined_metadata$sample, sample_to_gsm$sample_id)
matched_indices <- which(!is.na(sample_matches))

matched_gsm_ids <- sample_to_gsm$gsm_id[sample_matches[matched_indices]]
series_matches <- match(matched_gsm_ids, rownames(series_meta_df))
valid_matches <- !is.na(series_matches)

for (col in colnames(series_meta_df)) {
  combined_metadata[[col]] <- NA
  combined_metadata[[col]][matched_indices[valid_matches]] <-
    series_meta_df[[col]][series_matches[valid_matches]]
}

old_rownames <- rownames(combined_metadata)
needs_cleaning <- !grepl("^GSM[0-9]+_", old_rownames)
if (any(needs_cleaning)) {
  new_rownames <- old_rownames
  cell_barcodes_cleaning <- gsub(
    "^GSE217511_[^.]+\\.csv\\.", "", new_rownames[needs_cleaning]
  )
  cell_barcodes_cleaning <- gsub(
    "^GSM[0-9]+_", "", cell_barcodes_cleaning
  )

  if ("sample" %in% colnames(combined_metadata)) {
    unmatched_samples <- combined_metadata$sample[needs_cleaning]
    sample_matches <- match(unmatched_samples, sample_to_gsm$sample_id)
    matched_gsm <- sample_to_gsm$gsm_id[sample_matches[!is.na(sample_matches)]]
    if (length(matched_gsm) > 0) {
      matched_indices <- which(needs_cleaning)[!is.na(sample_matches)]
      new_rownames[matched_indices] <- paste0(
        matched_gsm, "_", cell_barcodes_cleaning[!is.na(sample_matches)]
      )
    }
    unmatched_indices <- which(needs_cleaning)[is.na(sample_matches)]
    if (length(unmatched_indices) > 0) {
      unmatched_cell_barcodes <- cell_barcodes_cleaning[is.na(sample_matches)]
      for (i in seq_along(unmatched_indices)) {
        idx <- unmatched_indices[i]
        parts <- strsplit(old_rownames[idx], "_")[[1]]
        if (length(parts) > 0 && grepl("^GSM", parts[1])) {
          prefix <- parts[1]
        } else {
          prefix <- gsub("^GSE217511_([^_]+)_.*", "\\1", old_rownames[idx])
          if (prefix == old_rownames[idx]) prefix <- "Unknown"
        }
        new_rownames[idx] <- paste0(prefix, "_", unmatched_cell_barcodes[i])
      }
    }
  }
  new_rownames <- make.unique(new_rownames, sep = "_")
  rownames(combined_metadata) <- new_rownames
}

combined_metadata$Brain_Region <- gsub(
  "^Brain, ", "", combined_metadata$Brain_Region
)
combined_metadata$Brain_Region <- gsub(
  "^([a-z])", "\\U\\1", combined_metadata$Brain_Region,
  perl = TRUE
)

combined_metadata$Region <- gsub(
  "^Brain, ", "", combined_metadata$Region
)
combined_metadata$Region <- gsub(
  "^([a-z])", "\\U\\1", combined_metadata$Region,
  perl = TRUE
)

before_count <- nrow(combined_metadata)
combined_metadata <- combined_metadata[combined_metadata$Region != "Svz+caudate", , drop = FALSE]
after_count <- nrow(combined_metadata)
log_message(
  "Filtered out svz+caudate data: {.val {before_count}} -> {.val {after_count}} cells"
)

combined_metadata$Age <- gsub(
  " weeks gestation", " PCW", combined_metadata$Age
)
combined_metadata$Age <- gsub(
  "^(\\d+(?:\\.\\d+)?) years$", "\\1", combined_metadata$Age,
  perl = TRUE
)

combined_metadata$cell_barcode_original <- rownames(combined_metadata)
combined_metadata$cell_barcode_clean <- sub(
  "_(\\d+)$", "", combined_metadata$cell_barcode_original
)

filtered_samples <- list()
if ("sample" %in% colnames(combined_metadata) && nrow(combined_metadata) > 0) {
  unique_samples <- unique(combined_metadata$sample)
  for (sample_id in unique_samples) {
    sample_rows <- combined_metadata[combined_metadata$sample == sample_id, , drop = FALSE]
    if (nrow(sample_rows) > 0) {
      first_rownames <- rownames(sample_rows)[1]
      gsm_match <- regmatches(
        first_rownames, regexec("^(GSM[0-9]+)_", first_rownames)
      )
      if (length(gsm_match[[1]]) > 1) {
        gsm_id <- gsm_match[[1]][2]
        filtered_samples[[paste0(gsm_id, "_", sample_id)]] <- list(
          gsm_id = gsm_id, sample_id = sample_id
        )
      }
    }
  }
}

read_10x_sample <- function(gsm_id, gsm_dir, sample_id) {
  pattern_barcodes <- paste0(
    "^", gsm_id, "_", sample_id, "_barcodes\\.tsv\\.gz$"
  )
  pattern_features <- paste0(
    "^", gsm_id, "_", sample_id, "_features\\.tsv\\.gz$"
  )
  pattern_matrix <- paste0(
    "^", gsm_id, "_", sample_id, "_matrix\\.mtx\\.gz$"
  )

  all_files <- list.files(gsm_dir, full.names = TRUE)

  barcodes_file <- grep(pattern_barcodes, basename(all_files), value = FALSE)
  features_file <- grep(pattern_features, basename(all_files), value = FALSE)
  matrix_file <- grep(pattern_matrix, basename(all_files), value = FALSE)

  barcodes_file <- all_files[barcodes_file[1]]
  features_file <- all_files[features_file[1]]
  matrix_file <- all_files[matrix_file[1]]

  temp_dir <- file.path(
    tempdir(),
    paste0("GSE217511_", gsm_id, "_", sample_id, "_", Sys.getpid())
  )
  dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)

  file.copy(
    barcodes_file, file.path(temp_dir, "barcodes.tsv.gz"),
    overwrite = TRUE
  )
  file.copy(
    features_file, file.path(temp_dir, "features.tsv.gz"),
    overwrite = TRUE
  )
  file.copy(
    matrix_file, file.path(temp_dir, "matrix.mtx.gz"),
    overwrite = TRUE
  )

  counts <- Seurat::Read10X(
    data.dir = temp_dir, gene.column = 2
  )

  unlink(temp_dir, recursive = TRUE)

  if (!is.null(counts)) {
    colnames(counts) <- paste0(gsm_id, "_", colnames(counts))
  }

  return(counts)
}


combo_df <- data.frame(
  key = names(filtered_samples),
  gsm_id = sapply(filtered_samples, function(x) x$gsm_id),
  sample_id = sapply(filtered_samples, function(x) x$sample_id),
  stringsAsFactors = FALSE
)

log_message("Processing {.val {nrow(combo_df)}} GSM-sample combinations...")

objects_list <- thisutils::parallelize_fun(
  seq_len(nrow(combo_df)),
  fun = function(i) {
    gsm_id <- combo_df$gsm_id[i]
    sample_id <- combo_df$sample_id[i]
    key <- combo_df$key[i]

    log_message("Processing {.val {gsm_id}} sample {.val {sample_id}}")
    gsm_dir <- file.path(data_dir, gsm_id)

    counts <- read_10x_sample(gsm_id, gsm_dir, sample_id)

    sample_meta <- combined_metadata[combined_metadata$sample == sample_id, , drop = FALSE]
    if (nrow(sample_meta) == 0) {
      return(NULL)
    }

    keep_cells <- intersect(colnames(counts), sample_meta$cell_barcode_clean)
    if (length(keep_cells) == 0) {
      return(NULL)
    }

    counts <- counts[, keep_cells, drop = FALSE]

    obj <- CreateSeuratObject(
      counts = counts
    )

    return(list(key = key, obj = obj))
  },
  cores = 17
)


objects_list <- objects_list[!sapply(objects_list, is.null)]
objects_list_named <- list()
for (item in objects_list) {
  objects_list_named[[item$key]] <- item$obj
}
objects_list <- objects_list_named


object <- merge(
  objects_list[[1]],
  y = objects_list[-1],
  project = "GSE217511"
)
object <- JoinLayers(object)

common_cells <- intersect(
  combined_metadata$cell_barcode_clean, colnames(object)
)

matched_indices <- match(common_cells, combined_metadata$cell_barcode_clean)
metadata <- combined_metadata[common_cells, ]
rownames(metadata) <- combined_metadata$cell_barcode_clean[matched_indices]

metadata$Cells <- rownames(metadata)
metadata$Dataset <- "GSE217511"
metadata$Technology <- "10X Genomics"
metadata$Sequence <- "snRNA-seq"
metadata$Sample <- metadata$sample
metadata$Sex <- ifelse(metadata$Sex == "female", "Female", "Male")
metadata$CellType_raw <- metadata$celltypes

column_order <- c(
  "Cells", "Dataset", "Technology", "Sequence", "Sample",
  "Sample_ID", "CellType_raw", "Brain_Region", "Region", "Age", "Sex"
)
metadata <- metadata[, column_order]
metadata <- na.omit(metadata)

counts <- GetAssayData(object, layer = "counts")
counts <- counts[, metadata$Cells]
object <- CreateSeuratObject(
  counts = counts,
  meta.data = metadata
)

log_message("Save data...")
saveRDS(
  object,
  file.path(res_dir, "GSE217511_processed.rds")
)
