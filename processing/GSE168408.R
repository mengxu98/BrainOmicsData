source("functions/data_paths.R")

repo_dir <- brainomics_repo_root()
raw_dir <- brainomics_data_path("raw", "GSE168408")
library_metadata_file <- file.path(
  raw_dir,
  "source_library_metadata.tsv"
)

read_geo_series_field <- function(lines, field) {
  matches <- lines[startsWith(lines, paste0(field, "\t"))]
  if (length(matches) != 1L) {
    stop("GSE168408 GEO series matrix has an invalid field: ", field)
  }
  values <- strsplit(matches, "\t", fixed = TRUE)[[1L]][-1L]
  sub('^"(.*)"$', "\\1", values)
}

read_geo_platform_map <- function(file) {
  if (!file.exists(file)) {
    stop("Missing GSE168408 GEO series matrix: ", file)
  }
  connection <- gzfile(file, "rt")
  lines <- tryCatch(
    readLines(connection, warn = FALSE),
    finally = close(connection)
  )
  gsm <- read_geo_series_field(lines, "!Sample_geo_accession")
  platform <- read_geo_series_field(lines, "!Sample_platform_id")
  title <- read_geo_series_field(lines, "!Sample_title")
  if (length(gsm) != length(platform) || length(gsm) != length(title)) {
    stop("GSE168408 GEO series fields have different lengths")
  }
  data.frame(
    Source_GEO_RNA_Record_ID = gsm,
    Source_GEO_Sample_Title = title,
    Source_GEO_Platform_ID = platform,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

build_source_library_metadata <- function() {
  series_files <- file.path(
    raw_dir,
    c(
      "GSE168408-GPL21697_series_matrix.txt.gz",
      "GSE168408-GPL24676_series_matrix.txt.gz"
    )
  )
  platform_map <- do.call(
    rbind,
    lapply(series_files, read_geo_platform_map)
  )
  if (anyDuplicated(platform_map$Source_GEO_RNA_Record_ID)) {
    stop("GSE168408 GEO records occur on more than one platform")
  }

  filtered_files <- list.files(
    raw_dir,
    pattern = paste0(
      "^GSM[0-9]+_RL[0-9]+_.*filtered.*[.]h5([.]gz)?$"
    ),
    recursive = TRUE,
    full.names = TRUE
  )
  source_name <- basename(filtered_files)
  matched <- grepl(
    paste0(
      "^GSM[0-9]+_RL[0-9]+_.*filtered.*[.]h5([.]gz)?$"
    ),
    source_name
  )
  if (length(source_name) < 32L || any(!matched)) {
    stop("GSE168408 filtered RNA source-file inventory is incomplete")
  }
  source_file_map <- data.frame(
    Source_GEO_RNA_Record_ID = sub(
      "^(GSM[0-9]+)_.*$",
      "\\1",
      source_name
    ),
    Source_RL_ID = sub(
      "^GSM[0-9]+_(RL[0-9]+)_.*$",
      "\\1",
      source_name
    ),
    Source_Filtered_Matrix_File = source_name,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  source_file_map <- unique(source_file_map)
  if (anyDuplicated(source_file_map$Source_GEO_RNA_Record_ID) ||
    anyDuplicated(source_file_map$Source_RL_ID)) {
    stop("GSE168408 RNA source file does not map one-to-one to a donor")
  }
  source_file_map <- merge(
    source_file_map,
    platform_map,
    by = "Source_GEO_RNA_Record_ID",
    all.x = TRUE,
    sort = FALSE
  )
  if (anyNA(source_file_map$Source_GEO_Platform_ID)) {
    stop("GSE168408 RNA source file lacks a GEO platform assignment")
  }
  platform_label <- c(
    GPL21697 = "Illumina NextSeq 550",
    GPL24676 = "Illumina NovaSeq 6000"
  )
  source_file_map$Source_Sequencing_Platform <- unname(
    platform_label[source_file_map$Source_GEO_Platform_ID]
  )
  if (anyNA(source_file_map$Source_Sequencing_Platform)) {
    stop("GSE168408 contains an unreviewed GEO platform")
  }
  source_file_map <- source_file_map[order(
    source_file_map$Source_GEO_RNA_Record_ID
  ), , drop = FALSE]
  write.table(
    source_file_map,
    library_metadata_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

build_source_library_metadata()

PROCESSING_DATASET <- "GSE168408"
source("processing/process_h5ad.R")
