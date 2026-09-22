process_egas <- function() {
  region_codes <- c("Cer", "FC", "GE", "Hipp", "Thal")
  region_labels <- c(
    Cer = "Cerebellum",
    FC = "Frontal cortex",
    GE = "Ganglionic eminence",
    Hipp = "Hippocampus",
    Thal = "Thalamus"
  )
  layers <- list()
  metadata_rows <- list()
  for (code in region_codes) {
    matrix_file <- file.path(
      raw_dir,
      paste0("cameron_2022_snRNAseq_", code, "_raw_count_gEX_matrix.txt.gz")
    )
    metadata_file <- file.path(
      raw_dir,
      paste0("cameron_2022_snRNAseq_", code, "_metadata.txt.gz")
    )
    required_source_files(c(matrix_file, metadata_file))
    table <- read_space_count_table_with_gene_column(matrix_file)
    genes <- as.character(table[[1L]])
    cells <- names(table)[-1L]
    if (anyDuplicated(genes)) genes <- make.unique(genes)
    values <- as.matrix(table[, -1L, drop = FALSE])
    storage.mode(values) <- "numeric"
    counts <- methods::as(Matrix::Matrix(values, sparse = TRUE), "dgCMatrix")
    rownames(counts) <- genes
    colnames(counts) <- cells
    rm(values, table)
    source_meta <- read_source_table(metadata_file)
    index <- match(cells, source_meta$cells)
    if (anyNA(index) || anyDuplicated(source_meta$cells) ||
      nrow(source_meta) != length(cells)) {
      stop("EGAS matrix and metadata differ for region ", code)
    }
    source_meta <- source_meta[index, , drop = FALSE]
    donor <- sub("_.*$", "", cells)
    sample <- as.character(source_meta$sample)
    meta <- base_source_metadata(
      cells,
      donor,
      sample,
      sample,
      sample,
      "14-15 PCW",
      "Female",
      region_labels[[code]],
      "Karyotypically normal",
      source_meta$cellIDs
    )
    meta$Source_CellType_Original <- as.character(source_meta$cellIDs)
    meta$Source_CellType_Annotation_Origin <- "original_study_author"
    meta$Original_Technical_Batch_ID <- sub("^.*_(B[0-9]+)$", "\\1", sample)
    meta$Library_Chemistry <- "10x Chromium 3' v3"
    meta$Sequencing_Platform <- "Illumina NovaSeq 6000"
    meta$Analysis_Include <- TRUE
    meta$Analysis_Role <- "reference"
    meta$Exclusion_Reason <- NA_character_
    meta <- age_provenance(
      meta,
      "14-15 PCW",
      "14-15 postconception weeks",
      "postconceptional weeks",
      "postconceptional",
      "Cameron et al. 2023, DOI 10.1016/j.biopsych.2022.06.033",
      "study-level 14-15 PCW range retained; donor-level mapping unavailable",
      "study-level range; exact donor age unresolved",
      FALSE
    )
    layers[[code]] <- counts
    metadata_rows[[code]] <- meta
  }
  layers <- align_layer_features(layers)
  metadata <- bind_rows_by_name(metadata_rows)
  create_complete_source_object(layers, metadata, dataset)
}

if (sys.nframe() == 0L) {
  PROCESSING_DATASET <- "EGAS00001006537"
  source("processing/process_original_source.R")
}
