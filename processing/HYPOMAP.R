source("functions/prepare_env.R")

data_dir <- brainomics_data_path("raw/HYPOMAP")
res_dir <- check_dir(brainomics_data_path("processed/HYPOMAP/"))

thisutils::log_message("Start loading data...")

rds_file <- file.path(data_dir, "human_HYPOMAP_snRNASeq.rds")
h5ad_file <- file.path(data_dir, "human_HYPOMAP_snRNASeq.h5ad")

if (file.exists(rds_file)) {
  object <- readRDS(rds_file)
} else if (file.exists(h5ad_file)) {
  object <- scop::h5ad_to_srt(h5ad_file, verbose = TRUE)
  saveRDS(object, rds_file)
} else {
  stop(
    "No HYPOMAP input found. Run download/HYPOMAP.sh or provide ",
    rds_file
  )
}

metadata <- object@meta.data
required_source_columns <- c(
  "Dataset", "Sample_ID", "Donor_ID", "Technology",
  "celltype_annotation", "C0", "C1", "C2", "C3", "C4",
  "C0_named", "C1_named", "C2_named", "C3_named", "C4_named",
  "region", "age_years", "sex"
)
missing_source_columns <- setdiff(required_source_columns, names(metadata))
if (length(missing_source_columns) > 0L) {
  stop(
    "HYPOMAP source metadata is missing columns: ",
    paste(missing_source_columns, collapse = ", ")
  )
}

source_dataset <- as.character(metadata$Dataset)
source_dataset_counts <- table(source_dataset)
expected_source_counts <- c(Siletti = 121405L, Tadross = 311964L)
if (!identical(
  as.integer(source_dataset_counts[names(expected_source_counts)]),
  as.integer(expected_source_counts)
)) {
  stop(
    "HYPOMAP source composition differs from the audited public object: ",
    paste(names(source_dataset_counts), source_dataset_counts, collapse = "; ")
  )
}

# The public HYPOMAP object is an atlas assembled from the newly generated
# Tadross data and an already published Siletti whole-brain dataset. Only the
# Tadross arm is reconstructed here so that Siletti nuclei are not imported a
# second time as newly generated HYPOMAP data.
metadata <- metadata[source_dataset == "Tadross", , drop = FALSE]
metadata$Source_Study_Dataset <- as.character(metadata$Dataset)
metadata$Cells <- rownames(metadata)
metadata$Dataset <- "HYPOMAP"
metadata$Technology <- "10x Genomics Chromium 3' v3.1"
metadata$Sequence <- "snRNA-seq"
metadata$Sample <- metadata$Donor_ID
metadata$Source_Sample_ID <- metadata$Sample_ID
metadata$Assay_Type <- "single-nucleus RNA-seq"
metadata$Sequencing_Platform <- "Illumina NovaSeq 6000"
metadata$Library_Chemistry <- "10x Genomics Chromium Single Cell 3' v3.1"
metadata$CellType_raw <- metadata$C4_named
metadata$Source_Fine_Region <- as.character(metadata$region)
fine_region_is_anatomical <-
  !is.na(metadata$Source_Fine_Region) &
    nzchar(metadata$Source_Fine_Region) &
    metadata$Source_Fine_Region != "Vascular"
metadata$Brain_Region <- ifelse(
  fine_region_is_anatomical,
  metadata$Source_Fine_Region,
  "Hypothalamus"
)
metadata$Region <- metadata$Brain_Region
metadata$Brain_Region_Assignment_Level <- ifelse(
  fine_region_is_anatomical,
  "source fine hypothalamic region",
  "study-wide tissue fallback"
)
metadata$Brain_Region_Assignment_Evidence <- ifelse(
  fine_region_is_anatomical,
  "retained from source HYPOMAP region field",
  paste(
    "Tadross et al. 2025, DOI 10.1038/s41586-024-08504-8;",
    "all retained Tadross nuclei originate from human hypothalamus,",
    "while the source fine-region field is missing or non-anatomical"
  )
)
metadata$Age <- metadata$age_years
metadata$Sex <- metadata$sex

column_order <- c(
  "Cells", "Dataset", "Technology", "Sequence", "Sample",
  "Sample_ID", "Source_Study_Dataset", "Source_Sample_ID",
  "Assay_Type", "Sequencing_Platform", "Library_Chemistry",
  "celltype_annotation", "C0", "C1", "C2", "C3", "C4",
  "C0_named", "C1_named", "C2_named", "C3_named", "C4_named",
  "CellType_raw",
  "Brain_Region", "Region", "Age", "Sex"
)
metadata <- retain_source_metadata(metadata, column_order)

if (
  nrow(metadata) != 311964L ||
    length(unique(metadata$Sample)) != 8L ||
    length(unique(metadata$Sample_ID)) != 58L ||
    any(vapply(
      split(metadata$Sample, metadata$Sample_ID),
      function(value) length(unique(value)) != 1L,
      logical(1)
    ))
) {
  stop(
    "HYPOMAP Tadross arm must contain 311,964 nuclei, 8 donors and ",
    "58 donor-specific source samples"
  )
}

counts <- GetAssayData(object, layer = "counts")
counts <- counts[, metadata$Cells]
object <- CreateSeuratObject(
  counts = counts,
  meta.data = metadata
)

thisutils::log_message("Save data...")
saveRDS(
  object,
  file.path(res_dir, "HYPOMAP_processed.rds")
)
