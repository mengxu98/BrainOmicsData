source("functions/prepare_env.R")

data_dir <- brainomics_data_path("raw/GSE207334")
res_dir <- check_dir(brainomics_data_path("processed/Ma_et_al_2022/"))


thisutils::log_message("Start loading data...")

object2 <- readRDS(
  file.path(data_dir, "Ma_Sestan_mat.rds")
)
object2 <- UpdateSeuratObject(object2)


sample_info <- data.frame(
  Sample = c("HSB106", "HSB189", "HSB340", "HSB628"),
  Age = c(64, 36, 19, 50),
  Sex = c("Male", "Male", "Male", "Female"),
  stringsAsFactors = FALSE
)
sample_info$Sex_Source_Raw <- sample_info$Sex
sample_info$Sex_Source_Standardized <- sample_info$Sex
sample_info$Sex_Source_Column <- "source publication donor table"
sample_info$Sex_Assignment_Method <-
  "curated from source publication donor table; not inferred"

metadata <- object2@meta.data
if (!"repname" %in% names(metadata) || anyNA(metadata$repname)) {
  stop("Ma et al. 2022 source repname field is incomplete")
}
metadata$tech_rep <- as.character(metadata$repname)
metadata$Technical_Replicate_Source_Column <- "repname"
metadata$Technical_Replicate_Derivation_Rule <-
  "identity: retain the exact source repname value"
metadata$Cells <- rownames(metadata)
metadata$Dataset <- "Ma_et_al_2022"
metadata$Technology <- "10X Genomics"
metadata$Sequence <- "snRNA-seq"
metadata$Sample <- metadata$samplename
metadata$Sample_ID <- metadata$samplename

metadata <- merge(
  metadata, sample_info,
  by = "Sample", all.x = TRUE, sort = FALSE
)
rownames(metadata) <- metadata$Cells
metadata <- metadata[colnames(object2), , drop = FALSE]
if (!identical(rownames(metadata), colnames(object2)) ||
  !"tech_rep" %in% names(metadata) || anyNA(metadata$tech_rep)) {
  stop("Ma et al. 2022 technical-replicate metadata is incomplete")
}
replicate_map <- unique(data.frame(
  donor = as.character(metadata$Sample),
  technical_replicate = as.character(metadata$tech_rep),
  stringsAsFactors = FALSE
))
if (nrow(replicate_map) != 16L ||
  any(table(replicate_map$donor) != 4L) ||
  any(sub("_[1-4]$", "", replicate_map$technical_replicate) !=
    replicate_map$donor)) {
  stop("Ma et al. 2022 must contain four technical replicates per donor")
}

metadata$CellType_raw <- metadata$subtype
metadata$Brain_Region <- "Dorsolateral prefrontal cortex"
metadata$Region <- "Dorsolateral prefrontal cortex"
column_order <- c(
  "Cells", "Dataset", "Technology", "Sequence", "Sample",
  "Sample_ID", "CellType_raw", "class", "subclass", "subtype",
  "Brain_Region", "Region", "Age", "Sex"
)
metadata <- retain_source_metadata(metadata, column_order)

counts <- GetAssayData(object2, layer = "counts")
counts <- counts[, metadata$Cells]
object2 <- CreateSeuratObject(
  counts = counts,
  meta.data = metadata
)

thisutils::log_message("Save data...")
saveRDS(
  object2,
  file.path(res_dir, "Ma_et_al_2022_processed.rds")
)
