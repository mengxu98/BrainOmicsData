source("functions/prepare_env.R")

data_dir <- "../../data/BrainOmicsData/raw/HYPOMAP"
res_dir <- check_dir("../../data/BrainOmicsData/processed/HYPOMAP/")

log_message("Start loading data...")

object <- readRDS(
  file.path(
    data_dir, "human_HYPOMAP_snRNASeq.rds"
  )
)

metadata <- object@meta.data
metadata <- na.omit(metadata)
metadata$Cells <- rownames(metadata)
metadata$Dataset <- "HYPOMAP"
metadata$Technology <- "10X Genomics"
metadata$Sequence <- "snRNA-seq"
metadata$Sample <- metadata$Donor_ID
metadata$CellType_raw <- metadata$celltype_annotation
metadata$Brain_Region <- metadata$region
metadata$Region <- metadata$region
metadata$Age <- metadata$age_years
metadata$Sex <- metadata$sex

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
  file.path(res_dir, "HYPOMAP_processed.rds")
)
