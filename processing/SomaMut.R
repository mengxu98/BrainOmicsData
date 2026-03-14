source("functions/prepare_env.R")

data_dir <- "../../data/BrainOmicsData/raw/SomaMut"
res_dir <- check_dir("../../data/BrainOmicsData/processed/SomaMut/")

log_message("Start loading data...")

object <- readRDS(
  file.path(
    data_dir, "pfc.clean.rds"
  )
)

metadata <- object@meta.data
metadata$Cells <- rownames(metadata)
metadata$Dataset <- "SomaMut"
metadata$Technology <- "10X Genomics"
metadata$Sequence <- "snRNA-seq"
metadata$Sample <- metadata$case
metadata$Sample_ID <- metadata$case
metadata$CellType_raw <- metadata$new_clusters3
metadata$Brain_Region <- "Prefrontal cortex"
metadata$Region <- "Prefrontal cortex"
metadata$Age <- metadata$age
metadata$Sex <- metadata$sex
metadata$Sex <- ifelse(metadata$Sex == "female", "Female", "Male")

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
  file.path(res_dir, "SomaMut_processed.rds")
)
