# data: https://brainscope.gersteinlab.org/


source("functions/prepare_env.R")

data_dir <- "../../data/BrainOmicsData/raw/BICCN"
res_dir <- check_dir("../../data/BrainOmicsData/processed/BICCN/")

log_message("Start loading data...")

counts <- readRDS(
  file.path(data_dir, "BICCN_mat.RDS")
)
metadata <- readRDS(
  file.path(data_dir, "BICCN_meta_share.RDS")
)
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$sample_id

metadata$Cells <- rownames(metadata)
metadata$Dataset <- "BICCN"
metadata$Technology <- "10X Genomics"
metadata$Sequence <- "snRNA-seq"
metadata$Sample <- metadata$donor
metadata$Sample_ID <- metadata$layer
metadata$Brain_Region <- metadata$region
metadata$Region <- "Dorsolateral Prefrontal Cortex"
metadata$Age <- "Unknown"
metadata$Sex <- metadata$sex
metadata$CellType_raw <- metadata$within_area_subclass

column_order <- c(
  "Cells", "Dataset", "Technology", "Sequence", "Sample",
  "Sample_ID", "CellType_raw", "Brain_Region", "Region", "Age", "Sex"
)

metadata <- metadata[, column_order]
metadata <- na.omit(metadata)

counts <- counts[, metadata$Cells]

object <- Seurat::CreateSeuratObject(
  counts = counts,
  meta.data = metadata
)
saveRDS(
  object,
  file.path(res_dir, "BICCN_processed.rds")
)
