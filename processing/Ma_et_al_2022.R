source("functions/prepare_env.R")

data_dir <- "../../data/BrainOmicsData/raw/GSE207334"
res_dir <- check_dir("../../data/BrainOmicsData/processed/Ma_et_al_2022/")


log_message("Start loading data...")

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

metadata <- object2@meta.data
metadata$Cells <- rownames(metadata)
metadata$Dataset <- "Ma_et_al_2022"
metadata$Technology <- "10X Genomics"
metadata$Sequence <- "snRNA-seq"
metadata$Sample <- metadata$samplename
metadata$Sample_ID <- metadata$samplename

metadata <- merge(
  metadata, sample_info,
  by = "Sample", all.x = TRUE
)
rownames(metadata) <- metadata$Cells

metadata$CellType_raw <- metadata$subclass
metadata$Brain_Region <- "Dorsolateral prefrontal cortex"
metadata$Region <- "Dorsolateral prefrontal cortex"
column_order <- c(
  "Cells", "Dataset", "Technology", "Sequence", "Sample",
  "Sample_ID", "CellType_raw", "Brain_Region", "Region", "Age", "Sex"
)
metadata <- metadata[, column_order]

metadata <- na.omit(metadata)

counts <- GetAssayData(object2, layer = "counts")
counts <- counts[, metadata$Cells]
object2 <- CreateSeuratObject(
  counts = counts,
  meta.data = metadata
)

log_message("Save data...")
saveRDS(
  object2,
  file.path(res_dir, "Ma_et_al_2022_processed.rds")
)
