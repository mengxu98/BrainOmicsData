source("functions/prepare_env.R")

data_dir <- "../../data/BrainOmicsData/raw/GSE207334"
res_dir <- check_dir("../../data/BrainOmicsData/processed/GSE207334/")

log_message("Start loading data...")

rna_counts_file <- file.path(data_dir, "GSE207334_Multiome_rna_counts.mtx.gz")
rna_genes_file <- file.path(data_dir, "GSE207334_Multiome_rna_genes.txt.gz")
cell_meta_file <- file.path(data_dir, "GSE207334_Multiome_cell_meta.txt.gz")

rna_counts <- Matrix::readMM(gzfile(rna_counts_file))

rna_genes <- read.table(
  gzfile(rna_genes_file),
  header = FALSE, stringsAsFactors = FALSE
)[, 1]

rownames(rna_counts) <- rna_genes

metadata <- read.table(
  gzfile(cell_meta_file),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1
)

colnames(rna_counts) <- rownames(metadata)

sample_info <- data.frame(
  Sample = c("2RT00374N", "RT00382N", "RT00383N", "RT00385N", "RT00390N"),
  Sample_ID = c("HSB6195", "HSB5871", "HSB8050", "HSB6154", "HSB8073"),
  Age = c(45, 60, 43, 68, 51),
  Sex = c("Male", "Male", "Male", "Female", "Female"),
  stringsAsFactors = FALSE
)

metadata$Cells <- rownames(metadata)
metadata$Dataset <- "GSE207334"
metadata$Technology <- "10X Genomics"
metadata$Sequence <- "snRNA-seq"
metadata$Sample <- metadata$samplename

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

rna_counts <- rna_counts[, metadata$Cells]
object1 <- CreateSeuratObject(
  counts = rna_counts,
  meta.data = metadata
)

log_message("Save data...")

saveRDS(
  object1,
  file.path(res_dir, "GSE207334_processed.rds")
)
