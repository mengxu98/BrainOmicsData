source("functions/prepare_env.R")

data_dir <- "../../data/BrainOmicsData/raw/GSE81475/"
res_dir <- check_dir("../../data/BrainOmicsData/processed/GSE81475/")

log_message("Start loading data...")
counts <- read.table(
  file.path(data_dir, "counts.txt"),
  header = TRUE,
  row.names = 1
)

metadata <- read.csv(
  file.path(res_dir, "metadata.csv"),
  row.names = 1
)

common_cells <- intersect(rownames(metadata), colnames(counts))
counts <- counts[, common_cells, drop = FALSE]

log_message("Processing gene names and aggregating duplicates...")
gene_names <- sub("^[^|]+\\|", "", rownames(counts))

counts_matrix <- as.matrix(counts)
gene_counts <- table(gene_names)
counts_sum <- rowsum(
  counts_matrix,
  group = gene_names, reorder = FALSE
)
gene_counts_vec <- as.numeric(gene_counts[rownames(counts_sum)])
counts <- sweep(counts_sum, 1, gene_counts_vec, "/")

metadata <- metadata[common_cells, , drop = FALSE]

metadata$Dataset <- "GSE81475"
metadata$Technology <- "Fluidigm C1"
metadata$Sequence <- "scRNA-seq"
metadata$Sample <- metadata$Cells
metadata$Sample_ID <- metadata$Cells
metadata$Brain_Region <- metadata$Region <- "Neocortex"
metadata$CellType_raw <- metadata$CellType
column_order <- c(
  "Cells", "Dataset", "Technology", "Sequence", "Sample",
  "Sample_ID", "CellType_raw", "Brain_Region", "Region", "Age", "Sex"
)
metadata <- metadata[, column_order]
metadata <- na.omit(metadata)
counts <- counts[, metadata$Cells]

log_message("Creating Seurat object...")
object <- CreateSeuratObject(
  counts = counts,
  meta.data = metadata
)

log_message("Save data...")
saveRDS(
  object,
  file.path(res_dir, "GSE81475_processed.rds")
)
