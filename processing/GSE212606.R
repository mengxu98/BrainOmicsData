source("functions/prepare_env.R")

data_dir <- "../../data/BrainOmicsData/raw/GSE212606"
res_dir <- check_dir("../../data/BrainOmicsData/processed/GSE212606/")

log_message("Start loading data...")
counts <- Matrix::readMM(
  file.path(data_dir, "GSM6657986_gene_count.txt.gz")
)
gene_meta <- read.csv(
  file.path(data_dir, "GSM6657986_gene_annotation.csv")
)

log_message("Counts matrix dimensions: {.val {dim(counts)}}")
log_message("Gene metadata dimensions: {.val {dim(gene_meta)}}")

gene_meta$gene_name_unique <- ifelse(
  duplicated(gene_meta$gene_short_name) | duplicated(gene_meta$gene_short_name, fromLast = TRUE),
  paste0(gene_meta$gene_short_name, "_", gene_meta$gene_id),
  gene_meta$gene_short_name
)

rownames(counts) <- gene_meta$gene_name_unique

metadata <- read.csv(
  file.path(data_dir, "GSM6657986_cell_annotation.csv"),
  header = TRUE,
  row.names = 1
)
colnames(counts) <- rownames(metadata)
metadata$Cells <- rownames(metadata)
metadata <- metadata[metadata$Condition == "WT", ]
counts <- counts[, rownames(metadata)]

sample_info <- data.frame(
  Individual_ID = c("5459", "5356", "1311", "1306", "1304", "1247"),
  Age = c(70, 94, 85, 83, 81, 94),
  Sex = c("Female", "Male", "Female", "Female", "Male", "Male"),
  stringsAsFactors = FALSE
)
metadata <- merge(metadata, sample_info, by = "Individual_ID")

rownames(metadata) <- metadata$Cells
metadata$Dataset <- "GSE212606"
metadata$Technology <- "EasySci-RNA"
metadata$Sequence <- "scRNA-seq"
metadata$Sample <- metadata$Individual_ID
metadata$Sample_ID <- metadata$Individual_ID
metadata$CellType_raw <- metadata$Cell_type
metadata$Brain_Region <- metadata$Region
column_order <- c(
  "Cells", "Dataset", "Technology", "Sequence", "Sample",
  "Sample_ID", "CellType_raw", "Brain_Region", "Region", "Age", "Sex"
)
metadata <- metadata[, column_order]
metadata <- na.omit(metadata)

common_cells <- intersect(
  rownames(metadata),
  colnames(counts)
)
metadata <- metadata[common_cells, ]
counts <- counts[, common_cells]
object <- CreateSeuratObject(
  counts = counts,
  meta.data = metadata
)

log_message("Save data...")
saveRDS(
  object,
  file.path(res_dir, "GSE212606_processed.rds")
)
