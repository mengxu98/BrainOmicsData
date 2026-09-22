source("functions/prepare_env.R")

data_dir <- "../../data/BrainOmicsData/raw/GSE207334"
res_dir <- check_dir("../../data/BrainOmicsData/processed/GSE207334/")

thisutils::log_message("Start loading data...")

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
  Source_GEO_RNA_Record_ID = c(
    "GSM6284664", "GSM6284665", "GSM6284666", "GSM6284667",
    "GSM6284668"
  ),
  Age = c(45, 60, 43, 68, 51),
  Sex = c("Male", "Male", "Male", "Female", "Female"),
  stringsAsFactors = FALSE
)
sample_info$Sex_Source_Raw <- sample_info$Sex
sample_info$Sex_Source_Standardized <- sample_info$Sex
sample_info$Sex_Source_Column <- "source publication donor table"
sample_info$Sex_Assignment_Method <-
  "curated from source publication donor table; not inferred"

metadata$Cells <- rownames(metadata)
metadata$Dataset <- "GSE207334"
metadata$Technology <- "10X Genomics"
metadata$Sequence <- "snRNA-seq"
metadata$Sample <- metadata$samplename

metadata <- merge(
  metadata, sample_info,
  by = "Sample", all.x = TRUE, sort = FALSE
)
rownames(metadata) <- metadata$Cells
metadata <- metadata[colnames(rna_counts), , drop = FALSE]
if (!identical(rownames(metadata), colnames(rna_counts)) ||
  anyNA(metadata$Source_GEO_RNA_Record_ID)) {
  stop("GSE207334 RNA library records do not cover every retained cell")
}

metadata$CellType_raw <- metadata$subclass
metadata$Brain_Region <- "Dorsolateral prefrontal cortex"
metadata$Region <- "Dorsolateral prefrontal cortex"

column_order <- c(
  "Cells", "Dataset", "Technology", "Sequence", "Sample",
  "Sample_ID", "Source_GEO_RNA_Record_ID", "CellType_raw",
  "Brain_Region", "Region", "Age", "Sex"
)
metadata <- retain_source_metadata(metadata, column_order)

rna_counts <- rna_counts[, metadata$Cells]
object1 <- CreateSeuratObject(
  counts = rna_counts,
  meta.data = metadata
)

thisutils::log_message("Save data...")

saveRDS(
  object1,
  file.path(res_dir, "GSE207334_processed.rds")
)
