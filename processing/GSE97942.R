source("functions/prepare_env.R")

data_dir <- "../../data/BrainOmicsData/raw/GSE97942/"
res_dir <- check_dir("../../data/BrainOmicsData/processed/GSE97942/")

log_message("Start loading data...")
counts_cerebellar_hemisphere <- read.table(
  file.path(data_dir, "CerebellarHem_counts.txt"),
  header = TRUE,
  row.names = 1,
  sep = "\t"
)
counts_frontal_cortex <- read.table(
  file.path(data_dir, "FrontalCortex_counts.txt"),
  header = TRUE,
  row.names = 1,
  sep = "\t"
)

counts_visual_cortex <- read.table(
  file.path(data_dir, "VisualCortex_counts.txt"),
  header = TRUE,
  row.names = 1,
  sep = "\t"
)
common_genes <- Reduce(
  intersect,
  list(
    rownames(counts_cerebellar_hemisphere),
    rownames(counts_frontal_cortex),
    rownames(counts_visual_cortex)
  )
)
counts_cerebellar_hemisphere <- counts_cerebellar_hemisphere[common_genes, ]
counts_frontal_cortex <- counts_frontal_cortex[common_genes, ]
counts_visual_cortex <- counts_visual_cortex[common_genes, ]
counts <- Reduce(
  cbind,
  list(
    counts_cerebellar_hemisphere,
    counts_frontal_cortex,
    counts_visual_cortex
  )
)
rm(
  counts_cerebellar_hemisphere,
  counts_frontal_cortex,
  counts_visual_cortex
)

original_colnames <- colnames(counts)
new_colnames <- sub("^[^_]+_", "", original_colnames)
colnames(counts) <- new_colnames

metadata <- read.csv(
  file.path(data_dir, "metadata.csv"),
  row.names = 1
)
common_cells <- intersect(colnames(counts), rownames(metadata))

counts <- counts[, common_cells, drop = FALSE]
metadata <- metadata[common_cells, , drop = FALSE]

metadata$Cells <- rownames(metadata)
metadata$Dataset <- "GSE97942"
metadata$Technology <- metadata$protocal
metadata$Sequence <- "snRNA-seq"
metadata$Sample <- metadata$orig.ident
metadata$Sample <- gsub("D7_", "", metadata$Sample)
metadata$Sample_ID <- metadata$Sample
metadata$CellType_raw <- metadata$priCluster
metadata$Brain_Region <- metadata$Area
metadata$Region <- metadata$Area
metadata$Age <- gsub(" years", "", metadata$Age)

column_order <- c(
  "Cells", "Dataset", "Technology", "Sequence", "Sample",
  "Sample_ID", "CellType_raw", "Brain_Region", "Region", "Age", "Sex"
)
metadata <- metadata[, column_order]
metadata <- na.omit(metadata)

counts <- counts[, metadata$Cells]
object <- CreateSeuratObject(
  counts = counts,
  meta.data = metadata
)

log_message("Save data...")
saveRDS(
  object,
  file.path(res_dir, "GSE97942_processed.rds")
)
