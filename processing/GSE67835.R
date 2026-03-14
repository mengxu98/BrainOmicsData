source("functions/prepare_env.R")

data_dir <- "../../data/BrainOmicsData/raw/GSE67835/GSE67835"
res_dir <- check_dir("../../data/BrainOmicsData/processed/GSE67835/")

log_message("Start loading data...")

csv_files <- list.files(
  data_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)
log_message("Found {.val {length(csv_files)}} CSV files")

read_cell_data <- function(csv_file) {
  cell_name <- gsub(
    "^(GSM[0-9]+)_.*\\.csv$", "\\1", basename(csv_file)
  )
  cell_data <- read.table(
    csv_file,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    col.names = c("Gene", "Expression")
  )

  cell_data$Gene <- trimws(cell_data$Gene)

  expr_vec <- cell_data$Expression
  names(expr_vec) <- cell_data$Gene

  return(
    list(cell_name = cell_name, expression = expr_vec)
  )
}

log_message("Reading CSV files...")
cell_data_list <- thisutils::parallelize_fun(
  csv_files,
  fun = read_cell_data,
  cores = 10
)

all_genes <- unique(
  unlist(lapply(cell_data_list, function(x) names(x$expression)))
)
log_message("Found {.val {length(all_genes)}} unique genes")

expr_matrix <- matrix(
  0,
  nrow = length(all_genes),
  ncol = length(cell_data_list),
  dimnames = list(all_genes, NULL)
)

cell_names <- character(length(cell_data_list))
for (i in seq_along(cell_data_list)) {
  cell_info <- cell_data_list[[i]]
  cell_names[i] <- cell_info$cell_name
  expr_matrix[names(cell_info$expression), i] <- cell_info$expression
}

colnames(expr_matrix) <- cell_names
log_message("Creating Seurat object...")
object <- CreateSeuratObject(
  counts = expr_matrix,
  project = "GSE67835"
)

log_message("Save data...")
saveRDS(
  object,
  file.path(res_dir, "GSE67835.rds")
)

object <- readRDS(
  file.path(res_dir, "GSE67835.rds")
)

metadata <- read.csv(
  file.path(res_dir, "metadata.csv"),
  row.names = 1
)
metadata$Cells <- rownames(metadata)
metadata$Dataset <- "GSE67835"
metadata$Technology <- "Fluidigm C1"
metadata$Sequence <- "scRNA-seq"
metadata$Sample <- rownames(metadata)
metadata$Sample_ID <- rownames(metadata)
metadata$Brain_Region <- metadata$Region
metadata$Region <- metadata$Region
metadata$CellType_raw <- metadata$cluster
region_full <- c(
  "CTX" = "Cerebral cortex",
  "ITC" = "Intra-temporal cortex"
)

metadata$Region_fullname <- region_full[metadata$Region]
metadata$Brain_Region <- metadata$Region_fullname
metadata$Region <- metadata$Region_fullname


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
  file.path(res_dir, "GSE67835_processed.rds")
)
