source("functions/prepare_env.R")

data_dir <- "../../data/BrainOmicsData/raw/GSE261983/"
res_dir <- check_dir("../../data/BrainOmicsData/processed/GSE261983/")

log_message("Start loading data...")

sample_dirs <- list.dirs(data_dir, full.names = FALSE, recursive = FALSE)
rna_dirs <- sample_dirs[grepl("_RNA$", sample_dirs)]
rna_dirs <- rna_dirs[grepl("^GSM", rna_dirs)]

delete_samples <- c(
  "RT00373N", "RT00374N", "RT00382N",
  "RT00383N", "RT00385N", "RT00389N", "RT00390N"
)
delete_indices <- which(substring(rna_dirs, 12, 19) %in% delete_samples)
rna_dirs <- rna_dirs[-delete_indices]
samples_info <- read.csv(
  file.path(data_dir, "samples_info.csv"),
  header = TRUE, stringsAsFactors = FALSE
)
gene_list <- read.table(
  file.path(data_dir, "GSE261983_RNA_gene_list.txt.gz"),
  header = FALSE, stringsAsFactors = FALSE
)[, 1]

log_message("Found {.val {length(rna_dirs)}} RNA sample directories")
objects_list <- thisutils::parallelize_fun(
  rna_dirs,
  fun = function(sample_dir) {
    sample_path <- file.path(data_dir, sample_dir)
    counts <- Matrix::readMM(file.path(sample_path, "matrix.mtx.gz"))
    cells <- readLines(gzfile(file.path(sample_path, "barcodes.tsv.gz")))
    cells <- gsub("#", "_", cells)
    rownames(counts) <- gene_list
    colnames(counts) <- cells
    object <- CreateSeuratObject(
      counts = counts,
      project = sample_dir
    )
    object$Sample_ID <- gsub("_RNA$", "", sample_dir)
    object$Sample <- substring(sample_dir, 12, 19)
    return(object)
  },
  cores = length(rna_dirs)
)

log_message("Merging {.val {length(objects_list)}} samples...")

object <- merge(
  objects_list[[1]],
  y = objects_list[-1],
  project = "GSE261983"
)
object <- JoinLayers(object)

metadata <- object@meta.data
metadata$Cells <- rownames(metadata)
metadata$Dataset <- "GSE261983"
metadata$Technology <- "10X Genomics"
metadata$Sequence <- "snRNA-seq"

metadata <- merge(
  metadata,
  samples_info,
  by.x = "Sample",
  by.y = "Individual.ID",
  all.x = TRUE
)
rownames(metadata) <- metadata$Cells

metadata$CellType_raw <- "Unknown"
metadata$Brain_Region <- "Prefrontal Cortex"
metadata$Region <- "Prefrontal Cortex"

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
  file.path(res_dir, "GSE261983_processed.rds")
)
