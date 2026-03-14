source("functions/prepare_env.R")

data_dir <- "../../data/BrainOmicsData/raw/GSE199762/"
res_dir <- check_dir("../../data/BrainOmicsData/processed/GSE199762/")

log_message("Start loading data...")

all_files <- list.files(
  data_dir,
  pattern = "\\.(tsv|mtx|csv)\\.gz$", full.names = FALSE
)
counts_files <- grep("_counts\\.mtx\\.gz$", all_files, value = TRUE)

gsm_samples <- gsub("_counts\\.mtx\\.gz$", "", counts_files)
gsm_ids <- gsub("_.*$", "", gsm_samples)
sample_names <- gsub("^[^_]+_", "", gsm_samples)

gsm_groups <- split(sample_names, gsm_ids)
unique_gsm_ids <- unique(gsm_ids)

log_message("Found {.val {length(unique_gsm_ids)}} GSM samples")

read_sample <- function(gsm_id, sample_name) {
  barcodes_file <- file.path(
    data_dir, paste0(gsm_id, "_", sample_name, "_barcodes.tsv.gz")
  )
  genes_file <- file.path(
    data_dir, paste0(gsm_id, "_", sample_name, "_genes.tsv.gz")
  )
  counts_file <- file.path(
    data_dir, paste0(gsm_id, "_", sample_name, "_counts.mtx.gz")
  )
  metadata_file <- file.path(
    data_dir, paste0(gsm_id, "_", sample_name, "_metadata.csv.gz")
  )
  counts <- Matrix::readMM(gzfile(counts_file))
  barcodes <- read.table(
    gzfile(barcodes_file),
    header = FALSE, stringsAsFactors = FALSE
  )[, 1]
  genes <- read.table(
    gzfile(genes_file),
    header = FALSE, stringsAsFactors = FALSE
  )[, 1]
  rownames(counts) <- genes
  colnames(counts) <- barcodes
  metadata <- NULL
  if (file.exists(metadata_file)) {
    metadata <- read.csv(
      gzfile(metadata_file),
      header = TRUE, stringsAsFactors = FALSE, row.names = 1
    )
    common_cells <- intersect(colnames(counts), rownames(metadata))
    counts <- counts[, common_cells]
    metadata <- metadata[common_cells, , drop = FALSE]
  }
  return(
    list(counts = counts, metadata = metadata)
  )
}

gsm_objects <- thisutils::parallelize_fun(
  unique_gsm_ids,
  fun = function(gsm_id) {
    log_message("Processing GSM ID: {.val {gsm_id}}")
    sample_names_gsm <- gsm_groups[[gsm_id]]
    sample_data_list <- list()
    for (sample_name in sample_names_gsm) {
      log_message("  Loading sample: {.val {sample_name}}")
      sample_data <- read_sample(gsm_id, sample_name)
      sample_data_list[[sample_name]] <- sample_data
    }
    all_counts_list <- lapply(
      sample_data_list, function(x) x$counts
    )
    common_genes <- Reduce(
      intersect, lapply(all_counts_list, rownames)
    )

    all_counts_list <- lapply(
      all_counts_list, function(x) x[common_genes, , drop = FALSE]
    )
    combined_counts <- do.call(cbind, all_counts_list)

    all_metadata_list <- lapply(
      sample_data_list, function(x) x$metadata
    )
    all_metadata_list <- all_metadata_list[!sapply(all_metadata_list, is.null)]

    sample_name_vec <- rep(
      sample_names_gsm,
      times = sapply(all_counts_list, ncol)
    )
    names(sample_name_vec) <- colnames(combined_counts)
    for (i in seq_along(all_metadata_list)) {
      if (!is.null(all_metadata_list[[i]])) {
        all_metadata_list[[i]]$sample_name <- names(all_metadata_list)[i]
      }
    }
    combined_metadata <- do.call(rbind, all_metadata_list)
    rownames(combined_metadata) <- colnames(combined_counts)
    combined_metadata$sample_name <- sample_name_vec[rownames(combined_metadata)]
    CreateSeuratObject(
      counts = combined_counts,
      project = gsm_id,
      meta.data = combined_metadata
    )
  },
  cores = 10
)

log_message("Merging all GSM samples into one Seurat object...")
object <- merge(
  gsm_objects[[1]],
  y = gsm_objects[-1],
  project = "GSE199762"
)

object <- JoinLayers(object)

metadata <- object@meta.data
metadata$Cells <- rownames(metadata)
metadata$Dataset <- "GSE199762"
metadata$Technology <- "10X Genomics"
metadata$Sequence <- "snRNA-seq"
metadata$Sample <- metadata$donor
metadata$Sample_ID <- metadata$orig.ident
metadata$CellType_raw <- metadata$all.exp_type
metadata$Brain_Region <- metadata$sample
metadata$Region <- metadata$region

metadata$age <- gsub("GW", " PCW", metadata$age)
days_to_years <- function(age_str) {
  if (grepl("d$", age_str)) {
    days <- as.numeric(gsub("d$", "", age_str))
    years <- round(days / 365, 2)
    return(paste0(years, "y"))
  }
  return(age_str)
}
metadata$age <- sapply(metadata$age, days_to_years)
metadata$age <- gsub("y$", "", metadata$age)
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
  file.path(res_dir, "GSE199762_processed.rds")
)
