source("functions/prepare_env.R")

data_dir <- "../../data/BrainOmicsData/raw/GSE204683"
res_dir <- check_dir("../../data/BrainOmicsData/processed/GSE204683/")

log_message("Start loading data...")
counts <- readRDS(
  file.path(data_dir, "GSE204683_count_matrix.RDS")
)

metadata <- read.csv(
  file.path(data_dir, "GSE204683_barcodes.tsv"),
  header = TRUE,
  row.names = NULL,
  sep = "\t"
)

donor_id_mapping <- c(
  "4" = "LaFet1", "8" = "LaFet2", "11" = "EaFet1", "16" = "EaFet2",
  "150656" = "Adult2", "150666" = "Adult1", "4413" = "Inf1", "4422" = "Inf2",
  "5936" = "Adol2", "5977" = "Child2", "6007" = "Adol1", "6032" = "Child1"
)

donor_id_reverse_mapping <- names(donor_id_mapping)
names(donor_id_reverse_mapping) <- donor_id_mapping

metadata$Donor_ID_numeric <- donor_id_reverse_mapping[metadata$Donor.ID]
metadata$Cells <- paste0(metadata$Donor_ID_numeric, "_", metadata$Barcode)

counts_cells <- colnames(counts)
metadata <- metadata[metadata$Cells %in% counts_cells, ]

missing_cells <- setdiff(counts_cells, metadata$Cells)
if (length(missing_cells) > 0) {
  log_message(
    "Found {.val {length(missing_cells)}} cells in counts but not in metadata, creating metadata for them..."
  )
  donor_ids_missing <- sub("_.*", "", missing_cells)
  barcodes_missing <- sub("^[^_]+_", "", missing_cells)
  missing_meta <- data.frame(
    Barcode = barcodes_missing,
    Donor.ID = donor_id_mapping[donor_ids_missing],
    Cell.type = NA_character_,
    UMI.counts = NA_real_,
    gene.count = NA_real_,
    raw.reads.count = NA_real_,
    mapped.reads.count = NA_real_,
    Donor_ID_numeric = donor_ids_missing,
    Cells = missing_cells,
    stringsAsFactors = FALSE,
    row.names = missing_cells
  )
  metadata <- rbind(metadata, missing_meta)
}

rownames(metadata) <- metadata$Cells
metadata <- metadata[counts_cells, ]

sample_info <- data.frame(
  Sample_Name = c(
    "Adol1", "Child1", "EaFet2", "Adult2", "Child2", "Adult1",
    "Inf1", "LaFet1", "LaFet2", "Adol2", "Inf2", "EaFet1"
  ),
  Brain_Region = c(
    "BA 9/46", "BA 9/46", "cortical plate",
    "BM_9/10/46", "BA 9/46", "BM_9/10/46",
    "BA 9/46", "cortical plate", "cortical plate",
    "BA 9/46", "BA 9/46", "cortical plate"
  ),
  Region = c(
    "Cerebral cortex", "Cerebral cortex", "Cerebral cortex",
    "Cerebral cortex", "Cerebral cortex", "Cerebral cortex",
    "Cerebral cortex", "Cerebral cortex", "Cerebral cortex",
    "Cerebral cortex", "Cerebral cortex", "Cerebral cortex"
  ),
  Age = c(
    "14", "4", "19 PCW", "39", "6", "20", "1",
    "23 PCW", "24 PCW", "14", "1", "18 PCW"
  ),
  Sex = c(
    "Male", "Male", "Female", "Male", "Female", "Female",
    "Female", "Male", "Female", "Female", "Male", "Female"
  ),
  stringsAsFactors = FALSE
)

metadata <- merge(
  metadata, sample_info,
  by.x = "Donor.ID", by.y = "Sample_Name", all.x = TRUE
)
rownames(metadata) <- metadata$Cells
common_cells <- intersect(rownames(metadata), colnames(counts))
metadata$Dataset <- "GSE204683"
metadata$Technology <- "10X Genomics"
metadata$Sequence <- "snRNA-seq"
metadata$Sample <- metadata$Donor.ID
metadata$Sample_ID <- metadata$Donor_ID_numeric
metadata$CellType_raw <- metadata$Cell.type
column_order <- c(
  "Cells", "Dataset", "Technology", "Sequence", "Sample",
  "Sample_ID", "CellType_raw", "Brain_Region", "Region", "Age", "Sex"
)
metadata <- metadata[common_cells, column_order]
counts <- counts[, common_cells]
object <- CreateSeuratObject(
  counts = counts,
  meta.data = metadata
)

log_message("Save data...")
saveRDS(
  object,
  file.path(res_dir, "GSE204683_processed.rds")
)
