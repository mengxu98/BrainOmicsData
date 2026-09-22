source("functions/prepare_env.R")

data_dir <- "../../data/BrainOmicsData/raw/GSE212606"
res_dir <- check_dir("../../data/BrainOmicsData/processed/GSE212606/")

thisutils::log_message("Start loading data...")
counts <- Matrix::readMM(
  file.path(data_dir, "GSM6657986_gene_count.txt.gz")
)
gene_meta <- read.csv(
  file.path(data_dir, "GSM6657986_gene_annotation.csv")
)

thisutils::log_message("Counts matrix dimensions: {.val {dim(counts)}}")
thisutils::log_message("Gene metadata dimensions: {.val {dim(gene_meta)}}")

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
source_cells <- nrow(metadata)
source_ad_cells <- sum(metadata$Condition == "AD", na.rm = TRUE)
source_wt_smtg_missing_donor_cells <- sum(
  metadata$Condition == "WT" &
    metadata$Region == "SMTG" &
    is.na(metadata$Individual_ID)
)
metadata <- metadata[
  metadata$Condition == "WT" & !is.na(metadata$Individual_ID), ,
  drop = FALSE
]
counts <- counts[, rownames(metadata), drop = FALSE]

sample_info <- data.frame(
  Individual_ID = c("5459", "5356", "1311", "1306", "1304", "1247"),
  Age = c(70, 94, 85, 83, 81, 94),
  Sex = c("Female", "Male", "Female", "Female", "Male", "Male"),
  stringsAsFactors = FALSE
)
sample_info$Sex_Source_Raw <- sample_info$Sex
sample_info$Sex_Source_Standardized <- sample_info$Sex
sample_info$Sex_Source_Column <- "source publication donor table"
sample_info$Sex_Assignment_Method <-
  "curated from source publication donor table; not inferred"
metadata <- merge(
  metadata,
  sample_info,
  by = "Individual_ID",
  sort = FALSE
)

rownames(metadata) <- metadata$Cells
metadata <- metadata[colnames(counts), , drop = FALSE]
if (!identical(rownames(metadata), colnames(counts))) {
  stop("GSE212606 metadata and count matrix cell order differ")
}
if (nrow(metadata) != 30801L || source_cells != 118240L ||
  source_ad_cells != 50658L ||
  source_wt_smtg_missing_donor_cells != 36781L) {
  stop("GSE212606 source or retained attrition counts differ from audit")
}

metadata$EasySci_PCR_Group_ID <- sub(
  "[.][^.]+$",
  "",
  metadata$Cells
)
cell_barcode <- sub("^.*[.]", "", metadata$Cells)
if (any(!grepl("^[ACGT]{20}$", cell_barcode)) ||
  any(!grepl(
    "^Hippocampus_[12]_[0-9]{2}$",
    metadata$EasySci_PCR_Group_ID
  ))) {
  stop("GSE212606 Cell_ID does not match the audited EasySci PCR-group pattern")
}
pcr_group_donors <- tapply(
  as.character(metadata$Individual_ID),
  metadata$EasySci_PCR_Group_ID,
  function(value) length(unique(value))
)
if (length(pcr_group_donors) != 120L ||
  any(pcr_group_donors != 6L)) {
  stop("GSE212606 retained PCR groups do not each span all six control donors")
}

metadata$Dataset <- "GSE212606"
metadata$Technology <- "EasySci-RNA"
metadata$Sequence <- "snRNA-seq"
metadata$Sequencing_Platform <- "Illumina NovaSeq 6000"
metadata$Library_Chemistry <- "EasySci-RNA combinatorial indexing workflow"
metadata$Sample <- metadata$Individual_ID
metadata$Sample_ID <- metadata$Individual_ID
metadata$CellType_raw <- metadata$Cell_type
metadata$Brain_Region <- metadata$Region
metadata$Diagnosis_raw <- metadata$Condition
metadata$Source_Preprocessing_Retained_Object_Scope <- paste(
  "complete 30,801-cell WT hippocampus cohort with source-reported",
  "Individual_ID from six control donors"
)
metadata$Source_Preprocessing_Attrition_Status <- paste(
  "source processed matrix: 118,240 cells; excluded 50,658 AD cells;",
  "excluded 36,781 WT SMTG cells because source Individual_ID is missing;",
  "retained 30,801 WT hippocampus cells"
)
metadata$Source_Preprocessing_Diagnosis_Audit_Status <- paste(
  "source Condition retained as WT for all 30,801 cells;",
  "all 50,658 source AD cells excluded before object construction"
)
column_order <- c(
  "Cells", "Dataset", "Technology", "Sequence", "Sample",
  "Sample_ID", "CellType_raw", "Brain_Region", "Region", "Age", "Sex",
  "Diagnosis_raw", "EasySci_PCR_Group_ID"
)
metadata <- retain_source_metadata(metadata, column_order)

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

thisutils::log_message("Save data...")
saveRDS(
  object,
  file.path(res_dir, "GSE212606_processed.rds")
)
