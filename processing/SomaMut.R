source("functions/prepare_env.R")

data_dir <- brainomics_data_path("raw/SomaMut")
res_dir <- check_dir(brainomics_data_path("processed/SomaMut/"))

thisutils::log_message("Start loading data...")

object <- readRDS(
  file.path(
    data_dir, "pfc.clean.rds"
  )
)

metadata <- object@meta.data
metadata$Cells <- rownames(metadata)
metadata$Dataset <- "SomaMut"
metadata$Technology <- "10X Genomics"
metadata$Sequence <- "snRNA-seq"
metadata$Sample <- as.character(metadata$orig.ident)
metadata$Sample_ID <- metadata$case
metadata$Original_Donor_ID <- as.character(metadata$orig.ident)
metadata$Original_Specimen_ID <- as.character(metadata$orig.ident)
metadata$Original_Library_ID <- metadata$case
metadata$Original_Source_Record_ID <- metadata$case
metadata$CellType_raw <- as.character(metadata$new_clusters)
metadata$Source_Reference_Cell_Type_Label <- metadata$predicted.id
metadata$Source_Reference_Cell_Type_Confidence <-
  metadata$prediction.score.max
metadata$Source_Reference_Cell_Type_Method <- paste(
  "RPCA mapping to Velmeshev et al. reference; public pfc.clean.rds",
  "retains all 17 per-cell prediction.score.* probabilities"
)
metadata$Brain_Region <- "Prefrontal cortex"
metadata$Region <- "Prefrontal cortex"
metadata$Age <- metadata$age
metadata$Age_Source_Raw <- as.character(metadata$age)
metadata$Age_Source_Unit <- "years"
metadata$Age_Source_Basis <- "postnatal age at death"
metadata$Age_Source_Reference <- "public pfc.clean.rds:age"
metadata$Sex_Source_Raw <- as.character(metadata$sex)
metadata$Sex_Assignment_Method <- "source reported"
metadata$Sex <- standardize_source_sex_value(metadata$Sex_Source_Raw)
metadata$Assay_Type <- "single-nucleus RNA-seq"
metadata$Sequencing_Platform <- "Illumina NovaSeq 6000"
metadata$Library_Chemistry <- "10x Next GEM Single Cell 3-prime v3 or v3.1"
metadata$Original_File_Type <- "public processed Seurat RDS with raw counts"
metadata$Expression_Units <- "raw integer UMI counts including introns"
metadata$Genome_Transcriptome_Build <- "GRCh38 / GENCODE v32"
metadata$Quantification_Software <- "Cell Ranger 6.0.2"
metadata$Raw_Integer_Counts_Available <- TRUE
metadata$Diagnosis_raw <-
  "no neurological disease history or neuropathology evidence"
metadata$Diagnosis_Inclusion_Status <- "retained neurologically normal"
metadata$Source_Full_Public_Clean_Nuclei <- 367317L
metadata$Source_Retained_Nuclei <- 367317L
metadata$Source_Donor_Count <- 19L
metadata$Source_Library_Count <- 42L
metadata$Source_Preparation_Batch_Count <- 13L
metadata$Source_Previously_Published_Donor <-
  as.character(metadata$orig.ident) == "1465"
metadata$Source_Preprocessing_Retained_Object_Scope <- paste(
  "complete public pfc.clean.rds analysis object after source-study",
  "sample, cell, doublet/artefact and cluster QC"
)
metadata$Source_Preprocessing_Attrition_Status <- paste(
  "the public pfc.clean.rds contains all 367,317 source-retained nuclei;",
  "pre-QC nuclei counts are not reconstructible from this object alone;",
  "the paper reports removing replicates 5817_200102, 5288_200128 and",
  "5887_PFC_210601 plus ambiguous artefact clusters"
)
metadata$Source_Preprocessing_Diagnosis_Audit_Status <- paste(
  "paper-verified selection by absence of neurological disease history",
  "or neuropathology evidence"
)

if (nrow(metadata) != 367317L ||
  length(unique(as.character(metadata$orig.ident))) != 19L ||
  length(unique(metadata$case)) != 42L ||
  length(unique(metadata$batch)) != 13L ||
  length(unique(metadata$new_clusters3)) != 7L ||
  length(unique(metadata$new_clusters2)) != 13L ||
  length(unique(metadata$new_clusters)) != 31L ||
  length(unique(metadata$predicted.id)) != 17L) {
  stop("SomaMut retained donor/library/batch/annotation design differs")
}
library_design <- unique(metadata[, c("orig.ident", "case", "batch")])
if (nrow(library_design) != 42L || anyDuplicated(library_design$case)) {
  stop("SomaMut source library mapping is not one-to-one")
}
library_donor <- sub("_.*$", "", library_design$case)
if (!identical(library_donor, as.character(library_design$orig.ident))) {
  stop("SomaMut source library prefix differs from donor identity")
}

column_order <- c(
  "Cells", "Dataset", "Technology", "Sequence", "Sample",
  "Sample_ID", "Original_Donor_ID", "Original_Specimen_ID",
  "Original_Library_ID", "Original_Source_Record_ID", "CellType_raw",
  "new_clusters3", "new_clusters2", "new_clusters", "predicted.id",
  "Source_Reference_Cell_Type_Label",
  "Source_Reference_Cell_Type_Confidence",
  "Source_Reference_Cell_Type_Method",
  "batch", "Brain_Region", "Region", "Age", "Sex"
)
metadata <- retain_source_metadata(metadata, column_order)

counts <- GetAssayData(object, layer = "counts")
counts <- counts[, metadata$Cells]
object <- CreateSeuratObject(
  counts = counts,
  meta.data = metadata
)

thisutils::log_message("Save data...")
saveRDS(
  object,
  file.path(res_dir, "SomaMut_processed.rds")
)
