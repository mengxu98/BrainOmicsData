source("functions/prepare_env.R")

data_dir <- "../../data/BrainOmicsData/raw/ROSMAP"
res_dir <- check_dir("../../data/BrainOmicsData/processed/ROSMAP/")

if (!file.exists(file.path(data_dir, "individual_metadata_deidentified.tsv"))) {
  download.file(
    "https://personal.broadinstitute.org/cboix/ad427_data/Data/Metadata/individual_metadata_deidentified.tsv",
    file.path(data_dir, "individual_metadata_deidentified.tsv")
  )
}

thisutils::log_message("Start loading data...")

metadata1 <- read.csv(
  file.path(
    data_dir, "individual_metadata_deidentified.tsv"
  ),
  sep = "\t"
)

metadata1$Sex_Source_Raw <- as.character(metadata1$msex)
metadata1$Sex_Source_Standardized <- ifelse(
  metadata1$msex == "0",
  "Female",
  ifelse(metadata1$msex == "1", "Male", NA_character_)
)
metadata1$Sex_Source_Column <- "msex"
metadata1$Sex_Assignment_Method <-
  "source reported code; 0=Female and 1=Male"
metadata1$Sex <- metadata1$Sex_Source_Standardized

metadata2 <- read.csv(
  file.path(
    data_dir, "RNA/meta.tsv"
  ),
  sep = "\t"
)

if (nrow(metadata2) != 414964L ||
  sum(metadata2$ADdiag3types == "nonAD") != 213472L ||
  length(unique(metadata2$Individual)) != 92L) {
  stop("ROSMAP complete RNA metadata design differs from source audit")
}

metadata <- merge(
  metadata2,
  metadata1,
  by.x = "Individual",
  by.y = "subject",
  all.x = TRUE
)
metadata$Diagnosis_RNA_Group_raw <- metadata$ADdiag3types
metadata$Diagnosis_Pathology_raw <- metadata$Pathologic_diagnosis_of_AD
metadata$Diagnosis_raw <- paste0(
  "ADdiag3types=", metadata$Diagnosis_RNA_Group_raw,
  "; Pathologic_diagnosis_of_AD=", metadata$Diagnosis_Pathology_raw
)
metadata$Diagnosis_Inclusion_Rule <- paste(
  "ADdiag3types == nonAD AND Pathologic_diagnosis_of_AD == no"
)
metadata$Diagnosis_Inclusion_Status <- ifelse(
  metadata$ADdiag3types == "nonAD" &
    metadata$Pathologic_diagnosis_of_AD == "no",
  "retained non-AD/pathology-negative",
  "excluded"
)
# Preserve the project-reader selection that precedes the processed object.
source_diagnosis <- aggregate(
  rep(1L, nrow(metadata)),
  metadata[c(
    "Individual", "Diagnosis_RNA_Group_raw", "Diagnosis_Pathology_raw",
    "Diagnosis_Inclusion_Rule", "Diagnosis_Inclusion_Status"
  )],
  sum
)
names(source_diagnosis)[names(source_diagnosis) == "x"] <- "Source_Cells"
source_diagnosis$Dataset <- "ROSMAP"
if (nrow(source_diagnosis) != 92L ||
  sum(source_diagnosis$Source_Cells) != 414964L) {
  stop("ROSMAP source diagnosis donor ledger failed conservation checks")
}
write.table(
  source_diagnosis,
  file.path(res_dir, "source_diagnosis_selection.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
metadata <- metadata[
  metadata$Diagnosis_Inclusion_Status ==
    "retained non-AD/pathology-negative", ,
  drop = FALSE
]
if (nrow(metadata) != 172683L ||
  length(unique(metadata$Individual)) != 40L ||
  length(unique(metadata$Batch)) != 4L) {
  stop("ROSMAP non-AD/pathology-negative retained design differs from audit")
}
batch_donors <- vapply(
  split(metadata$Individual, metadata$Batch),
  function(value) length(unique(value)),
  integer(1)
)
if (!identical(sort(unname(batch_donors)), c(8L, 8L, 11L, 13L))) {
  stop("ROSMAP retained Batch-to-donor design differs from audit")
}

metadata$Age_Source_Raw <- metadata$age_death
metadata$Age_Harmonization_Input <- ifelse(
  metadata$age_death == "90+",
  "90+ years",
  paste0(
    sub("^\\(([^,]+),([^]]+)\\]$", "\\1-\\2", metadata$age_death),
    " years"
  )
)
metadata$Age <- metadata$Age_Source_Raw
metadata$Age_Source_Unit <- "reported age-at-death interval in years"
metadata$Age_Source_Basis <- "postnatal age at death"
metadata$Age_Source_Reference <- "individual_metadata_deidentified.tsv:age_death"
metadata$Age_Conversion_Formula <- ifelse(
  metadata$age_death == "90+",
  "retain 90-year lower bound and open upper bound",
  "retain reported lower and upper bounds; midpoint is representative"
)
metadata$Age_Conversion_Confidence <- ifelse(
  metadata$age_death == "90+",
  "open-ended reported range",
  "reported interval"
)
metadata$Age_Conversion_Applied <- TRUE

rownames(metadata) <- metadata$cellId
metadata$Cells <- metadata$cellId
metadata$Dataset <- "ROSMAP"
metadata$Technology <- "10X Genomics"
metadata$Sequence <- "snRNA-seq"
metadata$Sample <- metadata$Individual
metadata$Sample_ID <- metadata$Individual
metadata$Original_Donor_ID <- metadata$Individual
metadata$Original_Specimen_ID <- metadata$Individual
metadata$Original_Source_Record_ID <- metadata$Individual
metadata$CellType_raw <- metadata$Celltype
metadata$Brain_Region <- "Prefrontal cortex"
metadata$Region <- "Prefrontal cortex"
metadata$Assay_Type <- "single-nucleus RNA-seq"
metadata$Sequencing_Platform <- paste(
  "NextSeq 500/550 or NovaSeq 6000;",
  "library-level platform mapping unavailable"
)
metadata$Library_Chemistry <- "10x Chromium Single Cell 3-prime v3"
metadata$Original_File_Type <- "10x Genomics Matrix Market gene-count matrix"
metadata$Expression_Units <- "raw integer UMI counts including pre-mRNA"
metadata$Genome_Transcriptome_Build <- "GRCh38"
metadata$Quantification_Software <- "Cell Ranger 3.0.2"
metadata$Raw_Integer_Counts_Available <- TRUE
metadata$Source_Full_RNA_Nuclei <- 414964L
metadata$Source_RNA_nonAD_Nuclei <- 213472L
metadata$Source_Retained_Nuclei <- 172683L
metadata$Source_Excluded_Pathology_Positive_Nuclei <- 40789L
metadata$Source_Preprocessing_Retained_Object_Scope <- paste(
  "complete released QC-passed RNA count subset satisfying both",
  "ADdiag3types=nonAD and Pathologic_diagnosis_of_AD=no"
)
metadata$Source_Preprocessing_Attrition_Status <- paste(
  "414,964 released RNA nuclei; 213,472 are labelled nonAD in RNA/meta.tsv;",
  "40,789 nuclei from eight of those donors have pathology=yes in the",
  "individual table; 172,683 nuclei from 40 concordant non-AD/pathology-",
  "negative donors are retained"
)
metadata$Source_Preprocessing_Diagnosis_Audit_Status <- paste(
  "verified by exact intersection of RNA/meta.tsv ADdiag3types and",
  "individual_metadata_deidentified.tsv Pathologic_diagnosis_of_AD"
)

column_order <- c(
  "Cells", "Dataset", "Technology", "Sequence", "Sample",
  "Sample_ID", "Original_Donor_ID", "Original_Specimen_ID",
  "Original_Source_Record_ID", "CellType_raw", "Celltype", "Subcelltype",
  "Batch", "Brain_Region", "Region", "Age", "Sex"
)
metadata <- retain_source_metadata(metadata, column_order)

counts <- Read10X(
  file.path(
    data_dir, "RNA"
  ),
  gene.column = 1
)

intersect_cells <- intersect(colnames(counts), rownames(metadata))
if (length(intersect_cells) != 172683L ||
  !setequal(intersect_cells, rownames(metadata))) {
  stop("ROSMAP retained metadata does not match exactly 172,683 count columns")
}
counts <- counts[, intersect_cells]

metadata <- metadata[intersect_cells, ]

object <- CreateSeuratObject(
  counts = counts,
  meta.data = metadata
)

thisutils::log_message("Save data...")
saveRDS(
  object,
  file.path(res_dir, "ROSMAP_processed.rds")
)
