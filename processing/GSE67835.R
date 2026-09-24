source("functions/prepare_env.R")
source("functions/dataset_metadata.R")

data_dir <- brainomics_data_path("raw/GSE67835/GSE67835")
res_dir <- check_dir(brainomics_data_path("processed/GSE67835/"))

thisutils::log_message("Start loading data...")

csv_files <- list.files(
  data_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)
thisutils::log_message("Found {.val {length(csv_files)}} CSV files")
if (length(csv_files) != 466L) {
  stop("GSE67835 requires all 466 published single-cell count files")
}

parse_geo_soft_samples <- function(soft_file) {
  lines <- readLines(soft_file, warn = FALSE)
  starts <- grep("^\\^SAMPLE = ", lines)
  if (length(starts) == 0L) {
    stop("GSE67835 GEO SOFT file contains no sample records")
  }
  ends <- c(starts[-1L] - 1L, length(lines))
  field_value <- function(block, field) {
    prefix <- paste0("!Sample_", field, " = ")
    value <- block[startsWith(block, prefix)]
    if (length(value) == 0L) {
      return(NA_character_)
    }
    sub(prefix, "", value[[1L]], fixed = TRUE)
  }
  characteristic_value <- function(block, field) {
    prefix <- "!Sample_characteristics_ch1 = "
    values <- sub(
      prefix,
      "",
      block[startsWith(block, prefix)],
      fixed = TRUE
    )
    pattern <- paste0("^", field, ":\\s*")
    selected <- values[grepl(pattern, values, ignore.case = TRUE)]
    if (length(selected) == 0L) {
      return(NA_character_)
    }
    sub(pattern, "", selected[[1L]], ignore.case = TRUE)
  }
  records <- lapply(seq_along(starts), function(index) {
    block <- lines[starts[[index]]:ends[[index]]]
    data.frame(
      Source_GEO_Record_ID = sub(
        "^\\^SAMPLE = ",
        "",
        block[[1L]]
      ),
      Source_Title = field_value(block, "title"),
      Source_Tissue = characteristic_value(block, "tissue"),
      Source_Cell_Type = characteristic_value(block, "cell type"),
      Source_Age = characteristic_value(block, "age"),
      Source_Capture_Device_ID = characteristic_value(
        block,
        "c1 chip id"
      ),
      Source_Experiment_Sample_Name = characteristic_value(
        block,
        "experiment_sample_name"
      ),
      Source_Platform_ID = field_value(block, "platform_id"),
      Source_Instrument_Model = field_value(block, "instrument_model"),
      Source_Supplementary_File = field_value(
        block,
        "supplementary_file"
      ),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, records)
}

source_metadata <- parse_geo_soft_samples(
  file.path(dirname(data_dir), "GSE67835_gsm.soft.txt")
)
required_source_fields <- c(
  "Source_GEO_Record_ID", "Source_Tissue", "Source_Cell_Type",
  "Source_Age", "Source_Capture_Device_ID",
  "Source_Experiment_Sample_Name", "Source_Platform_ID",
  "Source_Instrument_Model"
)
if (
  nrow(source_metadata) != 466L ||
    anyDuplicated(source_metadata$Source_GEO_Record_ID) ||
    anyNA(source_metadata[, required_source_fields, drop = FALSE])
) {
  stop("GSE67835 GEO SOFT sample metadata is incomplete")
}

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

thisutils::log_message("Reading CSV files...")
cell_data_list <- thisutils::parallelize_fun(
  csv_files,
  fun = read_cell_data,
  cores = 10
)

all_genes <- unique(
  unlist(lapply(cell_data_list, function(x) names(x$expression)))
)
thisutils::log_message("Found {.val {length(all_genes)}} unique genes")

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
if (
  anyDuplicated(cell_names) ||
    !setequal(cell_names, source_metadata$Source_GEO_Record_ID)
) {
  stop("GSE67835 count files do not match all 466 GEO sample records")
}
thisutils::log_message("Creating Seurat object...")
object <- CreateSeuratObject(
  counts = expr_matrix,
  project = "GSE67835"
)

thisutils::log_message("Save data...")
saveRDS(
  object,
  file.path(res_dir, "GSE67835.rds")
)

source_index <- match(cell_names, source_metadata$Source_GEO_Record_ID)
metadata <- source_metadata[source_index, , drop = FALSE]
rownames(metadata) <- cell_names
metadata$Cells <- rownames(metadata)
metadata$Dataset <- "GSE67835"
metadata$Technology <- "Fluidigm C1"
metadata$Sequence <- "scRNA-seq"
metadata$Sample <- metadata$Source_Experiment_Sample_Name
metadata$Sample_ID <- metadata$Source_GEO_Record_ID
metadata$Original_Donor_ID <- metadata$Source_Experiment_Sample_Name
metadata$Original_Sample_ID <- metadata$Source_Experiment_Sample_Name
metadata$Original_Specimen_ID <- paste(
  metadata$Source_Experiment_Sample_Name,
  metadata$Source_Tissue,
  sep = "::"
)
metadata$Original_Library_ID <- metadata$Source_GEO_Record_ID
metadata$Original_Source_Sample_ID <-
  metadata$Source_Experiment_Sample_Name
metadata$Original_Source_Record_ID <- metadata$Source_GEO_Record_ID
metadata$Original_Technical_Batch_ID <-
  metadata$Source_Instrument_Model
metadata$Source_Single_Cell_Library_ID <-
  metadata$Source_GEO_Record_ID
metadata$Source_Full_Series_Cell_Count <- 466L
metadata$Source_Retained_Cell_Count <- 466L
metadata$CellType_raw <- metadata$Source_Cell_Type
metadata$source_cell_type_original_label <- metadata$Source_Cell_Type
metadata$source_cell_type_label <- metadata$Source_Cell_Type
metadata$source_cell_type_level_1 <- metadata$Source_Cell_Type
metadata$source_cell_type_original_label_source_column <-
  "Source_Cell_Type"
metadata$source_cell_type_label_source_column <- "Source_Cell_Type"
metadata$source_cell_type_level_1_source_column <- "Source_Cell_Type"
metadata$Brain_Region <- dplyr::recode(
  metadata$Source_Tissue,
  cortex = "Cerebral cortex",
  hippocampus = "Hippocampus",
  .default = NA_character_
)
if (anyNA(metadata$Brain_Region)) {
  stop("GSE67835 contains an unsupported source tissue label")
}
metadata$Region <- metadata$Brain_Region

prenatal <- startsWith(tolower(metadata$Source_Age), "prenatal")
adult_age <- suppressWarnings(as.numeric(sub(
  "^postnatal ([0-9]+) years$",
  "\\1",
  metadata$Source_Age,
  perl = TRUE
)))
if (any(!prenatal & is.na(adult_age))) {
  stop("GSE67835 contains an unsupported postnatal age label")
}
metadata$Age_Source_Raw <- metadata$Source_Age
metadata$Age_Source_Unit <- ifelse(
  prenatal,
  "gestational weeks",
  "years"
)
metadata$Age_Source_Basis <- ifelse(
  prenatal,
  "gestational age",
  "postnatal age"
)
metadata$Age_Source_Reference <- paste(
  "Darmanis et al. 2015, DOI 10.1073/pnas.1507125112, reports",
  "16-18 gestational weeks for prenatal specimens; GEO GSE67835",
  "provides the cell-level source age strings"
)
metadata$Age_Harmonization_Input <- ifelse(
  prenatal,
  "14-16 PCW",
  paste(adult_age, "years")
)
metadata$Age_Conversion_Formula <- ifelse(
  prenatal,
  "source 16-18 gestational weeks - 2 weeks = 14-16 PCW",
  "reported postnatal years retained"
)
metadata$Age_Conversion_Confidence <-
  "high: publication-defined age basis and cell-level GEO age"
metadata$Age_Conversion_Applied <- prenatal
metadata$Age <- metadata$Age_Harmonization_Input
metadata$Age_Interval_Curation_Evidence <- ifelse(
  prenatal,
  paste(
    "study-level 14-16 PCW range represented by its 15 PCW midpoint and",
    "assigned to S4 for descriptive interval coverage"
  ),
  NA_character_
)
metadata$Age_Interval_Curation_Confidence <- ifelse(
  prenatal,
  "moderate: source range is published; exact prenatal donor ages are unavailable",
  NA_character_
)
metadata$Continuous_Age_Eligibility <- ifelse(
  prenatal,
  "not eligible: exact donor age is not available",
  "eligible when the source-reported exact postnatal age is retained"
)
metadata$Sex <- NA_character_
metadata$Sex_Source_Raw <- NA_character_
metadata$Sex_Source_Column <-
  "not reported in the public GSE67835 sample records"
metadata$Sex_Assignment_Method <- "not reported; not inferred"
metadata$Sequencing_Platform <- metadata$Source_Instrument_Model
metadata$Diagnosis_raw <- ifelse(
  prenatal,
  "prenatal donor; diagnosis not reported in GEO",
  paste(
    "medically refractory seizures/mesial temporal sclerosis;",
    "resected tissue considered normal by EEG and pathology"
  )
)


column_order <- c(
  "Cells", "Dataset", "Technology", "Sequence", "Sample",
  "Sample_ID", "CellType_raw", "Brain_Region", "Region",
  "Age_Source_Raw", "Age_Source_Unit", "Age_Source_Basis",
  "Age_Source_Reference", "Age_Harmonization_Input",
  "Age_Conversion_Formula", "Age_Conversion_Confidence",
  "Age_Conversion_Applied", "Age", "Sex"
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
  file.path(res_dir, "GSE67835_processed.rds")
)
