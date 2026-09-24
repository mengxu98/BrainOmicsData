source("functions/prepare_env.R")
source("functions/dataset_metadata.R")

data_dir <- brainomics_data_path("raw/GSE97942/")
res_dir <- check_dir(brainomics_data_path("processed/GSE97942/"))

thisutils::log_message("Start loading data...")
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

source_library_map <- data.frame(
  Source_Library_Title = c(
    "occ1", "occ4", "occ5", "occ6", "occ7", "occ8",
    "occ9", "occ10", "occ11", "occ12", "occ2", "occ3",
    "fcx1", "fcx2", "occ13", "occ14", "occ15", "occ16",
    "fcx3", "fcx4", "fcx5", "occ17", "fcx6", "fcx7",
    "fcx8", "fcx9", "occ18", "occ19", "fcx10", "fcx11",
    "occ20", "occ21", "cbm1", "cbm2", "cbm3", "cbm4",
    "occ22", "occ23", "occ24", "cbm5", "cbm6", "cbm7",
    "fcx12", "fcx13", "cbm8", "cbm9"
  ),
  Source_GEO_Record_ID = c(
    "GSM2581267", "GSM2581268", "GSM2581269", "GSM2581270",
    "GSM2581271", "GSM2581272", "GSM2581273", "GSM2581274",
    "GSM2581275", "GSM2581276", "GSM2581277", "GSM2581278",
    "GSM2581279", "GSM2581280", "GSM2581281", "GSM2581282",
    "GSM2581283", "GSM2581284", "GSM2581285", "GSM2581286",
    "GSM2581287", "GSM2581288", "GSM2581289", "GSM2581290",
    "GSM2581291", "GSM2581292", "GSM2581293", "GSM2581294",
    "GSM2581295", "GSM2581296", "GSM2581297", "GSM2581298",
    "GSM2734253", "GSM2734254", "GSM2734255", "GSM2734256",
    "GSM2734257", "GSM2734258", "GSM2734259", "GSM2734260",
    "GSM2734261", "GSM2734262", "GSM2734263", "GSM2734264",
    "GSM2734265", "GSM2734266"
  ),
  Source_Experiment_ID = c(
    "20161027", rep("20161130B", 6L),
    rep("20161207A", 2L), rep("20161207B", 3L),
    rep("20170214A", 2L), rep("20170215A", 2L),
    rep("20170308A", 2L), rep("20170308B", 3L),
    "20170310A", rep("20170310B", 2L),
    rep("20170314A", 2L), rep("20170314B", 2L),
    rep("20170315A", 2L), rep("20170315B", 2L),
    rep("20170321A", 2L), rep("20170321B", 2L),
    rep("20170321C", 3L), rep("20170322A", 3L),
    rep("20170322B", 2L), rep("20170322C", 2L)
  ),
  Source_Donor_ID = c(
    rep("5342", 14L), rep("4590", 4L), rep("5342", 3L),
    "4590", rep("5342", 2L), rep("0006-YO", 4L),
    rep("0007-OX", 4L), rep("0006-YO", 2L),
    rep("0007-OX", 2L), rep("GEO_anonymous_35F", 8L),
    rep("GEO_anonymous_49M", 2L)
  ),
  stringsAsFactors = FALSE
)
if (nrow(source_library_map) != 46L ||
  anyDuplicated(source_library_map$Source_Library_Title) ||
  anyDuplicated(source_library_map$Source_GEO_Record_ID) ||
  length(unique(source_library_map$Source_Experiment_ID)) != 20L ||
  length(unique(source_library_map$Source_Donor_ID)) != 6L) {
  stop("GSE97942 GEO library map does not match the published design")
}
source_library_title <- tolower(sub("^D7_", "", metadata$orig.ident))
source_index <- match(
  source_library_title,
  source_library_map$Source_Library_Title
)
if (anyNA(source_index) ||
  length(unique(source_library_title)) != 46L) {
  stop("GSE97942 retained cells do not cover the 46 GEO libraries")
}
metadata$Source_Library_Title <- source_library_title
metadata$Source_GEO_Record_ID <-
  source_library_map$Source_GEO_Record_ID[source_index]
metadata$Source_Experiment_ID <-
  source_library_map$Source_Experiment_ID[source_index]
metadata$Source_Donor_ID <-
  source_library_map$Source_Donor_ID[source_index]
metadata$Source_Donor_ID_Verification_Status <- ifelse(
  grepl("^GEO_anonymous_", metadata$Source_Donor_ID),
  paste(
    "curated anonymous donor group from GEO age, sex, region and the",
    "published six-donor design"
  ),
  "source GEO patient identifier"
)

metadata$Cells <- rownames(metadata)
metadata$Dataset <- "GSE97942"
metadata$Technology <- metadata$protocal
metadata$Sequence <- "snRNA-seq"
metadata$Sample <- metadata$Source_Donor_ID
metadata$Sample_ID <- metadata$Source_GEO_Record_ID
metadata$CellType_raw <- metadata$priCluster
metadata$Brain_Region <- metadata$Area
metadata$Region <- metadata$Area
metadata <- standardize_source_age_metadata(
  metadata,
  source_age = as.character(metadata$Age),
  source_unit = "years",
  source_reference = "Lake et al. 2018, DOI 10.1038/nbt.4038"
)
metadata$Sex_Source_Raw <- as.character(metadata$Sex)
metadata$Sex_Assignment_Method <- "source reported"
metadata$Sex <- standardize_source_sex_value(metadata$Sex_Source_Raw)

column_order <- c(
  "Cells", "Dataset", "Technology", "Sequence", "Sample",
  "Sample_ID", "CellType_raw", "Brain_Region", "Region",
  "Age_Source_Raw", "Age_Source_Unit", "Age_Source_Basis",
  "Age_Source_Reference", "Age_Harmonization_Input",
  "Age_Conversion_Formula", "Age_Conversion_Confidence",
  "Age_Conversion_Applied", "Age", "Sex"
)
metadata <- retain_source_metadata(metadata, column_order)

counts <- counts[, metadata$Cells]
object <- CreateSeuratObject(
  counts = counts,
  meta.data = metadata
)

thisutils::log_message("Save data...")
saveRDS(
  object,
  file.path(res_dir, "GSE97942_processed.rds")
)
