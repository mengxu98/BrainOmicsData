#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(SeuratObject)
})

source("functions/data_paths.R")
source("functions/metadata_schema.R")
source("functions/sample_schema.R")
source("functions/dataset_metadata.R")
source("functions/processed_object.R")
source("functions/original_source_readers.R")
source("functions/integration.R")
source("functions/dataset_config.R")
source("processing/process_layered_h5ad.R")
source("sciencedb/package_metadata.R")

stopifnot(identical(brainomics_metadata_schema_version(), "1.6.1"))

canonical_sample_fixture <- data.frame(
  Dataset = c("Fixture", "Fixture"),
  Original_Sample_ID = c("source-library-1", "source-library-2"),
  Donor_ID = c("Fixture:D1", "Fixture:D2"),
  Sample_ID = c("Fixture:S1", "Fixture:S2"),
  Library_ID = c("Fixture:L1", "Fixture:L2"),
  stringsAsFactors = FALSE
)
canonical_sample_result <- add_sample_schema(canonical_sample_fixture)
stopifnot(
  identical(
    canonical_sample_result$Original_Sample_ID,
    canonical_sample_fixture$Original_Sample_ID
  ),
  identical(
    canonical_sample_result$Donor_ID,
    canonical_sample_fixture$Donor_ID
  ),
  identical(
    canonical_sample_result$Sample_ID,
    canonical_sample_fixture$Sample_ID
  ),
  identical(
    canonical_sample_result$Library_ID,
    canonical_sample_fixture$Library_ID
  ),
  all(canonical_sample_result$sample_schema_rule ==
    "canonical_identifiers_preserved")
)

repo_dir <- normalizePath(".", mustWork = TRUE)
test_dir <- tempfile("brainomics-clean-room-")
dir.create(test_dir)
on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)
# Source-recovery tests must use synthetic inputs, never the user's raw data.
brainomics_data_root <- function() file.path(test_dir, "data")

soft_fixture <- file.path(test_dir, "fixture_family.soft.gz")
soft_lines <- c(
  "^SAMPLE = GSM000001",
  "!Sample_title = 1: cortex",
  "!Sample_geo_accession = GSM000001",
  "!Sample_characteristics_ch1 = group: Control",
  "!Sample_characteristics_ch1 = Sex: Female",
  "^SAMPLE = GSM000002",
  "!Sample_title = 2: cortex",
  "!Sample_geo_accession = GSM000002",
  "!Sample_characteristics_ch1 = group: Case",
  "!Sample_characteristics_ch1 = Sex: Male"
)
soft_connection <- gzfile(soft_fixture, "wt")
writeLines(soft_lines, soft_connection)
close(soft_connection)
soft_metadata <- read_geo_family_soft_metadata(soft_fixture)
stopifnot(
  identical(soft_metadata$Sample_geo_accession, c("GSM000001", "GSM000002")),
  identical(soft_metadata$characteristic_group, c("Control", "Case")),
  identical(soft_metadata$characteristic_sex, c("Female", "Male"))
)

source_features <- c("gene-a", "gene-b", "gene-c")
source_cells <- paste0("source-cell-", 1:4)
source_layers <- list(
  donor_a = Matrix::sparseMatrix(
    i = c(1L, 2L),
    j = c(1L, 2L),
    x = c(1, 2),
    dims = c(3L, 2L),
    dimnames = list(source_features, source_cells[1:2])
  ),
  donor_b = Matrix::sparseMatrix(
    i = c(2L, 3L),
    j = c(1L, 2L),
    x = c(3, 4),
    dims = c(3L, 2L),
    dimnames = list(source_features, source_cells[3:4])
  )
)
source_matrix_fixture <- do.call(cbind, source_layers)
content_verification <- verify_materialized_matrix_content(
  source_matrix_fixture,
  source_matrix_fixture,
  block_cells = 2L
)
stopifnot(
  isTRUE(content_verification$verified),
  content_verification$blocks == 2L
)
altered_matrix_fixture <- source_matrix_fixture
altered_matrix_fixture@x[[1L]] <- altered_matrix_fixture@x[[1L]] + 1
stopifnot(inherits(
  try(
    verify_materialized_matrix_content(
      source_matrix_fixture,
      altered_matrix_fixture,
      block_cells = 2L
    ),
    silent = TRUE
  ),
  "try-error"
))
source_metadata <- data.frame(
  Cells = source_cells,
  Dataset = "fixture",
  row.names = source_cells,
  stringsAsFactors = FALSE
)
source_object <- create_complete_source_object(
  source_layers,
  source_metadata,
  "fixture"
)
stopifnot(
  ncol(source_object) == 4L,
  nrow(source_object) == 3L,
  identical(sort(colnames(source_object)), sort(source_cells)),
  identical(processed_object_matrix_class(source_object), "dgCMatrix_layers")
)

stopifnot(
  identical(
    normalize_gse202210_donor_key(c(
      "HSDG07HC_HHT", "hsDG101HC_HHT", "HSDG10HC"
    )),
    c("HSDG07HC", "HSDG101HC", "HSDG10HC")
  )
)

new_bundle_dir <- file.path(test_dir, "new-bundle")
dir.create(new_bundle_dir)
new_bundle_file <- file.path(new_bundle_dir, "artifact.tsv")
writeLines("fixture", new_bundle_file)
stopifnot(refresh_processed_bundle_sha256(
  new_bundle_dir,
  basename(new_bundle_file)
))
new_bundle_manifest <- file.path(
  new_bundle_dir,
  ".processed_bundle.sha256"
)
stopifnot(
  file.exists(new_bundle_manifest),
  identical(
    readLines(new_bundle_manifest, warn = FALSE),
    paste(
      processed_file_sha256(new_bundle_file),
      basename(new_bundle_file),
      sep = "  "
    )
  )
)

wang_config <- dataset_config("Wang_2025")
velmeshev_config <- dataset_config("Velmeshev_2023")
gse294786_config <- dataset_config("GSE294786")
stopifnot(
  identical(wang_config$field_map$library_id, "Sample_ID"),
  identical(wang_config$field_map$technical_batch_id, "Sample_ID"),
  identical(
    wang_config$field_map$cell_type_original,
    "Source_Cell_Type_Type_Updated"
  ),
  identical(wang_config$field_map$brain_region_source, "tissue"),
  grepl(
    "cell_type_ontology_term_id",
    wang_config$constants$cell_type_source_column,
    fixed = TRUE
  ),
  grepl("^record only:", wang_config$constants$technical_batch_use_status),
  identical(velmeshev_config$field_map$library_id, "sample"),
  identical(velmeshev_config$field_map$technical_batch_id, "sample"),
  identical(
    velmeshev_config$field_map$cell_type_original,
    "Source_Original_Lineage"
  ),
  identical(
    velmeshev_config$field_map$cell_type_cluster_id,
    "Source_Original_Seurat_Cluster"
  ),
  grepl(
    "Seurat_clusters",
    velmeshev_config$constants$cell_type_source_column,
    fixed = TRUE
  ),
  grepl(
    "^record only:",
    velmeshev_config$constants$technical_batch_use_status
  ),
  identical(
    gse294786_config$constants$brain_region,
    "Dorsolateral prefrontal cortex"
  ),
  identical(
    gse294786_config$constants$brain_region_ontology_id,
    "UBERON:0009834"
  ),
  identical(
    gse294786_config$constants$sequencing_platform,
    "Illumina HiSeq 2500"
  ),
  identical(
    gse294786_config$field_map$sex,
    "Source_Paper_Table_S1_Sex"
  )
)

gse294786_source_fixture <- data.frame(
  SampleID = c("1452", "1490"),
  Species = c("Hg", "Hg"),
  Age2 = c("150", "120"),
  CellType = c("OPC", "L2/L3"),
  CellType2 = c("OPC", "L2/L3"),
  cell_type = c("OPC", "L2/L3"),
  cell_type_final = c("OPC", "L2/L3"),
  cell_type_final_broad = c("OPC", NA_character_),
  integrated_snn_res.0.8 = c("1", "2"),
  seurat_clusters = c("1", "2"),
  leiden = c("0", "1"),
  leiden_0_5 = c("0", NA_character_),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
gse294786_transformed_fixture <-
  gse294786_config$raw_metadata_transform(gse294786_source_fixture)
stopifnot(
  identical(
    gse294786_transformed_fixture$Age_Harmonization_Input,
    c("2 months", "4 months")
  ),
  identical(
    gse294786_transformed_fixture$Source_Age_Metadata_Paper_Agreement,
    c("conflict", "agree")
  ),
  identical(
    gse294786_transformed_fixture$Source_Paper_Table_S1_Sex,
    c("Male", "Male")
  )
)

gse294786_filter_fixture <- data.frame(
  Species = c("Homo sapiens", "Macaca mulatta", NA_character_),
  Analysis_Include = rep(TRUE, 3L),
  Analysis_Role = rep("reference", 3L),
  Exclusion_Reason = rep(NA_character_, 3L),
  stringsAsFactors = FALSE
)
gse294786_filtered_fixture <- apply_required_filters(
  gse294786_filter_fixture,
  gse294786_config
)
stopifnot(
  identical(
    gse294786_filtered_fixture$Analysis_Include,
    c(TRUE, FALSE, FALSE)
  ),
  identical(
    gse294786_filtered_fixture$Analysis_Role,
    c("reference", "reference_excluded", "reference_excluded")
  ),
  grepl(
    "non-human species excluded",
    gse294786_filtered_fixture$Exclusion_Reason[[2L]],
    fixed = TRUE
  ),
  grepl(
    "no matching released cell metadata",
    gse294786_filtered_fixture$Exclusion_Reason[[3L]],
    fixed = TRUE
  )
)

geo_fixture <- file.path(test_dir, "GSE_fixture_series_matrix.txt.gz")
geo_lines <- c(
  '!Sample_title\t"cell-1"\t"cell-2"',
  '!Sample_geo_accession\t"GSM1"\t"GSM2"',
  '!Sample_characteristics_ch1\t"donor: D1"\t"donor: D2"',
  '!Sample_characteristics_ch1\t"c1_chip: C1-1"\t"c1_chip: C1-2"',
  '!Sample_instrument_model\t"HiSeq 2000"\t"HiSeq 2000"',
  "!series_matrix_table_begin"
)
geo_connection <- gzfile(geo_fixture, "wt")
writeLines(geo_lines, geo_connection)
close(geo_connection)
geo_metadata <- read_geo_series_matrix_sample_metadata(geo_fixture)
stopifnot(
  identical(geo_metadata$Sample_geo_accession, c("GSM1", "GSM2")),
  identical(geo_metadata$characteristic_donor, c("D1", "D2")),
  identical(geo_metadata$characteristic_c1_chip, c("C1-1", "C1-2"))
)

identifier_fixture <- file.path(test_dir, "identifier_fixture.tsv.gz")
identifier_connection <- gzfile(identifier_fixture, "wt")
write.table(
  data.frame(
    Original_Donor_ID = c("0950", "1465"),
    orig.ident = c("0950", "1465"),
    case = c("0950_240109", "1465_dapi_ADpaper"),
    Age = c(53, 17),
    stringsAsFactors = FALSE
  ),
  identifier_connection,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
close(identifier_connection)
identifier_roundtrip <- read_brainomics_table(identifier_fixture)
stopifnot(
  identical(identifier_roundtrip$Original_Donor_ID, c("0950", "1465")),
  identical(identifier_roundtrip$orig.ident, c("0950", "1465")),
  identical(identifier_roundtrip$case, c("0950_240109", "1465_dapi_ADpaper")),
  is.numeric(identifier_roundtrip$Age)
)

# Age intervals use lower-inclusive and upper-exclusive boundaries. The final
# prenatal interval therefore contains 39 PCW and stops at 40 PCW.
age_values <- c(
  4, 8, 10, 13, 16, 19, 24, 39, 39.999, 40,
  0, 0.5, 1, 6, 12, 20, 40, 60
)
age_units <- c(rep("PCW", 10L), rep("years", 8L))
stopifnot(identical(
  dataset_age_interval(age_values, age_units),
  c(
    "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S7", "S7", NA,
    "S8", "S9", "S10", "S11", "S12", "S13", "S14", "S15"
  )
))
age_fixture <- add_age_schema(data.frame(
  Age = c("39 PCW", "40 PCW", "7-8 PCW", "12-20 years"),
  stringsAsFactors = FALSE
))
age_fixture <- add_metadata_schema(age_fixture)
open_age_fixture <- add_age_schema(data.frame(
  Age = "90+ years",
  Age_Source_Raw = "90+",
  stringsAsFactors = FALSE
))
stopifnot(
  age_fixture$Stage[[1L]] == "S7",
  age_fixture$AgeRange[[1L]] == "[24, 40) PCW",
  is.na(age_fixture$Stage[[2L]]),
  is.na(age_fixture$Stage[[3L]]),
  is.na(age_fixture$Stage[[4L]]),
  open_age_fixture$Stage[[1L]] == "S15",
  open_age_fixture$Age_Lower[[1L]] == 90,
  is.infinite(open_age_fixture$Age_Upper[[1L]]),
  open_age_fixture$Age_Representative_Method[[1L]] ==
    "reported open-ended lower bound",
  all(
    age_fixture$Age_Interval_Boundary_Convention ==
      "lower inclusive; upper exclusive"
  )
)

gestational_fixture <- data.frame(Age = c("8w", "30 weeks gestation"))
gestational_fixture <- standardize_source_age_metadata(
  gestational_fixture,
  source_age = gestational_fixture$Age,
  source_basis = "gestational",
  source_unit = "gestational weeks",
  source_reference = "synthetic gestational-age fixture"
)
postconception_fixture <- data.frame(Age = c("22w", "23w"))
postconception_fixture <- standardize_source_age_metadata(
  postconception_fixture,
  source_age = postconception_fixture$Age,
  source_basis = "postconceptional",
  source_unit = "postconceptional weeks",
  source_reference = "synthetic postconception-age fixture"
)
stopifnot(
  identical(gestational_fixture$Age, c("6 PCW", "28 PCW")),
  all(gestational_fixture$Age_Conversion_Applied),
  identical(postconception_fixture$Age, c("22 PCW", "23 PCW"))
)

region_fixture <- add_metadata_schema(data.frame(
  Dataset = c("Fixture_A", "Fixture_A", "Fixture_A", "HYPOMAP"),
  BrainRegion = c(
    "Prefrontal Cortex", "Cerebellar cortex",
    "Fornix/Optic tract/Anterior commissure",
    "Periventricular region"
  ),
  stringsAsFactors = FALSE
))
stopifnot(
  metadata_value_examples(
    data.frame(region = "Prefrontal cortex", stringsAsFactors = FALSE),
    "region"
  ) == "Prefrontal cortex"
)
stopifnot(
  identical(
    region_fixture$BrainRegion,
    c(
      "Prefrontal cortex", "Cerebellum",
      "Fornix/Optic tract/Anterior commissure",
      "Periventricular region"
    )
  ),
  region_fixture$brain_region_ontology_id[[1L]] == "UBERON:0000451",
  region_fixture$brain_region_ontology_id[[2L]] == "UBERON:0002037",
  is.na(region_fixture$brain_region_ontology_id[[3L]]),
  region_fixture$brain_region_ontology_mapping_status[[3L]] ==
    "composite region; no single ontology term",
  region_fixture$brain_region_ontology_id[[4L]] == "UBERON:0002271"
)

li_region_fixture <- add_metadata_schema(data.frame(
  Dataset = rep("Li_et_al_2018", 4L),
  BrainRegion = c("DFC", "NCX", "Pallium", "pallium"),
  stringsAsFactors = FALSE
))
stopifnot(
  identical(
    li_region_fixture$BrainRegion,
    c(
      "Dorsolateral prefrontal cortex", "Neocortex",
      "Pallium", "Pallium"
    )
  ),
  identical(
    li_region_fixture$brain_region_ontology_id,
    c(
      "UBERON:0009834", "UBERON:0001950",
      "UBERON:0000203", "UBERON:0000203"
    )
  ),
  all(li_region_fixture$brain_region_ontology_mapping_status ==
    "ontology mapped")
)

region_case_fixture <- add_metadata_schema(data.frame(
  Dataset = rep("Fixture_A", 9L),
  BrainRegion = c(
    "cerebral cortex", "frontal cortex", "ganglionic eminence",
    "neocortex", "prefrontal cortex", "primary motor cortex",
    "forebrain", "telencephalon", "temporal cortex"
  ),
  stringsAsFactors = FALSE
))
stopifnot(
  identical(
    region_case_fixture$BrainRegion,
    c(
      "Cerebral cortex", "Frontal cortex", "Ganglionic eminence",
      "Neocortex", "Prefrontal cortex", "Primary motor cortex",
      "Forebrain", "Telencephalon", "Temporal cortex"
    )
  ),
  identical(
    region_case_fixture$brain_region_raw,
    c(
      "cerebral cortex", "frontal cortex", "ganglionic eminence",
      "neocortex", "prefrontal cortex", "primary motor cortex",
      "forebrain", "telencephalon", "temporal cortex"
    )
  ),
  all(region_case_fixture$brain_region_ontology_mapping_status ==
    "ontology mapped")
)

expanded_region_fixture <- add_metadata_schema(data.frame(
  Dataset = rep("EGAD00001006049", 12L),
  BrainRegion = c(
    "Forebrain", "Striatum", "Brain", "Midbrain dorsal",
    "Cortex frontal", "Cortex temporal", "Cortex parietal",
    "Cortex occipital", "Caudate+Putamen", "Telencephalon",
    "Cortex hemisphere A", "Dorsolateral prefrontal cortex (BA9)"
  ),
  stringsAsFactors = FALSE
))
stopifnot(
  identical(
    expanded_region_fixture$brain_region_ontology_id,
    c(
      "UBERON:0001890", "UBERON:0002435", "UBERON:0000955",
      "UBERON:0002314", "UBERON:0001870", "UBERON:0016538",
      "UBERON:0016530", "UBERON:0016540", "UBERON:0005383",
      "UBERON:0001893", "UBERON:0000956", "UBERON:0009834"
    )
  ),
  all(expanded_region_fixture$brain_region_ontology_mapping_status ==
    "ontology mapped")
)

egad_region_cells <- c("egad-head", "egad-cortex")
egad_region_fixture <- canonicalize_existing_reference_metadata_frame(
  metadata = data.frame(
    Cells = egad_region_cells,
    Original_Donor_ID = "D1",
    Original_Sample_ID = c("S1", "S2"),
    Original_Specimen_ID = c("S1", "S2"),
    Original_Library_ID = c("L1", "L2"),
    Original_Source_Record_ID = c("R1", "R2"),
    Age = "8 PCW",
    Sex = "Unknown",
    BrainRegion = c("Head", "Cortex"),
    CellType_raw = c("Fibroblast", "Radial glia"),
    Analysis_Include = TRUE,
    row.names = egad_region_cells,
    stringsAsFactors = FALSE
  ),
  cells = egad_region_cells,
  dataset = "EGAD00001006049"
)
stopifnot(
  identical(egad_region_fixture$Analysis_Include, c(FALSE, TRUE)),
  egad_region_fixture$Analysis_Role[[1L]] == "reference_excluded",
  grepl(
    "cranial fibroblast",
    egad_region_fixture$Exclusion_Reason[[1L]],
    fixed = TRUE
  )
)

cells <- sprintf("cell-%03d", seq_len(120L))
features <- sprintf("gene-%02d", seq_len(40L))
counts_dense <- outer(
  seq_along(features),
  seq_along(cells),
  function(gene, cell) {
    ifelse((gene * 3L + cell * 5L) %% 11L < 2L, (gene + cell) %% 7L + 1L, 0L)
  }
)
dimnames(counts_dense) <- list(features, cells)
counts <- methods::as(Matrix::Matrix(counts_dense, sparse = TRUE), "dgCMatrix")

raw_meta <- data.frame(
  Original_Cell_ID = cells,
  donor = rep(sprintf("donor-%02d", 1:4), each = 30L),
  specimen = rep(sprintf("specimen-%02d", 1:8), each = 15L),
  library = rep(c("library-a", "library-b"), each = 60L),
  age = rep(c("8 PCW", "10 PCW", "6 years", "20 years"), each = 30L),
  sex = rep(c("Female", "Male", "Female", "Male"), each = 30L),
  region = rep(c("Prefrontal cortex", "Cerebellum"), each = 60L),
  source_cell_type = rep(c("Excitatory neuron", "Oligodendrocyte"), 60L),
  stringsAsFactors = FALSE
)
field_map <- list(
  cell_id = "Original_Cell_ID",
  donor_id = "donor",
  sample_id = "specimen",
  specimen_id = "specimen",
  library_id = "library",
  age = "age",
  sex = "sex",
  brain_region = "region",
  brain_region_source = "region",
  cell_type = "source_cell_type",
  cell_type_level_1 = "source_cell_type"
)
metadata <- build_dataset_metadata(
  raw_meta,
  dataset = "Fixture_A",
  field_map = field_map,
  constants = list(
    species = "Homo sapiens",
    technology = "10x Genomics",
    modality = "snRNA-seq",
    donor_id_semantics = "synthetic donor",
    specimen_id_semantics = "synthetic specimen",
    specimen_id_verification_status = "verified: synthetic fixture",
    library_id_semantics = "synthetic library",
    library_id_verification_status = "verified: synthetic fixture",
    library_id_evidence = "synthetic fixture library field",
    source_record_id_semantics = "not provided in synthetic fixture",
    technical_batch_semantics =
      "not provided in synthetic fixture",
    cell_type_source_column = "source_cell_type",
    cell_type_source_file = "synthetic fixture"
  ),
  default_analysis_role = "reference"
)
metadata <- add_source_publication_metadata(
  metadata,
  data.frame(
    source_accession = "SYNTHETIC",
    source_repository = "synthetic repository",
    publication_doi = "10.0000/synthetic.fixture",
    verified_title = "Synthetic metadata fixture",
    journal = "Synthetic Journal",
    publication_year = "2026",
    repository_record_url = "https://example.org/synthetic",
    stringsAsFactors = FALSE
  )
)
metadata <- add_age_schema(metadata)
metadata <- add_metadata_schema(metadata)
validate_dataset_metadata(metadata, matrix_cells = cells, expected_cells = 120L)
stopifnot(
  identical(unique(metadata$Sex), c("Female", "Male")),
  all(metadata$sex_provenance == "source reported"),
  !any(metadata$Sex_Donor_Conflict, na.rm = TRUE),
  all(
    metadata$source_cell_type_original_label_source_column ==
      "source_cell_type"
  ),
  all(
    metadata$source_cell_type_label_source_column ==
      "source_cell_type"
  ),
  all(
    metadata$source_cell_type_level_1_source_column ==
      "source_cell_type"
  ),
  all(metadata$Source_Publication_DOI == "10.0000/synthetic.fixture")
)
sex_fixture <- add_metadata_schema(data.frame(
  Global_Donor_ID = c("donor-a", "donor-a", "donor-b", NA),
  Specimen_ID = c("specimen-a", "specimen-b", "specimen-c", NA),
  Sex_Source_Raw = c("F", "male", "unknown", "non-binary"),
  Sex_Source_Column = "source_sex",
  Sex_Assignment_Method = "source reported",
  stringsAsFactors = FALSE
))
stopifnot(
  identical(sex_fixture$Sex, c("Female", "Male", NA, NA)),
  all(sex_fixture$Sex_Donor_Conflict[1:2]),
  sex_fixture$sex_standardized[[3L]] == "Not reported",
  sex_fixture$sex_standardized[[4L]] == "Other/unspecified"
)

dataset_batch <- assign_integration_batch(metadata, model = "dataset")
stopifnot(
  identical(
    unique(dataset_batch$Integration_Batch_ID),
    "dataset:Fixture_A"
  ),
  all(dataset_batch$Integration_Batch_Source == "source dataset")
)
technical_fixture <- metadata
technical_fixture$Original_Technical_Batch_ID <- rep(
  c("run-1", "run-2"),
  each = 60L
)
technical_fixture$Technical_Batch_ID <- paste0(
  "Fixture_A:technical_batch:",
  technical_fixture$Original_Technical_Batch_ID
)
technical_fixture$Technical_Batch_Verification_Status <-
  "verified: synthetic fixture"
technical_fixture$Technical_Batch_Evidence <-
  "each synthetic run spans two donors"
technical_fixture$Technical_Batch_Source_Column <- "synthetic_run"
technical_fixture$Technical_Batch_Derivation_Rule <-
  "direct source metadata field synthetic_run"
technical_fixture$Technical_Batch_Use_Status <-
  "eligible for technical-batch sensitivity model"
technical_fixture$Technical_Batch_Use_Evidence <-
  "synthetic runs are balanced across multiple donors"
technical_batch <- assign_integration_batch(
  technical_fixture,
  model = "verified_technical"
)
stopifnot(
  all(
    technical_batch$Integration_Batch_ID[1:60] ==
      "Fixture_A:technical_batch:run-1"
  ),
  all(
    technical_batch$Integration_Batch_ID[61:120] ==
      "Fixture_A:technical_batch:run-2"
  ),
  all(technical_batch$Technical_Batch_Dataset_Eligible),
  all(
    technical_batch$Integration_Batch_ID_Dataset == "dataset:Fixture_A"
  ),
  identical(
    technical_batch$Integration_Batch_ID,
    technical_batch$Integration_Batch_ID_Verified_Technical
  )
)
technical_audit <- integration_batch_audit(technical_batch)
technical_design_audit <- technical_batch_design_audit(technical_fixture)
technical_source_meta <- raw_meta
technical_source_meta$synthetic_run <- rep(
  c("run-1", "run-2"),
  each = 60L
)
technical_source_meta$lane <- rep(c("lane-1", "lane-2"), 60L)
technical_field_inventory <- source_technical_field_inventory(
  technical_fixture,
  source_meta = technical_source_meta
)
stopifnot(
  "Technical_Batch_Evidence" %in% names(technical_audit),
  any(
    technical_audit$Source_Technical_Batch_ID ==
      "Fixture_A:technical_batch:run-1"
  ),
  technical_design_audit$Technical_Batch_Design_Gate_Passed[[1L]],
  technical_design_audit$Technical_Batch_Coverage_Fraction[[1L]] == 1,
  technical_design_audit$Technical_Batch_Levels[[1L]] == 2L,
  all(technical_design_audit$Spans_Multiple_Donors_Or_Specimens),
  any(
    technical_field_inventory$Source_Field == "synthetic_run" &
      technical_field_inventory$Exact_Audited_Technical_Batch_Source
  ),
  any(
    technical_field_inventory$Source_Field == "lane" &
      !technical_field_inventory$Exact_Audited_Technical_Batch_Source
  )
)
analysis_confounded_technical <- technical_fixture
analysis_confounded_technical$Analysis_Include <- rep(FALSE, 120L)
analysis_confounded_technical$Analysis_Include[c(1:15, 61:75)] <- TRUE
analysis_confounded_technical$Exclusion_Reason[
  !analysis_confounded_technical$Analysis_Include
] <- "synthetic species filter"
analysis_confounded_audit <- technical_batch_design_audit(
  analysis_confounded_technical
)
stopifnot(
  analysis_confounded_audit$All_Batch_Levels_Span_Multiple_Donors_Or_Specimens[[1L]],
  !analysis_confounded_audit$All_Analysis_Batch_Levels_Span_Multiple_Donors_Or_Specimens[[1L]],
  !analysis_confounded_audit$Technical_Batch_Design_Gate_Passed[[1L]],
  !analysis_confounded_audit$Declared_Eligibility_Consistent_With_Design_Gate[[1L]],
  analysis_confounded_audit$Analysis_Technical_Batch_Levels[[1L]] == 2L,
  analysis_confounded_audit$Analysis_Technical_Batch_Coverage_Fraction[[1L]] == 1
)
incomplete_technical <- technical_fixture
incomplete_technical$Original_Technical_Batch_ID[120] <- NA_character_
incomplete_technical$Technical_Batch_ID[120] <- NA_character_
incomplete_technical$Technical_Batch_Verification_Status[120] <-
  "not available"
incomplete_technical$Technical_Batch_Evidence[120] <- NA_character_
incomplete_technical$Technical_Batch_Use_Status[120] <-
  "not available; record only"
incomplete_technical$Technical_Batch_Use_Evidence[120] <-
  "synthetic batch identifier intentionally missing"
incomplete_design_audit <- technical_batch_design_audit(
  incomplete_technical
)
incomplete_design_summary <- incomplete_design_audit[
  incomplete_design_audit$Record_Type == "dataset_summary", ,
  drop = FALSE
]
incomplete_design_levels <- incomplete_design_audit[
  incomplete_design_audit$Record_Type == "batch_level", ,
  drop = FALSE
]
incomplete_batch <- assign_integration_batch(
  incomplete_technical,
  model = "verified_technical"
)
stopifnot(
  incomplete_design_summary$Technical_Batch_Known_Cells[[1L]] == 119L,
  incomplete_design_summary$Technical_Batch_Missing_Cells[[1L]] == 1L,
  sum(incomplete_design_levels$Cells) == 119L,
  !anyNA(incomplete_design_levels$Cells),
  all(incomplete_batch$Integration_Batch_ID == "dataset:Fixture_A"),
  !any(incomplete_batch$Technical_Batch_Dataset_Eligible)
)
stopifnot(identical(
  metadata_count_values(c("0", "119"), "fixture counts"),
  c(0, 119)
))
invalid_count_error <- tryCatch(
  {
    metadata_count_values(c("1", "unknown"), "fixture counts")
    FALSE
  },
  error = function(error) {
    grepl(
      "missing or non-integer text values",
      conditionMessage(error),
      fixed = TRUE
    )
  }
)
stopifnot(invalid_count_error)

ineligible_technical <- technical_fixture
ineligible_technical$Technical_Batch_Use_Status <-
  "record only: ineligible for technical-batch correction"
ineligible_technical$Technical_Batch_Use_Evidence <-
  "synthetic technical factor is biologically confounded"
ineligible_batch <- assign_integration_batch(
  ineligible_technical,
  model = "verified_technical"
)
stopifnot(
  all(ineligible_batch$Integration_Batch_ID == "dataset:Fixture_A"),
  !any(ineligible_batch$Technical_Batch_Dataset_Eligible)
)

# Source annotation fields are coalesced cell by cell without overwriting the
# original columns retained in the preprocessing object.
annotation_fixture <- data.frame(
  primary = c("Astrocyte", NA, "Microglia"),
  backup = c("Astrocyte broad", "Oligodendrocyte", NA),
  stringsAsFactors = FALSE
)
annotation_value <- first_reported_metadata_field(
  annotation_fixture,
  c("primary", "backup")
)
stopifnot(
  identical(
    annotation_value$value,
    c("Astrocyte", "Oligodendrocyte", "Microglia")
  ),
  identical(annotation_value$name, c("primary", "backup", "primary"))
)

gse168408_annotation_fixture <- add_existing_source_cell_type_metadata(
  data.frame(
    cell_type = c("PN", "Non-Neu"),
    major_clust = c("PN_dev", "Astro"),
    sub_clust = c("PN_dev", "Astro_dev-1"),
    leiden = c("13", "2"),
    stringsAsFactors = FALSE
  ),
  "GSE168408"
)
stopifnot(
  identical(
    gse168408_annotation_fixture$source_cell_type_original_label,
    c("PN_dev", "Astro_dev-1")
  ),
  identical(
    gse168408_annotation_fixture$source_cell_type_label,
    c("PN_dev", "Astro_dev-1")
  ),
  all(
    gse168408_annotation_fixture$
      source_cell_type_original_label_source_column == "sub_clust"
  ),
  all(grepl(
    "original-author annotation",
    gse168408_annotation_fixture$source_cell_type_label_semantics
  ))
)

# GSE217511 is a dataset-specific exception: the paper states that every
# donor-region sample was sequenced as a separate library, and the deposited
# GSM matrices map one-to-one to those samples. ROSMAP Batch is independently
# admitted only with explicit paper/source evidence.
gse103723_donor_key <- rep(c("22WF", "23WF", "23WM"), each = 17L)
gse103723_cells <- paste0("gse103723-", seq_len(51L))
gse103723_fixture <- data.frame(
  Cells = gse103723_cells,
  Sample = paste0(
    "R",
    rep(seq_len(17L), 3L),
    "_",
    gse103723_donor_key,
    "_B",
    seq_len(51L)
  ),
  Sample_ID = paste0("GSM", seq_len(51L)),
  batch = as.character(seq_len(51L) - 1L),
  cell_type = rep(c("Neuron", "Glia", "Progenitor"), 17L),
  Age = ifelse(gse103723_donor_key == "22WF", "22 PCW", "23 PCW"),
  Brain_Region = "Cerebral cortex",
  stringsAsFactors = FALSE,
  row.names = gse103723_cells
)
gse103723_meta <- upgrade_existing_reference_metadata(
  gse103723_fixture,
  "GSE103723"
)
stopifnot(
  identical(
    sort(unique(gse103723_meta$Original_Donor_ID)),
    c("22WF", "23WF", "23WM")
  ),
  length(unique(gse103723_meta$Original_Library_ID)) == 51L,
  all(gse103723_meta$Library_ID_Source_Column == "Sample_ID"),
  all(grepl(
    "^record only:",
    gse103723_meta$Technical_Batch_Use_Status
  ))
)

gse186538_cells <- paste0("gse186538-", seq_len(25L))
gse186538_donor <- c(
  rep("HSB179", 7L), rep("HSB181", 7L), rep("HSB282", 8L),
  "HSB231", "HSB237", "HSB628"
)
gse186538_fixture <- data.frame(
  Cells = gse186538_cells,
  Sample = gse186538_donor,
  Sample_ID = "GSE186538",
  batch = paste0(
    gse186538_donor,
    "_",
    ave(
      seq_along(gse186538_donor),
      gse186538_donor,
      FUN = seq_along
    )
  ),
  original_name = "Astrocyte",
  Age = "50 years",
  Brain_Region = "Hippocampus",
  Region = "Dentate gyrus",
  stringsAsFactors = FALSE,
  row.names = gse186538_cells
)
gse186538_meta <- upgrade_existing_reference_metadata(
  gse186538_fixture,
  "GSE186538"
)
stopifnot(
  length(unique(gse186538_meta$Original_Library_ID)) == 25L,
  all(gse186538_meta$Library_ID_Source_Column == "batch"),
  all(grepl(
    "^record only:",
    gse186538_meta$Technical_Batch_Use_Status
  ))
)

gse217511_fixture <- data.frame(
  Cells = c("cell-a", "cell-b"),
  Sample = c("2C", "2G"),
  Sample_ID = c("GSM1", "GSM2"),
  celltypes = c("ExN", "ExN"),
  Age = c("30 PCW", "30 PCW"),
  BrainRegion = c("Cortex", "Ganglionic eminence"),
  stringsAsFactors = FALSE,
  row.names = c("cell-a", "cell-b")
)
gse217511_meta <- upgrade_existing_reference_metadata(
  gse217511_fixture,
  "GSE217511"
)
stopifnot(
  identical(
    gse217511_meta$Library_ID,
    c("GSE217511:library:GSM1", "GSE217511:library:GSM2")
  ),
  all(gse217511_meta$Library_ID_Source_Column == "Sample_ID"),
  identical(
    gse217511_meta$Source_Record_ID,
    c("GSE217511:source_record:GSM1", "GSE217511:source_record:GSM2")
  ),
  all(grepl(
    "^record only:",
    gse217511_meta$Technical_Batch_Use_Status
  )),
  all(gse217511_meta$source_cell_type_label == "ExN")
)

gse204683_source_dir <- brainomics_data_path("raw", "GSE204683")
dir.create(gse204683_source_dir, recursive = TRUE)
write.table(data.frame(
  "Donor ID" = c("Child2", "Child1"),
  Barcode = c("bc-b", "bc-a"),
  "Cell type" = c("ExN", "ExN"),
  check.names = FALSE
), file.path(gse204683_source_dir, "GSE204683_barcodes.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)
gse204683_fixture <- data.frame(
  Cells = c("6032_bc-a", "5977_bc-b"),
  Sample = c("Child1", "Child2"),
  Sample_ID = c("6032", "5977"),
  CellType_raw = c("ExN", "ExN"),
  Age = c("4", "6"),
  Brain_Region = rep("Prefrontal cortex", 2L),
  stringsAsFactors = FALSE,
  row.names = c("6032_bc-a", "5977_bc-b")
)
gse204683_meta <- upgrade_existing_reference_metadata(
  gse204683_fixture,
  "GSE204683"
)
stopifnot(
  identical(gse204683_meta$Source_CellType, c("ExN", "ExN")),
  identical(
    gse204683_meta$Library_ID,
    c("GSE204683:library:6032", "GSE204683:library:5977")
  ),
  all(gse204683_meta$Library_ID_Source_Column == "Sample_ID"),
  all(grepl("one library", gse204683_meta$Library_ID_Evidence)),
  all(grepl(
    "^record only:",
    gse204683_meta$Technical_Batch_Use_Status
  )),
  all(grepl(
    "CELLxGENE Batch field",
    gse204683_meta$Technical_Batch_Use_Evidence
  ))
)

gse168408_fixture <- data.frame(
  Cells = paste0("gse168408-", letters[1:6]),
  Sample = paste0("RL", 1:6),
  Sample_ID = paste0("RL", 1:6, "_age_chem"),
  batch = paste0("RL", 1:6, "_age_chem"),
  Herring_Technical_Batch_Source = rep(
    c("GPL21697_v2", "GPL21697_v3", "GPL24676_v3"),
    each = 2L
  ),
  CellType_raw = rep("PN_dev", 6L),
  Age = rep("22 PCW", 6L),
  BrainRegion = "Prefrontal cortex",
  stringsAsFactors = FALSE,
  row.names = paste0("gse168408-", letters[1:6])
)
gse168408_meta <- upgrade_existing_reference_metadata(
  gse168408_fixture,
  "GSE168408"
)
gse168408_batch <- assign_integration_batch(
  gse168408_meta,
  model = "verified_technical"
)
stopifnot(
  identical(
    gse168408_meta$Original_Library_ID,
    gse168408_fixture$batch
  ),
  all(gse168408_meta$Library_ID_Source_Column == "batch"),
  all(grepl(
    "^eligible for technical-batch",
    gse168408_meta$Technical_Batch_Use_Status
  )),
  setequal(
    unique(gse168408_batch$Integration_Batch_ID),
    paste0(
      "GSE168408:technical_batch:",
      c("GPL21697_v2", "GPL21697_v3", "GPL24676_v3")
    )
  )
)

gse207334_fixture <- data.frame(
  Cells = c("gse207334-a", "gse207334-b"),
  Sample = c("2RT00374N", "RT00382N"),
  Sample_ID = c("HSB6195", "HSB5871"),
  Source_GEO_RNA_Record_ID = c("GSM6284664", "GSM6284665"),
  subclass = c("Astrocyte", "Microglia"),
  Age = c("45", "60"),
  BrainRegion = "Dorsolateral prefrontal cortex",
  stringsAsFactors = FALSE,
  row.names = c("gse207334-a", "gse207334-b")
)
gse207334_meta <- upgrade_existing_reference_metadata(
  gse207334_fixture,
  "GSE207334"
)
stopifnot(
  identical(
    gse207334_meta$Original_Library_ID,
    gse207334_fixture$Source_GEO_RNA_Record_ID
  ),
  all(
    gse207334_meta$Library_ID_Source_Column ==
      "Source_GEO_RNA_Record_ID"
  ),
  all(grepl("^record only:", gse207334_meta$Technical_Batch_Use_Status))
)

ma_fixture <- data.frame(
  Cells = paste0("ma-", 1:4),
  Sample = rep(c("HSB106", "HSB189"), each = 2L),
  Sample_ID = rep(c("HSB106", "HSB189"), each = 2L),
  tech_rep = c("HSB106_1", "HSB106_2", "HSB189_1", "HSB189_2"),
  class = rep("Non-neuronal", 4L),
  subclass = c("Astrocyte", "Microglia", "Astrocyte", "Microglia"),
  subtype = c("Astro-1", "Micro-1", "Astro-1", "Micro-1"),
  Age = rep(c("64", "36"), each = 2L),
  BrainRegion = "Dorsolateral prefrontal cortex",
  stringsAsFactors = FALSE,
  row.names = paste0("ma-", 1:4)
)
ma_meta <- upgrade_existing_reference_metadata(ma_fixture, "Ma_et_al_2022")
ma_batch <- assign_integration_batch(ma_meta, model = "verified_technical")
stopifnot(
  identical(ma_meta$Original_Library_ID, ma_fixture$tech_rep),
  all(ma_meta$Library_ID_Source_Column == "tech_rep"),
  identical(ma_meta$source_cell_type_level_1, ma_fixture$class),
  identical(ma_meta$source_cell_type_level_2, ma_fixture$subclass),
  identical(ma_meta$source_cell_type_level_3, ma_fixture$subtype),
  all(grepl("^record only:", ma_meta$Technical_Batch_Use_Status)),
  all(ma_batch$Integration_Batch_ID == "dataset:Ma_et_al_2022")
)

gse261983_fixture <- data.frame(
  Cells = c("gse261983-a", "gse261983-b"),
  Sample = c("RT00372N", "RT00374N"),
  Sample_ID = c("GSM1_RT00372N", "GSM2_RT00374N"),
  Source_GEO_Record_ID = c("GSM1_RT00372N", "GSM2_RT00374N"),
  anno = c("Astrocyte", "Microglia"),
  Age = c("89", "73"),
  Brain_Region = rep("Prefrontal cortex", 2L),
  stringsAsFactors = FALSE,
  row.names = c("gse261983-a", "gse261983-b")
)
gse261983_meta <- upgrade_existing_reference_metadata(
  gse261983_fixture,
  "GSE261983"
)
stopifnot(
  identical(
    gse261983_meta$Original_Specimen_ID,
    c("RT00372N", "RT00374N")
  ),
  identical(
    gse261983_meta$Original_Library_ID,
    c("GSM1_RT00372N", "GSM2_RT00374N")
  ),
  all(
    gse261983_meta$Library_ID_Source_Column == "Source_GEO_Record_ID"
  ),
  all(grepl(
    "^record only:",
    gse261983_meta$Technical_Batch_Use_Status
  ))
)

gse67835_donor_counts <- c(
  AB_S1 = 24L, AB_S11 = 63L, AB_S2 = 5L, AB_S3 = 4L,
  AB_S4 = 77L, AB_S5 = 44L, AB_S7 = 57L, AB_S8 = 58L,
  FB_S1 = 26L, FB_S2 = 46L, FB_S3 = 33L, FB_S6 = 29L
)
gse67835_donor <- rep(
  names(gse67835_donor_counts),
  gse67835_donor_counts
)
gse67835_tissue <- rep("cortex", length(gse67835_donor))
for (donor_value in c("AB_S5", "AB_S8")) {
  donor_index <- which(gse67835_donor == donor_value)
  gse67835_tissue[head(donor_index, length(donor_index) %/% 3L)] <-
    "hippocampus"
}
gse67835_library <- paste0("GSM", seq_along(gse67835_donor))
gse67835_fixture <- data.frame(
  Cells = paste0("gse67835-", seq_along(gse67835_donor)),
  Sample = gse67835_donor,
  Sample_ID = gse67835_library,
  Original_Donor_ID = gse67835_donor,
  Original_Specimen_ID = paste(
    gse67835_donor,
    gse67835_tissue,
    sep = "::"
  ),
  Source_GEO_Record_ID = gse67835_library,
  Source_Single_Cell_Library_ID = gse67835_library,
  Source_Capture_Device_ID = rep(
    paste0("chip-", seq_len(16L)),
    length.out = length(gse67835_donor)
  ),
  Source_Instrument_Model = ifelse(
    gse67835_donor %in% names(gse67835_donor_counts)[seq_len(6L)],
    "Illumina MiSeq",
    "Illumina NextSeq 500"
  ),
  Source_Full_Series_Cell_Count = 466L,
  Source_Retained_Cell_Count = 466L,
  Source_Cell_Type = "neurons",
  CellType_raw = "neurons",
  Age = ifelse(startsWith(gse67835_donor, "FB_"), "14-16 PCW", "50 years"),
  Brain_Region = ifelse(
    gse67835_tissue == "cortex",
    "Cerebral cortex",
    "Hippocampus"
  ),
  Sequencing_Platform = ifelse(
    gse67835_donor %in% names(gse67835_donor_counts)[seq_len(6L)],
    "Illumina MiSeq",
    "Illumina NextSeq 500"
  ),
  stringsAsFactors = FALSE,
  row.names = paste0("gse67835-", seq_along(gse67835_donor))
)
gse67835_meta <- add_age_schema(
  upgrade_existing_reference_metadata(
    gse67835_fixture,
    "GSE67835"
  )
)
gse67835_batch <- assign_integration_batch(
  gse67835_meta,
  model = "verified_technical"
)
stopifnot(
  length(unique(gse67835_meta$Original_Donor_ID)) == 12L,
  length(unique(gse67835_meta$Original_Specimen_ID)) == 14L,
  length(unique(gse67835_meta$Original_Library_ID)) == 466L,
  identical(
    unname(unique(gse67835_meta$Stage[startsWith(gse67835_meta$Original_Donor_ID, "FB_")])),
    "S4"
  ),
  all(gse67835_batch$Technical_Batch_Dataset_Eligible),
  setequal(
    unique(gse67835_batch$Integration_Batch_ID),
    paste0(
      "GSE67835:technical_batch:",
      c("Illumina MiSeq", "Illumina NextSeq 500")
    )
  )
)

gse296073_fixture <- data.frame(
  Cells = c("gse296073-a", "gse296073-b"),
  Sample = c("h2023001", "h2023002"),
  Sample_ID = c("h2023001whgw22_all1", "h2023002whgw24_all1"),
  libraryID = c("h2023001whgw22_all1", "h2023002whgw24_all1"),
  CellType_raw = c("Microglia", "Microglia"),
  Age = c("22 PCW", "24 PCW"),
  BrainRegion = "Cortex",
  stringsAsFactors = FALSE,
  row.names = c("gse296073-a", "gse296073-b")
)
gse296073_meta <- upgrade_existing_reference_metadata(
  gse296073_fixture,
  "GSE296073"
)
stopifnot(
  identical(
    gse296073_meta$Original_Library_ID,
    gse296073_fixture$libraryID
  ),
  all(gse296073_meta$Library_ID_Source_Column == "libraryID"),
  all(grepl(
    "^record only:",
    gse296073_meta$Technical_Batch_Use_Status
  ))
)

gse97942_cells <- paste0("gse97942-", seq_len(46L))
gse97942_experiment_index <- rep(seq_len(20L), length.out = 46L)
gse97942_donor <- paste0(
  "donor-",
  ((gse97942_experiment_index - 1L) %% 6L) + 1L
)
gse97942_fixture <- data.frame(
  Cells = gse97942_cells,
  Sample = gse97942_donor,
  Sample_ID = paste0("GSM", seq_len(46L)),
  Source_GEO_Record_ID = paste0("GSM", seq_len(46L)),
  Source_Experiment_ID = paste0("experiment-", gse97942_experiment_index),
  Source_Donor_ID = gse97942_donor,
  priCluster = "Neuron",
  Age = "35 years",
  Brain_Region = "Visual cortex",
  Region = "Visual cortex",
  stringsAsFactors = FALSE,
  row.names = gse97942_cells
)
gse97942_meta <- upgrade_existing_reference_metadata(
  gse97942_fixture,
  "GSE97942"
)
stopifnot(
  length(unique(gse97942_meta$Original_Donor_ID)) == 6L,
  length(unique(gse97942_meta$Original_Library_ID)) == 46L,
  length(unique(gse97942_meta$Original_Technical_Batch_ID)) == 20L,
  all(gse97942_meta$Library_ID_Source_Column == "Source_GEO_Record_ID"),
  all(grepl(
    "^record only:",
    gse97942_meta$Technical_Batch_Use_Status
  ))
)

prjca_cells <- paste0("prjca-", seq_len(6L))
prjca_donors <- paste0("HM", seq_len(6L))
prjca_tech <- c(rep("scMultiome", 4L), rep("snRNA-seq", 2L))
prjca_fixture <- data.frame(
  Cells = prjca_cells,
  Sample = prjca_donors,
  Sample_ID = paste(prjca_donors, prjca_tech, sep = "::"),
  Source_RNA_Library_ID = paste(prjca_donors, prjca_tech, sep = "::"),
  Source_Scrublet_Score = rep(0.01, 6L),
  Source_Scrublet_Class = rep("singlet", 6L),
  Tech = prjca_tech,
  BigCellType = "EX",
  Annotation = "EX IT",
  Age = "50 years",
  Brain_Region = "Anterior cingulate cortex",
  Region = "Anterior cingulate cortex",
  stringsAsFactors = FALSE,
  row.names = prjca_cells
)
prjca_meta <- upgrade_existing_reference_metadata(
  prjca_fixture,
  "PRJCA015229"
)
stopifnot(
  length(unique(prjca_meta$Original_Donor_ID)) == 6L,
  length(unique(prjca_meta$Original_Library_ID)) == 6L,
  all(prjca_meta$Library_ID_Source_Column == "Source_RNA_Library_ID"),
  all(grepl("^record only:", prjca_meta$Technical_Batch_Use_Status))
)

gse212606_cells <- c(
  paste0("Hippocampus_1_01.", c(
    "AAAAAAAAAAAAAAAAAAAA", "CCCCCCCCCCCCCCCCCCCC"
  )),
  paste0("Hippocampus_1_02.", c(
    "GGGGGGGGGGGGGGGGGGGG", "TTTTTTTTTTTTTTTTTTTT"
  ))
)
gse212606_fixture <- data.frame(
  Cells = gse212606_cells,
  Sample = rep(c("1247", "1304"), 2L),
  Sample_ID = rep(c("1247", "1304"), 2L),
  EasySci_PCR_Group_ID = sub("[.][^.]+$", "", gse212606_cells),
  CellType_raw = c("Astrocytes", "Microglia", "Astrocytes", "Microglia"),
  Diagnosis_raw = "WT",
  Technology = "EasySci-RNA",
  Sequence = "scRNA-seq",
  Assay_Type = "scRNA-seq",
  Sequencing_Platform = NA_character_,
  Library_Chemistry = NA_character_,
  Age = rep(c("94", "81"), 2L),
  Brain_Region = "Hippocampus",
  stringsAsFactors = FALSE,
  row.names = gse212606_cells
)
gse212606_meta <- upgrade_existing_reference_metadata(
  gse212606_fixture,
  "GSE212606"
)
stopifnot(
  all(gse212606_meta$Original_Library_ID == "GSM6657986"),
  all(gse212606_meta$Library_ID_Source_Column == "GEO accession constant"),
  identical(
    gse212606_meta$Original_Technical_Batch_ID,
    gse212606_fixture$EasySci_PCR_Group_ID
  ),
  all(
    gse212606_meta$Technical_Batch_Source_Column ==
      "EasySci_PCR_Group_ID"
  ),
  all(grepl(
    "^eligible",
    gse212606_meta$Technical_Batch_Use_Status
  )),
  all(gse212606_meta$Technology == "EasySci-RNA"),
  all(gse212606_meta$Sequence == "snRNA-seq"),
  all(gse212606_meta$Assay_Type == "snRNA-seq"),
  all(gse212606_meta$Sequencing_Platform == "Illumina NovaSeq 6000"),
  all(
    gse212606_meta$Library_Chemistry ==
      "EasySci-RNA combinatorial indexing workflow"
  )
)

rosmap_fixture <- data.frame(
  Cells = sprintf("rosmap-%d", 1:4),
  Sample = c("D1", "D2", "D3", "D4"),
  Sample_ID = c("D1", "D2", "D3", "D4"),
  Batch = c("run-a", "run-a", "run-b", "run-b"),
  Diagnosis_RNA_Group_raw = rep("nonAD", 4L),
  Diagnosis_Pathology_raw = rep("no", 4L),
  Celltype = c("Ast", "Ast", "Mic", "Mic"),
  Subcelltype = c("Ast1", "Ast2", "Mic1", "Mic2"),
  Age = c("70", "71", "72", "73"),
  BrainRegion = rep("Prefrontal cortex", 4L),
  stringsAsFactors = FALSE,
  row.names = sprintf("rosmap-%d", 1:4)
)
rosmap_meta <- upgrade_existing_reference_metadata(
  rosmap_fixture,
  "ROSMAP"
)
stopifnot(
  identical(
    rosmap_meta$Technical_Batch_ID,
    c(
      "ROSMAP:technical_batch:run-a", "ROSMAP:technical_batch:run-a",
      "ROSMAP:technical_batch:run-b", "ROSMAP:technical_batch:run-b"
    )
  ),
  all(grepl(
    "^verified:",
    rosmap_meta$Technical_Batch_Verification_Status
  )),
  all(grepl(
    "^eligible",
    rosmap_meta$Technical_Batch_Use_Status
  )),
  all(is.na(rosmap_meta$Library_ID)),
  identical(
    rosmap_meta$source_cell_type_level_2,
    c("Ast1", "Ast2", "Mic1", "Mic2")
  )
)

gse199762_fixture <- data.frame(
  Cells = paste0("gse199762-", 1:9),
  Sample = c("H69", "H71", "H31", "H37", "H48", "H39", "H46", "H29", "H33"),
  Sample_ID = paste0("GSM", 1:9),
  all.exp_type = "Neuron",
  Age = c("21 PCW", "0.04 years", "13 years", "27 years", "2 years", "3 years", "0.09 years", "0.15 years", "0.15 years"),
  BrainRegion = "Entorhinal cortex",
  stringsAsFactors = FALSE,
  row.names = paste0("gse199762-", 1:9)
)
gse199762_meta <- upgrade_existing_reference_metadata(
  gse199762_fixture,
  "GSE199762"
)
gse199762_batch <- assign_integration_batch(
  gse199762_meta,
  model = "verified_technical"
)
stopifnot(
  length(unique(gse199762_meta$Technical_Batch_ID)) == 2L,
  all(grepl("^verified:", gse199762_meta$Technical_Batch_Verification_Status)),
  all(grepl("^record only:", gse199762_meta$Technical_Batch_Use_Status)),
  all(gse199762_batch$Integration_Batch_ID == "dataset:GSE199762"),
  !any(gse199762_batch$Technical_Batch_Dataset_Eligible)
)

soma_source <- data.frame(
  Cells = paste0("soma-", seq_len(6L)),
  Dataset = "SomaMut",
  Sample = c(
    "0950_240109", "5572_240109", "3848_PFC_210601",
    "5087_240109", "1465_NeuN_broad_ADpaper",
    "1465_dapi_ADpaper"
  ),
  Sample_ID = c(
    "0950_240109", "5572_240109", "3848_PFC_210601",
    "5087_240109", "1465_NeuN_broad_ADpaper",
    "1465_dapi_ADpaper"
  ),
  orig.ident = c("0950", "5572", "3848", "5087", "1465", "1465"),
  case = c(
    "0950_240109", "5572_240109", "3848_PFC_210601",
    "5087_240109", "1465_NeuN_broad_ADpaper",
    "1465_dapi_ADpaper"
  ),
  batch = c(
    "240109", "240109", "210601", "240430", "ADpaper", "ADpaper"
  ),
  Age = c(53, 70, 38, 57, 17, 17),
  Region = "Prefrontal cortex",
  CellType_raw = "Oligodendrocytes-1",
  new_clusters3 = "Oligodendrocytes",
  new_clusters2 = "Oligodendrocytes",
  new_clusters = "Oligodendrocytes-1",
  predicted.id = "Oligodendrocytes",
  prediction.score.max = 0.9,
  stringsAsFactors = FALSE
)
soma_meta <- upgrade_existing_reference_metadata(soma_source, "SomaMut")
stopifnot(
  identical(
    soma_meta$Original_Donor_ID,
    c("0950", "5572", "3848", "5087", "1465", "1465")
  ),
  identical(soma_meta$Original_Library_ID, soma_source$Sample_ID),
  identical(
    soma_meta$Original_Technical_Batch_ID,
    c(
      "240109", "240109", "210601", "240430",
      "ADpaper", "ADpaper"
    )
  ),
  all(soma_meta$Library_ID_Source_Column == "case"),
  all(soma_meta$Technical_Batch_Source_Column == "batch"),
  all(soma_meta$source_cell_type_level_1 == "Oligodendrocytes"),
  all(soma_meta$source_cell_type_level_2 == "Oligodendrocytes"),
  all(soma_meta$source_cell_type_level_3 == "Oligodendrocytes-1")
)
soma_batch <- assign_integration_batch(
  soma_meta,
  model = "verified_technical"
)
stopifnot(
  all(soma_batch$Integration_Batch_ID == "dataset:SomaMut"),
  !any(soma_batch$Technical_Batch_Dataset_Eligible)
)

donor_crosswalk <- file.path(test_dir, "donor_crosswalk.tsv")
write.table(
  data.frame(
    Dataset = "Fixture_A",
    Original_Donor_ID = "donor-04",
    Global_Donor_ID = "Fixture_A:donor:donor-04",
    Duplicate_Group = "fixture-duplicate",
    Duplicate_Evidence = "synthetic exact donor match",
    Primary_Analysis_Include = "false",
    Analysis_Role = "reference_excluded",
    Exclusion_Reason = "duplicate donor",
    stringsAsFactors = FALSE
  ),
  donor_crosswalk,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
crosswalk_metadata <- apply_donor_crosswalk(metadata, donor_crosswalk)
stopifnot(
  sum(crosswalk_metadata$Analysis_Include) == 90L,
  all(!crosswalk_metadata$Analysis_Include[91:120]),
  all(crosswalk_metadata$Duplicate_Group[91:120] == "fixture-duplicate")
)

external_matrix_fixture <- data.frame(
  Original_Cell_ID = paste0("external-cell-", 1:4),
  SampleID = rep(c("donor-1", "donor-2"), each = 2L),
  stringsAsFactors = FALSE
)
external_donor_fixture <- data.frame(
  Source_Donor_ID = c("donor-1", "donor-2", "source-only-donor"),
  Source_Age = c(15, 40, 60),
  stringsAsFactors = FALSE
)
external_joined_fixture <- join_external_metadata(
  matrix_meta = external_matrix_fixture,
  external_meta = external_donor_fixture,
  key_column = "Source_Donor_ID",
  matrix_key_column = "SampleID",
  match_column = "Donor_Metadata_Matched",
  require_all_matrix_cells = TRUE,
  require_all_metadata_cells = FALSE
)
external_join_audit <- attr(
  external_joined_fixture,
  "external_metadata_audit"
)
stopifnot(
  identical(external_joined_fixture$Source_Age, c(15, 15, 40, 40)),
  all(external_joined_fixture$Donor_Metadata_Matched),
  external_join_audit$matrix_key_column[[1L]] == "SampleID",
  external_join_audit$external_key_column[[1L]] == "Source_Donor_ID",
  external_join_audit$metadata_cells_without_matrix[[1L]] == 1L
)

# A complete two-shard fixture exercises the no-loss invariant: every source
# index appears exactly once and every shard is a self-contained dgCMatrix.
if (!requireNamespace("BPCells", quietly = TRUE)) {
  stop("BPCells is required for the complete layered-processing fixture")
}
layered_dir <- file.path(test_dir, "layered")
dir.create(layered_dir)
counts_h5 <- BPCells::write_matrix_memory(
  BPCells::convert_matrix_type(counts, "uint32_t")
)
layered_config <- list(
  matrix_storage = "dgCMatrix_shards",
  layer_by = "Original_Library_ID",
  expected_layers = 2L,
  expected_cells = 120L,
  expected_features = 40L,
  expected_analysis_cells = 120L,
  expected_analysis_donors = 4L,
  matrix_group = "X"
)
metadata$Float_Roundtrip <- rep(
  c(0.123456701, 0.987654328),
  length.out = nrow(metadata)
)
log_message <- function(message, message_type = "info") {
  invisible(list(message = message, message_type = message_type))
}
process_layered_h5ad(
  dataset = "Fixture_A",
  config = layered_config,
  counts_h5 = counts_h5,
  metadata = metadata,
  feature_meta = data.frame(Original_Feature_ID = features),
  input_cells = cells,
  input_features = features,
  processed_dir = layered_dir,
  canonical_metadata_file = file.path(layered_dir, "metadata_canonical.tsv.gz"),
  processing_audit_file = file.path(layered_dir, "processing_audit.tsv"),
  object_format_file = file.path(layered_dir, "processed_object_format.tsv"),
  full_object_file = file.path(layered_dir, "unused_full.rds"),
  analysis_object_file = file.path(layered_dir, "unused_analysis.rds"),
  source_record = data.frame(
    source_accession = "Fixture_A",
    source_repository = "synthetic fixture",
    publication_doi = "10.0000/synthetic.fixture",
    verified_title = "Synthetic clean-room fixture",
    journal = "Synthetic",
    publication_year = "2026",
    repository_record_url = "https://example.org/synthetic-fixture",
    stringsAsFactors = FALSE
  )
)
layered_summary <- validate_processed_shard_bundle(
  layered_dir,
  expected_dataset = "Fixture_A"
)
layered_manifest <- read_processed_shard_manifest(layered_dir)
stopifnot(
  layered_summary$Full_Cells == 120L,
  layered_summary$Analysis_Cells == 120L,
  layered_summary$Features == 40L,
  layered_summary$Matrix_Class == "dgCMatrix_shards",
  sum(layered_manifest$Cells) == 120L,
  sum(layered_manifest$Nonzero_Values) == Matrix::nnzero(counts),
  isTRUE(all.equal(sum(layered_manifest$Total_Counts), sum(counts), tolerance = 0))
)
layered_metadata_revision <- data.table::fread(
  file.path(layered_dir, "metadata_canonical.tsv.gz"),
  data.table = FALSE
)
layered_metadata_revision$Metadata_Revision <- "fixture metadata refresh"
layered_refresh <- refresh_processed_shard_metadata(
  processed_dir = layered_dir,
  metadata = layered_metadata_revision,
  expected_dataset = "Fixture_A"
)
layered_refreshed_objects <- load_processed_shards(
  layered_dir,
  verify_sha256 = TRUE
)
stopifnot(
  layered_refresh$Shards == 2L,
  layered_refresh$Full_Cells == 120L,
  !layered_refresh$Matrix_Recomputed,
  all(vapply(
    layered_refreshed_objects,
    function(object) {
      all(object$Metadata_Revision == "fixture metadata refresh")
    },
    logical(1)
  )),
  sum(vapply(
    layered_refreshed_objects,
    function(object) sum(processed_counts(object)),
    numeric(1)
  )) == sum(counts)
)
rm(layered_refreshed_objects)

# Seurat v5 keeps merged per-dataset count layers as ordered feature subsets
# when source datasets do not share every feature. The processed-object gate
# must validate the union and order without requiring artificial zero padding.
heterogeneous_a <- Matrix::Matrix(
  matrix(c(1, 0, 2, 0, 3, 0), nrow = 3L),
  sparse = TRUE
)
rownames(heterogeneous_a) <- c("GENE1", "GENE2", "GENE3")
colnames(heterogeneous_a) <- c("heterogeneous_a_1", "heterogeneous_a_2")
heterogeneous_b <- Matrix::Matrix(
  matrix(c(4, 0, 0, 5), nrow = 2L),
  sparse = TRUE
)
rownames(heterogeneous_b) <- c("GENE2", "GENE4")
colnames(heterogeneous_b) <- c("heterogeneous_b_1", "heterogeneous_b_2")
heterogeneous_merged <- merge(
  SeuratObject::CreateSeuratObject(heterogeneous_a, project = "a"),
  SeuratObject::CreateSeuratObject(heterogeneous_b, project = "b"),
  merge.data = FALSE
)
validate_processed_object(
  heterogeneous_merged,
  expected_features = 4L,
  expected_cells = 4L
)
heterogeneous_summary <- processed_object_content_summary(
  heterogeneous_merged
)
stopifnot(
  heterogeneous_summary$count_layers == 2L,
  heterogeneous_summary$features == 4L,
  heterogeneous_summary$cells == 4L,
  heterogeneous_summary$total_counts ==
    sum(heterogeneous_a) + sum(heterogeneous_b)
)
rm(
  heterogeneous_a,
  heterogeneous_b,
  heterogeneous_merged,
  heterogeneous_summary
)

# Standardize one representative original dataset through the production
# entry point. This proves sidecar coverage without changing the retained RDS.
existing_root <- file.path(test_dir, "existing")
existing_dataset <- "GSE217511"
existing_dir <- file.path(existing_root, existing_dataset)
dir.create(existing_dir, recursive = TRUE)
existing_counts <- counts[seq_len(12L), seq_len(8L), drop = FALSE]
existing_meta <- data.frame(
  Cells = colnames(existing_counts),
  Dataset = existing_dataset,
  Technology = "10X Genomics",
  Sequence = "snRNA-seq",
  Sample = rep(c("2C", "2G", "3C", "3G"), each = 2L),
  Sample_ID = rep(sprintf("GSM%d", 1:4), each = 2L),
  CellType_raw = rep(c("RG", "Neuron"), 4L),
  Brain_Region = rep(c("Cortical plate", "Germinal matrix"), 4L),
  Region = rep(c("Cortical plate", "Germinal matrix"), 4L),
  Age = rep(c("16-18 PCW", "20 PCW"), each = 4L),
  Sex = rep(c("Female", "Male"), each = 4L),
  row.names = colnames(existing_counts),
  stringsAsFactors = FALSE
)
existing_object <- SeuratObject::CreateSeuratObject(
  counts = existing_counts,
  meta.data = existing_meta
)
existing_file <- file.path(
  existing_dir,
  paste0(existing_dataset, "_processed.rds")
)
save_processed_object(existing_object, existing_file)
existing_sha256 <- processed_file_sha256(existing_file)
stale_sidecar <- file.path(existing_dir, "metadata_raw.tsv.gz")
writeLines("stale metadata placeholder", stale_sidecar)
bundle_manifest <- file.path(existing_dir, ".processed_bundle.sha256")
writeLines(
  c(
    paste(existing_sha256, basename(existing_file), sep = "  "),
    paste(
      processed_file_sha256(stale_sidecar),
      basename(stale_sidecar),
      sep = "  "
    )
  ),
  bundle_manifest
)
empty_crosswalk <- file.path(test_dir, "empty_donor_crosswalk.tsv")
write.table(
  data.frame(
    Dataset = character(),
    Original_Donor_ID = character(),
    Global_Donor_ID = character(),
    Duplicate_Group = character(),
    Duplicate_Evidence = character(),
    Primary_Analysis_Include = character(),
    Analysis_Role = character(),
    Exclusion_Reason = character(),
    stringsAsFactors = FALSE
  ),
  empty_crosswalk,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
standardize_status <- system2(
  file.path(R.home("bin"), "Rscript"),
  c(
    "processing/standardize_datasets.R",
    "--dataset", existing_dataset,
    "--processed-root", existing_root,
    "--donor-crosswalk", empty_crosswalk
  )
)
stopifnot(
  standardize_status == 0L,
  processed_file_sha256(existing_file) == existing_sha256,
  all(file.exists(file.path(
    existing_dir,
    c(
      "metadata_raw.tsv.gz", "metadata_canonical.tsv.gz",
      "features_raw.tsv.gz", "age_crosswalk.tsv",
      "region_crosswalk.tsv", "source_cell_type_crosswalk.tsv",
      "technical_batch_crosswalk.tsv",
      "technical_batch_design_audit.tsv",
      "source_technical_field_inventory.tsv",
      "donor_specimen_library_crosswalk.tsv.gz",
      "processing_audit.tsv", "processed_object_format.tsv"
    )
  )))
)
existing_bundle <- validate_processed_object_bundle(
  existing_dir,
  expected_dataset = existing_dataset
)
bundle_lines <- readLines(bundle_manifest, warn = FALSE)
bundle_paths <- sub(
  "^[0-9a-f]{64}[[:space:]]+[*]?",
  "",
  bundle_lines
)
bundle_checksums <- sub(
  "^([0-9a-f]{64})[[:space:]]+[*]?.+$",
  "\\1",
  bundle_lines
)
names(bundle_checksums) <- bundle_paths
stopifnot(
  all(c(
    basename(existing_file),
    "metadata_raw.tsv.gz",
    "processed_object_format.tsv"
  ) %in% bundle_paths),
  all(vapply(
    bundle_paths,
    function(path) {
      identical(
        unname(bundle_checksums[[path]]),
        processed_file_sha256(file.path(existing_dir, path))
      )
    },
    logical(1)
  )),
  identical(
    unname(bundle_checksums[[basename(existing_file)]]),
    processed_file_sha256(existing_file)
  ),
  identical(
    unname(bundle_checksums[[basename(stale_sidecar)]]),
    processed_file_sha256(stale_sidecar)
  )
)
existing_canonical <- read_brainomics_table(file.path(
  existing_dir,
  "metadata_canonical.tsv.gz"
))
existing_source_cell_type <- read_brainomics_table(file.path(
  existing_dir,
  "source_cell_type_crosswalk.tsv"
))
existing_technical_design <- read_brainomics_table(file.path(
  existing_dir,
  "technical_batch_design_audit.tsv"
))
existing_technical_inventory <- read_brainomics_table(file.path(
  existing_dir,
  "source_technical_field_inventory.tsv"
))
stopifnot(
  existing_bundle$Full_Cells == 8L,
  existing_bundle$Features == 12L,
  existing_bundle$Matrix_Class == "dgCMatrix",
  all(existing_canonical$Source_Accession == "GSE217511"),
  all(existing_canonical$Source_Publication_DOI ==
    "10.1038/s41467-022-34975-2"),
  all(existing_canonical$source_cell_type_original_label %in%
    c("RG", "Neuron")),
  identical(
    unique(existing_canonical$Original_Library_ID),
    sprintf("GSM%d", 1:4)
  ),
  all(existing_canonical$Library_ID_Source_Column == "Sample_ID"),
  all(grepl(
    "one-to-one",
    existing_canonical$Library_ID_Derivation_Rule,
    fixed = TRUE
  )),
  all(c(
    "Source_Publication_Title", "source_cell_type_original_label",
    "source_cell_type_original_label_source_column"
  ) %in% names(existing_source_cell_type)),
  nrow(existing_technical_design) >= 1L,
  existing_technical_design$Record_Type[[1L]] == "dataset_summary",
  nrow(existing_technical_inventory) >= 1L,
  existing_technical_inventory$Record_Type[[1L]] == "dataset_summary"
)

# Explicit method identity prevents the historical Harmony-as-RPCA relabeling.
stopifnot(
  length(reference_datasets()) == 22L,
  identical(
    excluded_reference_datasets(),
    c(
      "GSE103723", "GSE199762", "GSE261983",
      "Nowakowski_et_al_2017"
    )
  ),
  identical(query_validation_datasets(), "EGAD00001006049"),
  length(formal_datasets()) == 23L,
  !"Catching_2026" %in% formal_datasets()
)

# Complete-reference RPCA uses future.apply even under a sequential plan. A
# tiny globals limit reproduces the production guard failure; the formal
# configuration must remove that guard only after enforcing one in-process
# worker.
future::plan(future::sequential)
options(future.globals.maxSize = 1)
future_fixture_global <- raw(1024L)
future_guard_failure <- try(
  future.apply::future_lapply(
    1L,
    function(index) length(future_fixture_global) + index
  ),
  silent = TRUE
)
stopifnot(inherits(future_guard_failure, "try-error"))
future_configuration <- configure_complete_integration_future()
future_guard_result <- future.apply::future_lapply(
  1L,
  function(index) length(future_fixture_global) + index
)
stopifnot(
  identical(future_configuration$Plan, "sequential"),
  identical(as.integer(future_configuration$Workers), 1L),
  is.infinite(getOption("future.globals.maxSize")),
  identical(future_guard_result, list(1025L))
)

# The R-to-Python scVI boundary must export every reference cell exactly once
# without sampling or materializing one global sparse matrix.
scvi_fixture_features <- paste0("scviGene", seq_len(10L))
scvi_fixture_objects <- lapply(
  seq_along(reference_datasets()),
  function(index) {
    fixture_counts <- Matrix::sparseMatrix(
      i = rep(seq_len(10L), 5L),
      j = rep(seq_len(5L), each = 10L),
      x = rep(as.numeric(index), 50L),
      dims = c(10L, 5L),
      dimnames = list(
        scvi_fixture_features,
        paste0(reference_datasets()[[index]], "_scvi_cell_", seq_len(5L))
      )
    )
    fixture_object <- SeuratObject::CreateSeuratObject(fixture_counts)
    fixture_object$Integration_Batch_ID <- paste0(
      "dataset:",
      reference_datasets()[[index]]
    )
    fixture_object
  }
)
names(scvi_fixture_objects) <- reference_datasets()
scvi_fixture_contract <- scvi_source_matrix_contract()
scvi_fixture_contract$Source_Matrix_Units[] <-
  "raw non-negative integer gene counts"
scvi_fixture_contract$Source_Provenance_Eligible[] <- TRUE
scvi_fixture_contract$Source_Evidence[] <- "synthetic integer-count fixture"
scvi_fixture_eligibility <- audit_scvi_input_eligibility(
  scvi_fixture_objects,
  source_contract = scvi_fixture_contract
)
scvi_fixture_dir <- file.path(test_dir, "scvi_input")
scvi_fixture_cells <- unlist(lapply(
  scvi_fixture_objects,
  colnames
), use.names = FALSE)
scvi_fixture_audit <- export_scvi_input_shards(
  objects_list = scvi_fixture_objects,
  features = scvi_fixture_features,
  output_dir = scvi_fixture_dir,
  eligibility_audit = scvi_fixture_eligibility,
  expected_cells = scvi_fixture_cells
)
scvi_fixture_manifest <- read.delim(
  file.path(scvi_fixture_dir, "manifest.tsv"),
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
stopifnot(
  scvi_fixture_audit$Datasets[[1L]] == 22L,
  scvi_fixture_audit$Cells[[1L]] == 110L,
  scvi_fixture_audit$Features[[1L]] == 10L,
  nrow(scvi_fixture_manifest) == 22L,
  sum(scvi_fixture_manifest$Cells) == 110L,
  all(scvi_fixture_manifest$Nonzero_Values == 50L)
)

scvi_fractional_objects <- scvi_fixture_objects
fractional_counts <- LayerData(
  scvi_fractional_objects[["AllenM1"]],
  assay = "RNA",
  layer = "counts"
)
fractional_counts@x[[1L]] <- 0.5
scvi_fractional_objects[["AllenM1"]][["RNA"]]$counts <-
  fractional_counts
scvi_fractional_eligibility <- audit_scvi_input_eligibility(
  scvi_fractional_objects,
  source_contract = scvi_fixture_contract
)
stopifnot(
  !scvi_fractional_eligibility$Numeric_Count_Eligible[
    scvi_fractional_eligibility$Dataset == "AllenM1"
  ],
  !scvi_fractional_eligibility$ScVI_Eligible[
    scvi_fractional_eligibility$Dataset == "AllenM1"
  ],
  scvi_fractional_eligibility$Fractional_Values[
    scvi_fractional_eligibility$Dataset == "AllenM1"
  ] == 1
)

method_identity <- integration_method_identity()
reduction_contract_file <- file.path(
  test_dir,
  "integration_reduction_audit.tsv"
)
fixture_cell_hash <- paste(rep("a", 64L), collapse = "")
reduction_contract_fixture <- transform(
  method_identity,
  Cells = 120L,
  Dimensions = ifelse(Space == "latent", 5L, 2L),
  Cell_Order_SHA256 = fixture_cell_hash,
  Finite = TRUE
)
write.table(
  reduction_contract_fixture,
  reduction_contract_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
reduction_contract <- read_integration_reduction_cell_contract(
  reduction_contract_file,
  expected_cells = 120L
)
stopifnot(
  identical(reduction_contract$Cells, 120),
  identical(reduction_contract$Cell_Order_SHA256, fixture_cell_hash),
  identical(
    reduction_contract$Audit$Output_Name,
    method_identity$Output_Name
  )
)
invalid_reduction_contract <- reduction_contract_fixture
invalid_reduction_contract$Cell_Order_SHA256[[1L]] <- "invalid"
write.table(
  invalid_reduction_contract,
  reduction_contract_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
stopifnot(inherits(
  try(
    read_integration_reduction_cell_contract(reduction_contract_file),
    silent = TRUE
  ),
  "try-error"
))
batch_method_contract <- integration_batch_method_contract()
validate_integration_batch_method_contract(batch_method_contract)
stopifnot(
  identical(
    batch_method_contract$Method,
    rep(method_levels[-1L], 2L)
  ),
  identical(
    batch_method_contract$Batch_Model,
    rep(c("dataset", "verified_technical"), each = 3L)
  ),
  all(
    batch_method_contract$Required_Metadata_Key ==
      "Integration_Batch_ID"
  )
)
embedding_data <- list(
  pca_raw = matrix(1, nrow = 120L, ncol = 5L, dimnames = list(cells, NULL)),
  pca_rpca = matrix(2, nrow = 120L, ncol = 5L, dimnames = list(cells, NULL)),
  pca_harmony = matrix(3, nrow = 120L, ncol = 5L, dimnames = list(cells, NULL)),
  pca_scvi = matrix(4, nrow = 120L, ncol = 5L, dimnames = list(cells, NULL)),
  umap_raw = matrix(11, nrow = 120L, ncol = 2L, dimnames = list(cells, NULL)),
  umap_rpca = matrix(12, nrow = 120L, ncol = 2L, dimnames = list(cells, NULL)),
  umap_harmony = matrix(13, nrow = 120L, ncol = 2L, dimnames = list(cells, NULL)),
  umap_scvi = matrix(14, nrow = 120L, ncol = 2L, dimnames = list(cells, NULL))
)
lisi_meta <- data.frame(
  Dataset = rep(c("Fixture_A", "Fixture_B"), each = 60L),
  row.names = cells,
  stringsAsFactors = FALSE
)
fake_lisi <- function(X, meta_data, label_colnames, perplexity, nn_eps) {
  stopifnot(
    identical(rownames(X), rownames(meta_data)),
    identical(label_colnames, "Dataset"),
    identical(perplexity, 30),
    identical(nn_eps, 0)
  )
  data.frame(Dataset = rowMeans(X))
}
lisi_results <- compute_lisi_results(
  embedding_data,
  lisi_meta,
  method_identity = method_identity,
  lisi_fun = fake_lisi
)
expected_lisi_contract <- paste(
  names(lisi_numeric_contract()),
  lisi_numeric_contract(),
  sep = "="
)
stopifnot(isTRUE(validate_lisi_numeric_contract(expected_lisi_contract)))
invalid_lisi_contract <- expected_lisi_contract
invalid_lisi_contract[[4L]] <- "numeric_patch=single_precision"
stopifnot(inherits(
  try(validate_lisi_numeric_contract(invalid_lisi_contract), silent = TRUE),
  "try-error"
))
bounded_lisi <- canonicalize_lisi_lower_bound(c(1 - 1e-14, 1, 2))
stopifnot(
  identical(as.numeric(bounded_lisi), c(1, 1, 2)),
  attr(bounded_lisi, "theoretical_lower_bound_adjustments") == 1L,
  attr(bounded_lisi, "pre_adjustment_minimum") == 1 - 1e-14,
  inherits(
    try(canonicalize_lisi_lower_bound(c(0.99, 1)), silent = TRUE),
    "try-error"
  )
)
stopifnot(
  identical(names(lisi_results), method_levels),
  identical(attr(lisi_results, "space"), "latent"),
  identical(attr(lisi_results, "nn_eps"), 0),
  all(lisi_results$Raw == 1),
  all(lisi_results$RPCA == 2),
  all(lisi_results$Harmony == 3),
  all(lisi_results$scVI == 4)
)

# Export the deterministic mini-atlas and read it back through the generated R
# reader. Python round-trip runs when its documented dependencies are present.
release_meta <- data.frame(
  Cells = cells,
  Dataset = lisi_meta$Dataset,
  Technology = "10x Genomics",
  Sequence = "snRNA-seq",
  Sample = rep(sprintf("donor-%02d", 1:4), each = 30L),
  Sample_ID = rep(sprintf("library-%02d", 1:8), each = 15L),
  Library_ID = c(
    rep(NA_character_, 10L),
    rep(sprintf("library-%02d", 1:8), each = 15L)[11:120]
  ),
  BrainRegion = rep(c("Prefrontal cortex", "Cerebellum"), each = 60L),
  Age = rep(c("8 PCW", "10 PCW", "6 years", "20 years"), each = 30L),
  Sex = rep(c("Female", "Male"), 60L),
  Stage = rep(c("S2", "S3", "S11", "S13"), each = 30L),
  CellType = rep(c("Excitatory neuron", "Oligodendrocyte"), 60L),
  row.names = cells,
  stringsAsFactors = FALSE
)
release_count_layers <- list(
  Fixture_A = counts[, seq_len(60L), drop = FALSE],
  Fixture_B = counts[, 60L + seq_len(60L), drop = FALSE]
)
release_object <- SeuratObject::CreateSeuratObject(
  counts = release_count_layers,
  meta.data = release_meta
)
make_reduction <- function(value, dimensions, key) {
  embedding <- matrix(
    value + seq_len(120L * dimensions) / 10000,
    nrow = 120L,
    ncol = dimensions,
    dimnames = list(cells, paste0(key, seq_len(dimensions)))
  )
  SeuratObject::CreateDimReducObject(
    embeddings = embedding,
    key = key,
    assay = "RNA"
  )
}
release_object[["integrated.rpca"]] <- make_reduction(2, 5L, "RPCA_")
release_object[["umap.rpca"]] <- make_reduction(3, 2L, "UMAP_")
release_object[["umap.unintegrated"]] <- make_reduction(1, 2L, "RAWUMAP_")

integration_dir <- file.path(test_dir, "integration")
dir.create(integration_dir)
object_file <- file.path(integration_dir, "objects_celltypes.rds")
lisi_file <- file.path(integration_dir, "lisi_results.rds")
saveRDS(release_object, object_file)
saveRDS(lisi_results, lisi_file)
fixture_reduction_dir <- file.path(
  integration_dir, "annotation", "reductions"
)
fixture_evaluation_dir <- file.path(integration_dir, "evaluation")
dir.create(fixture_reduction_dir, recursive = TRUE)
dir.create(fixture_evaluation_dir, recursive = TRUE)
fixture_rpca_latent <- matrix(
  2 + seq_len(120L * 50L) / 100000,
  nrow = 120L,
  ncol = 50L,
  dimnames = list(cells, paste0("RPCA_", seq_len(50L)))
)
fixture_umap_plot <- data.frame(
  Cell = cells,
  Dataset = lisi_meta$Dataset,
  Raw_1 = SeuratObject::Embeddings(
    release_object, "umap.unintegrated"
  )[, 1L],
  Raw_2 = SeuratObject::Embeddings(
    release_object, "umap.unintegrated"
  )[, 2L],
  RPCA_1 = SeuratObject::Embeddings(release_object, "umap.rpca")[, 1L],
  RPCA_2 = SeuratObject::Embeddings(release_object, "umap.rpca")[, 2L],
  stringsAsFactors = FALSE
)
saveRDS(
  fixture_rpca_latent,
  file.path(fixture_reduction_dir, "rpca_latent.rds")
)
saveRDS(
  fixture_umap_plot,
  file.path(fixture_evaluation_dir, "umap_plot_data.rds")
)
fixture_celltype_assignments <- data.frame(
  Cells = rev(cells),
  Cluster = rep(c("C00", "C01"), length.out = length(cells)),
  CellType = rep(
    c("Fixture excitatory neurons", "Fixture inhibitory neurons"),
    length.out = length(cells)
  ),
  stringsAsFactors = FALSE
)
fixture_celltype_assignment_file <- file.path(
  integration_dir, "main_celltype_assignments.rds"
)
saveRDS(fixture_celltype_assignments, fixture_celltype_assignment_file)
package_dir <- file.path(test_dir, "ScienceDB")
export_status <- system2(
  file.path(R.home("bin"), "Rscript"),
  c(
    "sciencedb/export.R",
    "--repo-dir", repo_dir,
    "--object-file", object_file,
    "--metadata-object-file", object_file,
    "--celltype-assignment-file", fixture_celltype_assignment_file,
    "--lisi-file", lisi_file,
    "--out-dir", package_dir,
    "--overwrite"
  )
)
stopifnot(export_status == 0L)
fixture_expression_files <- sort(list.files(
  file.path(package_dir, "expression"),
  recursive = TRUE,
  full.names = TRUE
))
fixture_expression_sha_before <- unname(tools::md5sum(fixture_expression_files))
reuse_export_status <- system2(
  file.path(R.home("bin"), "Rscript"),
  c(
    "sciencedb/export.R",
    "--repo-dir", repo_dir,
    "--object-file", object_file,
    "--metadata-object-file", object_file,
    "--celltype-assignment-file", fixture_celltype_assignment_file,
    "--lisi-file", lisi_file,
    "--out-dir", package_dir,
    "--reuse-expression",
    "--overwrite"
  )
)
fixture_expression_sha_after <- unname(tools::md5sum(fixture_expression_files))
stopifnot(
  reuse_export_status == 0L,
  identical(fixture_expression_sha_before, fixture_expression_sha_after)
)
release_readme <- readLines(
  file.path(package_dir, "README.md"),
  warn = FALSE
)
stopifnot(any(grepl(
  "`provenance/environment/`", release_readme,
  fixed = TRUE
)))

restored_file <- file.path(test_dir, "reader_roundtrip.rds")
r_reader_status <- system2(
  file.path(R.home("bin"), "Rscript"),
  c(
    file.path(package_dir, "scripts", "read_seurat.R"),
    package_dir,
    restored_file
  )
)
stopifnot(r_reader_status == 0L, file.exists(restored_file))
restored <- readRDS(restored_file)
restored_count_layers <- SeuratObject::Layers(
  restored[["RNA"]],
  search = "^counts"
)
restored_count_sum <- sum(vapply(
  restored_count_layers,
  function(layer) {
    sum(SeuratObject::LayerData(
      restored[["RNA"]],
      layer = layer, fast = TRUE
    ))
  },
  numeric(1L)
))
stopifnot(
  identical(dim(restored), dim(release_object)),
  identical(colnames(restored), cells),
  identical(rownames(restored), features),
  length(restored_count_layers) == 2L,
  isTRUE(all.equal(restored_count_sum, sum(counts), tolerance = 0)),
  all(c(
    "integrated_pca", "integrated_umap", "unintegrated_umap"
  ) %in% names(restored@reductions)),
  file.exists(file.path(package_dir, "provenance", "file_manifest.tsv")),
  file.exists(file.path(package_dir, "md5sum.txt")),
  file.exists(file.path(package_dir, "validation", "lisi.tsv.gz"))
)
release_metadata <- utils::read.delim(
  gzfile(file.path(package_dir, "metadata", "metadata.tsv.gz")),
  sep = "\t", check.names = FALSE, stringsAsFactors = FALSE
)
fixture_assignment_index <- match(
  release_metadata$Cells, fixture_celltype_assignments$Cells
)
stopifnot(
  identical(
    release_metadata$Cluster,
    fixture_celltype_assignments$Cluster[fixture_assignment_index]
  ),
  identical(
    release_metadata$CellType,
    fixture_celltype_assignments$CellType[fixture_assignment_index]
  )
)
release_manifest <- utils::read.delim(
  file.path(package_dir, "provenance", "file_manifest.tsv"),
  sep = "\t", check.names = FALSE, stringsAsFactors = FALSE
)
release_shard_manifest <- utils::read.delim(
  file.path(package_dir, "expression", "shard_manifest.tsv"),
  sep = "\t", check.names = FALSE, stringsAsFactors = FALSE
)
release_matrix_rows <- release_manifest[
  grepl("^expression/shards/.+/matrix[.]mtx[.]gz$", release_manifest$path), ,
  drop = FALSE
]
stopifnot(
  all(c(
    "path", "size_bytes", "rows", "columns", "nonzero_values",
    "schema_version", "license_scope", "md5", "sha256"
  ) %in% names(release_manifest)),
  all(nchar(release_manifest$md5) == 32L),
  all(nchar(release_manifest$sha256) == 64L),
  release_manifest$rows[
    release_manifest$path == "metadata/metadata.tsv.gz"
  ] == 120L,
  nrow(release_shard_manifest) == 2L,
  identical(release_shard_manifest$Dataset, c("Fixture_A", "Fixture_B")),
  sum(release_shard_manifest$Cells) == 120L,
  sum(release_shard_manifest$Nonzero_Values) == length(counts@x),
  nrow(release_matrix_rows) == 2L,
  all(release_matrix_rows$columns == 60L),
  sum(release_matrix_rows$nonzero_values) == length(counts@x)
)

release_metadata <- utils::read.delim(
  gzfile(file.path(package_dir, "metadata", "metadata.tsv.gz")),
  sep = "\t", check.names = FALSE, stringsAsFactors = FALSE
)
expected_release_metadata_columns <- c(
  "Cells", "Dataset", "Original_Sample_ID", "Donor_ID", "Sample_ID",
  "Library_ID", "Technology", "Modality", "Age", "AgeIntervalID",
  "Sex", "BrainRegion", "Source_CellType", "Cluster", "CellType"
)
release_dataset_summary <- utils::read.delim(
  file.path(package_dir, "metadata", "dataset_summary.tsv"),
  sep = "\t", check.names = FALSE, stringsAsFactors = FALSE
)
release_attrition_summary <- utils::read.delim(
  file.path(package_dir, "metadata", "dataset_attrition_summary.tsv"),
  sep = "\t", check.names = FALSE, stringsAsFactors = FALSE
)
release_verification_dictionary <- utils::read.delim(
  file.path(package_dir, "metadata", "verification_status_dictionary.tsv"),
  sep = "\t", check.names = FALSE, stringsAsFactors = FALSE
)
release_references <- utils::read.delim(
  file.path(package_dir, "provenance", "dataset_manifest.tsv"),
  sep = "\t", check.names = FALSE, stringsAsFactors = FALSE
)
source_input_contract <- brainomics_source_input_contract()
stopifnot(
  identical(names(release_metadata), expected_release_metadata_columns),
  identical(release_metadata$Cells, cells),
  setequal(unique(release_metadata$Age), c("8 PCW", "10 PCW", "6", "20")),
  !any(grepl("year|month|stage", release_metadata$Age, ignore.case = TRUE)),
  all(c(
    "First_Author", "Original_File_Type", "Genome_Transcriptome_Build"
  ) %in% names(release_dataset_summary)),
  "First_Author" %in% names(release_references),
  nrow(release_attrition_summary) == 2L,
  sum(release_attrition_summary$Source_Cells) == 120L,
  sum(release_attrition_summary$Exported_Cells) == 120L,
  identical(
    release_verification_dictionary$Verification_Status,
    c(
      "verified", "source-reported", "derived", "curated",
      "unresolved", "excluded", "not-applicable"
    )
  ),
  nrow(source_input_contract) == 22L,
  !anyDuplicated(source_input_contract$Dataset),
  all(c(
    "Original_File_Type", "Expression_Units",
    "Genome_Transcriptome_Build"
  ) %in% names(source_input_contract)),
  !anyNA(source_input_contract[, c(
    "Dataset", "Original_File_Type", "Expression_Units",
    "Genome_Transcriptome_Build"
  )])
)

python <- Sys.getenv("BRAINOMICS_PYTHON", "python3")
python_ready <- system2(
  python,
  c("-c", shQuote("import scanpy, pandas, anndata")),
  stdout = FALSE,
  stderr = FALSE
) == 0L
if (python_ready) {
  h5ad_dir <- file.path(test_dir, "reader_roundtrip_h5ad")
  python_status <- system2(
    python,
    c(
      file.path(package_dir, "scripts", "read_h5ad.py"),
      package_dir,
      h5ad_dir
    )
  )
  stopifnot(
    python_status == 0L,
    file.exists(file.path(h5ad_dir, "Fixture_A.h5ad")),
    file.exists(file.path(h5ad_dir, "Fixture_B.h5ad"))
  )
  stopifnot(system2(python, c("tests/test_reader_roundtrip.py", shQuote(package_dir), shQuote(h5ad_dir))) == 0L)
  message("Python_reader_status=passed")
} else {
  if (identical(Sys.getenv("BRAINOMICS_REQUIRE_PYTHON_READER"), "true")) {
    stop("Required Python reader round-trip cannot run: scanpy/pandas/anndata unavailable")
  }
  message("Python_reader_status=skipped; scanpy/pandas/anndata unavailable")
}

source_original_sample_ids <- release_metadata$Original_Sample_ID
synthetic_public_sidecar <- data.frame(
  Cells = cells,
  Global_Donor_ID = release_metadata$Donor_ID,
  Specimen_ID = release_metadata$Sample_ID,
  Library_ID = release_metadata$Library_ID,
  Original_Donor_ID = release_metadata$Donor_ID,
  stringsAsFactors = FALSE
)
synthetic_public_sidecar_file <- file.path(
  package_dir, "validation", "synthetic_identifier_sidecar.tsv.gz"
)
synthetic_public_sidecar_connection <- gzfile(
  synthetic_public_sidecar_file, "wt"
)
utils::write.table(
  synthetic_public_sidecar, synthetic_public_sidecar_connection,
  sep = "\t", quote = FALSE, row.names = FALSE
)
close(synthetic_public_sidecar_connection)
source_metadata_file <- file.path(integration_dir, "metadata_filtered.rds")
source_metadata_object <- readRDS(object_file)
source_metadata_fixture <- data.table::as.data.table(
  source_metadata_object@meta.data
)
rownames(source_metadata_fixture) <- colnames(source_metadata_object)
saveRDS(source_metadata_fixture, source_metadata_file)
anonymize_status <- system2(
  file.path(R.home("bin"), "Rscript"),
  c(
    "sciencedb/anonymize_public_package.R",
    "--repo-dir", repo_dir,
    "--package-dir", package_dir,
    "--source-object", object_file,
    "--source-metadata", source_metadata_file
  )
)
stopifnot(anonymize_status == 0L)
public_metadata <- utils::read.delim(
  gzfile(file.path(package_dir, "metadata", "metadata.tsv.gz")),
  sep = "\t", check.names = FALSE, stringsAsFactors = FALSE
)
public_object <- readRDS(
  file.path(package_dir, "objects", "objects_celltype_plot.rds")
)
public_sidecar <- utils::read.delim(
  gzfile(synthetic_public_sidecar_file),
  sep = "\t", check.names = FALSE, stringsAsFactors = FALSE
)
public_shard_barcodes <- unlist(lapply(
  file.path(
    package_dir, "expression", release_shard_manifest$Relative_Directory,
    "barcodes.tsv.gz"
  ),
  function(path) {
    connection <- gzfile(path, "rt")
    on.exit(close(connection), add = TRUE)
    readLines(connection, warn = FALSE)
  }
), use.names = FALSE)
stopifnot(
  identical(names(public_metadata), expected_release_metadata_columns),
  identical(public_metadata$Cells, sprintf("Cell%07d", seq_len(120L))),
  identical(public_shard_barcodes, public_metadata$Cells),
  identical(public_metadata$Original_Sample_ID, source_original_sample_ids),
  all(grepl("^D[0-9]+$", public_metadata$Donor_ID)),
  all(grepl("^S[0-9]+$", public_metadata$Sample_ID)),
  all(is.na(public_metadata$Library_ID) |
    grepl("^L[0-9]+$", public_metadata$Library_ID)),
  identical(colnames(public_object), public_metadata$Cells),
  identical(names(public_object@meta.data), expected_release_metadata_columns),
  identical(
    public_object@meta.data$Original_Sample_ID,
    source_original_sample_ids
  ),
  identical(public_sidecar$Cells, public_metadata$Cells),
  identical(public_sidecar$Global_Donor_ID, public_metadata$Donor_ID),
  identical(public_sidecar$Specimen_ID, public_metadata$Sample_ID),
  identical(public_sidecar$Library_ID, public_metadata$Library_ID),
  identical(
    public_sidecar$Original_Donor_ID,
    synthetic_public_sidecar$Original_Donor_ID
  )
)

# Refreshing audit sidecars on an already anonymized package must be
# deterministic and must not require the internal identifiers to reappear in
# the primary package tables.
anonymize_again_status <- system2(
  file.path(R.home("bin"), "Rscript"),
  c(
    "sciencedb/anonymize_public_package.R",
    "--repo-dir", repo_dir,
    "--package-dir", package_dir,
    "--source-object", file.path(
      package_dir, "objects", "objects_celltype_plot.rds"
    ),
    "--source-metadata", source_metadata_file
  )
)
stopifnot(anonymize_again_status == 0L)
public_metadata_again <- utils::read.delim(
  gzfile(file.path(package_dir, "metadata", "metadata.tsv.gz")),
  sep = "\t", check.names = FALSE, stringsAsFactors = FALSE
)
public_sidecar_again <- utils::read.delim(
  gzfile(synthetic_public_sidecar_file),
  sep = "\t", check.names = FALSE, stringsAsFactors = FALSE
)
public_manifest_again <- utils::read.delim(
  file.path(package_dir, "provenance", "file_manifest.tsv"),
  sep = "\t", check.names = FALSE, stringsAsFactors = FALSE
)
public_audit_index <- utils::read.delim(
  file.path(package_dir, "provenance", "audit_sidecar_index.tsv"),
  sep = "\t", check.names = FALSE, stringsAsFactors = FALSE
)
synthetic_manifest_row <- public_manifest_again[
  public_manifest_again$path ==
    "validation/synthetic_identifier_sidecar.tsv.gz", ,
  drop = FALSE
]
stopifnot(
  identical(public_metadata_again, public_metadata),
  identical(public_sidecar_again, public_sidecar),
  nrow(synthetic_manifest_row) == 1L,
  synthetic_manifest_row$rows == 120L,
  synthetic_manifest_row$columns == 5L,
  !"provenance/audit_sidecar_index.tsv" %in% public_audit_index$path,
  all(file.exists(file.path(package_dir, public_audit_index$path))),
  identical(
    as.character(public_audit_index$sha256),
    unname(vapply(
      file.path(package_dir, public_audit_index$path),
      processed_file_sha256,
      character(1L)
    ))
  )
)

package_summary_dir <- file.path(test_dir, "publication_tables")
package_summary_status <- system2(
  file.path(R.home("bin"), "Rscript"),
  c(
    "sciencedb/export_package_summary.R",
    "--package-dir", package_dir,
    "--output-dir", package_summary_dir
  )
)
stopifnot(package_summary_status == 0L)
package_summary <- utils::read.delim(
  file.path(package_summary_dir, "sciencedb_package_summary.tsv"),
  sep = "\t", check.names = FALSE, stringsAsFactors = FALSE
)
package_contract <- utils::read.delim(
  file.path(package_summary_dir, "sciencedb_package_contract.tsv"),
  sep = "\t", check.names = FALSE, stringsAsFactors = FALSE
)
package_contract_lookup <- stats::setNames(
  package_contract$Value, package_contract$Field
)
stopifnot(
  sum(package_summary$Payload_Files) == nrow(public_manifest_again),
  as.integer(package_contract_lookup[["Payload_Files"]]) ==
    nrow(public_manifest_again),
  as.integer(package_contract_lookup[["Manifest_and_Checksum_Files"]]) == 2L,
  as.integer(package_contract_lookup[["Total_Files"]]) ==
    nrow(public_manifest_again) + 2L,
  identical(
    package_contract_lookup[["Validation_Status"]],
    "complete post-anonymization payload validation"
  )
)

# The complete neuronal-lineage workflow keeps every broad-lineage cell in
# unsupervised integration and clustering, while the separately audited
# subtype boundary prevents a mixed global cluster from driving marker-based
# subtype names or biological claims.
lineage_integration_path <- file.path(repo_dir, "annotation", "integrate_neuronal_lineage.R")
lineage_integration_source <- if (file.exists(lineage_integration_path)) paste(readLines(lineage_integration_path, warn = FALSE), collapse = "\n") else ""
lineage_marker_path <- file.path(repo_dir, "annotation", "neuronal_lineage_markers.R")
lineage_marker_source <- if (file.exists(lineage_marker_path)) paste(readLines(lineage_marker_path, warn = FALSE), collapse = "\n") else ""
if (nzchar(lineage_integration_source)) stopifnot(
  grepl("excluded_sensitivity", lineage_integration_source, fixed = TRUE),
  grepl(
    "Exclude_From_Subtype_Inference",
    lineage_integration_source,
    fixed = TRUE
  ),
  grepl(
    "Subtype_Inference_Eligible",
    lineage_marker_source,
    fixed = TRUE
  ),
  grepl(
    "lineage_subtype_inference_exclusions.tsv",
    lineage_marker_source,
    fixed = TRUE
  )
)

cat("clean-room synthetic fixture passed\n")
