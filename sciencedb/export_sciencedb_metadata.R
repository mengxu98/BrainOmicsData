#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1]]
}

integration_dir <- normalizePath(
  value_after("--integration-dir", "../../data/BrainOmicsData/integration"),
  mustWork = TRUE
)
out_dir <- value_after(
  "--out-dir",
  "../../data/ScienceDB/HumanBrain_sc_snRNAseq_age_resource"
)
repo_dir <- normalizePath(value_after("--repo-dir", "."), mustWork = FALSE)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  readme = file.path(out_dir, "00_README"),
  source = file.path(out_dir, "01_source_datasets"),
  sample = file.path(out_dir, "02_sample_metadata"),
  cell = file.path(out_dir, "03_cell_metadata"),
  objects = file.path(out_dir, "04_processed_objects"),
  emb = file.path(out_dir, "05_embeddings_and_clusters"),
  qc = file.path(out_dir, "06_quality_control"),
  anno = file.path(out_dir, "07_celltype_annotation"),
  workflow = file.path(out_dir, "08_reusable_workflow"),
  value = file.path(out_dir, "09_value_added_annotations"),
  figtab = file.path(out_dir, "10_figures_and_tables")
)
invisible(lapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE))
dir.create(file.path(paths$workflow, "scripts"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(paths$objects, "by_celltype"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(paths$objects, "by_major_region"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(paths$value, "HAR_related_annotation_optional"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(paths$value, "network_edge_tables_optional"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(paths$figtab, "figure_source_data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(paths$figtab, "main_table_source_files"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(paths$figtab, "supplementary_table_source_files"), recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, file) {
  utils::write.table(
    x, file = file, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE, na = "missing"
  )
}

fill_missing <- function(x, value = "missing") {
  x <- as.character(x)
  x[is.na(x) | x == "" | x == "NA" | x == "Unknown"] <- value
  x
}

metadata_file <- file.path(integration_dir, "metadata_filtered.rds")
if (!file.exists(metadata_file)) {
  stop("Missing metadata file: ", metadata_file)
}
metadata <- readRDS(metadata_file)
metadata[] <- lapply(metadata, as.character)
sanitize_text <- function(x) {
  if (!is.character(x)) return(x)
  x <- gsub(intToUtf8(0x2013), "-", x, fixed = TRUE)
  x <- gsub(intToUtf8(0x2014), "-", x, fixed = TRUE)
  x <- gsub(intToUtf8(0x2212), "-", x, fixed = TRUE)
  x
}
metadata[] <- lapply(metadata, sanitize_text)

add_age_interval_schema <- function(x) {
  age_interval <- c(
    "S1" = "Embryonic",
    "S2" = "Early fetal",
    "S3" = "Early fetal",
    "S4" = "Early mid-fetal",
    "S5" = "Early mid-fetal",
    "S6" = "Late mid-fetal",
    "S7" = "Late fetal",
    "S8" = "Neonatal and early infancy",
    "S9" = "Late infancy",
    "S10" = "Early childhood",
    "S11" = "Middle and late childhood",
    "S12" = "Adolescence",
    "S13" = "Young adulthood",
    "S14" = "Middle adulthood",
    "S15" = "Late adulthood"
  )
  age_range <- c(
    "S1" = "4-8 PCW",
    "S2" = "8-10 PCW",
    "S3" = "10-13 PCW",
    "S4" = "13-16 PCW",
    "S5" = "16-19 PCW",
    "S6" = "19-24 PCW",
    "S7" = "24-38 PCW",
    "S8" = "0-0.5 years",
    "S9" = "0.5-1 years",
    "S10" = "1-6 years",
    "S11" = "6-12 years",
    "S12" = "12-20 years",
    "S13" = "20-40 years",
    "S14" = "40-60 years",
    "S15" = "60+ years"
  )
  if (!"AgeIntervalID" %in% names(x)) x$AgeIntervalID <- x$Stage
  if (!"AgeInterval" %in% names(x)) {
    x$AgeInterval <- unname(age_interval[as.character(x$AgeIntervalID)])
  }
  if (!"AgeRange" %in% names(x)) {
    x$AgeRange <- unname(age_range[as.character(x$AgeIntervalID)])
  }
  x
}
metadata <- add_age_interval_schema(metadata)

required <- c(
  "Cells", "Dataset", "Technology", "Sequence", "Sample", "Sample_ID",
  "CellType_raw", "BrainRegion", "Stage", "DevelopmentStage",
  "AgeIntervalID", "AgeInterval", "AgeRange", "Age", "Sex"
)
missing_required <- setdiff(required, names(metadata))
if (length(missing_required) > 0) {
  stop("metadata_filtered.rds is missing required columns: ", paste(missing_required, collapse = ", "))
}

stage_lookup <- unique(metadata[, c("AgeIntervalID", "AgeInterval", "AgeRange")])
stage_lookup$stage_order <- as.integer(sub("^S", "", stage_lookup$AgeIntervalID))
stage_lookup <- stage_lookup[order(stage_lookup$stage_order), ]
stage_lookup$reported_age_unit <- ifelse(grepl("PCW", stage_lookup$AgeRange), "PCW", "years")
stage_lookup <- stage_lookup[, c("AgeIntervalID", "AgeInterval", "AgeRange", "reported_age_unit")]
names(stage_lookup) <- c("age_interval_id", "age_interval", "age_range", "dominant_reported_age_unit")
write_tsv(stage_lookup, file.path(paths$sample, "age_interval_mapping.tsv"))

source_info <- data.frame(
  source_dataset = c(
    "AllenM1", "EGAD00001006049", "EGAS00001006537", "GSE103723",
    "GSE104276", "GSE144136", "GSE168408", "GSE178175", "GSE186538",
    "GSE199762", "GSE202210", "GSE204683", "GSE207334", "GSE212606",
    "GSE217511", "GSE261983", "GSE296073", "GSE67835", "GSE81475",
    "GSE97942", "HYPOMAP", "Li_et_al_2018", "Ma_et_al_2022",
    "Nowakowski_et_al_2017", "PRJCA015229", "ROSMAP", "SomaMut"
  ),
  source_repository = c(
    "Allen Brain Map", "EGA", "EGA", rep("GEO", 17), "CELLxGENE",
    "publication supplement", "Sestan lab / BrainSCOPE", "publication supplement",
    "CNGBdb", "ROSMAP / AD Knowledge Portal", "project website"
  ),
  source_accession = c(
    "AllenM1", "EGAD00001006049", "EGAS00001006537", "GSE103723",
    "GSE104276", "GSE144136", "GSE168408", "GSE178175", "GSE186538",
    "GSE199762; phs003509.v1.p1", "GSE202210", "GSE204683",
    "GSE207334", "GSE212606", "GSE217511", "GSE261983",
    "GSE296073", "GSE67835", "GSE81475", "GSE97942", "HYPOMAP",
    "Li et al. 2018", "Ma et al. 2022", "Nowakowski et al. 2017",
    "PRJCA015229", "ROSMAP", "SomaMut"
  ),
  publication_doi = c(
    NA, "10.1126/science.adf1226", NA, "10.1126/sciadv.adg3754",
    "10.1126/sciadv.adg3754", NA, NA, NA, "10.1126/sciadv.adg3754",
    "10.1038/s41586-023-06981-x", NA, "10.1126/sciadv.adg3754",
    "10.1126/science.abo7257", "10.1126/sciadv.adg3754",
    "10.1038/s41467-022-34975-2", "10.1126/science.adi5199",
    "10.1038/s41586-025-09362-8", "10.1073/pnas.1507125112",
    "10.1016/j.celrep.2016.08.038", "10.1038/nbt.4038",
    "10.1038/s41586-024-08504-8", "10.1126/science.aat7615",
    "10.1126/science.abo7257", "10.1126/science.aap8809",
    "10.1016/j.xgen.2024.100703", "10.1016/j.cell.2023.08.039",
    "10.1038/s41586-025-09435-8"
  ),
  access_level = c(
    "public", "controlled", "controlled", "public", "public", "public",
    "public", "public", "public", "mixed_public_controlled", "public",
    "public", "public", "public", "public", "public", "public", "public",
    "public", "public", "public", "public", "public", "public", "public",
    "controlled_or_restricted", "public"
  ),
  science_db_redistribution_plan = c(
    rep("derived_metadata_and_processed_outputs", 27)
  ),
  citation_status = c(
    "needs_source_publication_confirmation",
    "doi_from_processing_script",
    "needs_source_publication_confirmation",
    "doi_from_processing_script",
    "doi_from_processing_script",
    "needs_source_publication_confirmation",
    "needs_source_publication_confirmation",
    "needs_source_publication_confirmation",
    "doi_from_processing_script",
    "doi_from_processing_script",
    "needs_source_publication_confirmation",
    "doi_from_processing_script",
    "doi_from_processing_script",
    "doi_from_processing_script",
    "doi_from_processing_script",
    "doi_from_processing_script",
    "doi_from_processing_script",
    "doi_from_processing_script",
    "doi_from_processing_script",
    "doi_from_processing_script",
    "doi_from_processing_script",
    "doi_from_processing_script",
    "doi_from_processing_script",
    "doi_from_processing_script",
    "doi_from_processing_script",
    "doi_from_processing_script",
    "doi_from_processing_script"
  ),
  stringsAsFactors = FALSE
)

dataset_counts <- as.data.frame(table(metadata$Dataset), stringsAsFactors = FALSE)
names(dataset_counts) <- c("source_dataset", "number_of_cells_after_filtering")
sample_counts <- aggregate(
  Sample_ID ~ Dataset, metadata,
  function(x) length(unique(x))
)
names(sample_counts) <- c("source_dataset", "number_of_samples_or_donors")
region_counts <- aggregate(
  BrainRegion ~ Dataset, metadata,
  function(x) length(unique(x))
)
names(region_counts) <- c("source_dataset", "number_of_brain_regions")
stage_counts_by_dataset <- aggregate(
  AgeIntervalID ~ Dataset, metadata,
  function(x) length(unique(x))
)
names(stage_counts_by_dataset) <- c("source_dataset", "number_of_age_intervals")
source_summary <- Reduce(
  function(x, y) merge(x, y, by = "source_dataset", all.x = TRUE),
  list(source_info, dataset_counts, sample_counts, region_counts, stage_counts_by_dataset)
)
source_summary <- source_summary[order(source_summary$source_dataset), ]
write_tsv(source_summary, file.path(paths$source, "source_dataset_summary.tsv"))
write_tsv(
  source_summary[, c(
    "source_dataset", "source_repository", "source_accession",
    "publication_doi", "access_level", "science_db_redistribution_plan"
  )],
  file.path(paths$source, "source_accessions.tsv")
)
write_tsv(
  source_summary[, c("source_dataset", "publication_doi", "citation_status")],
  file.path(paths$source, "publication_references.tsv")
)
write_tsv(
  source_summary[, c(
    "source_dataset", "source_repository", "source_accession", "access_level",
    "science_db_redistribution_plan"
  )],
  file.path(paths$source, "access_level_public_controlled.tsv")
)

processing_summary <- data.frame(
  source_dataset = source_summary$source_dataset,
  processing_script = paste0("processing/", source_summary$source_dataset, ".R"),
  retained_cells_or_nuclei = source_summary$number_of_cells_after_filtering,
  sample_or_donor_groups = source_summary$number_of_samples_or_donors,
  brain_region_count = source_summary$number_of_brain_regions,
  age_interval_count = source_summary$number_of_age_intervals,
  celltype_source_field = "CellType_raw",
  special_handling = ifelse(
    source_summary$source_dataset == "GSE296073",
    "organoid-derived records excluded; only real human postmortem embryonic/perinatal samples retained",
    "standard BrainOmicsData processing"
  ),
  stringsAsFactors = FALSE
)
write_tsv(processing_summary, file.path(paths$source, "source_dataset_processing.tsv"))

organoid_check <- data.frame(
  source_dataset = source_info$source_dataset,
  organoid_excluded = ifelse(source_info$source_dataset == "GSE296073", "yes_partial", "not_applicable_or_not_detected"),
  notes = ifelse(
    source_info$source_dataset == "GSE296073",
    "Processing script uses postmortem embryonic and perinatal human samples and excludes organoid-derived records from the integrated resource.",
    "No organoid exclusion rule is currently recorded for this dataset."
  ),
  stringsAsFactors = FALSE
)
write_tsv(organoid_check, file.path(paths$source, "organoid_exclusion_check.tsv"))

sample_metadata <- aggregate(
  Cells ~ Dataset + Sample + Sample_ID + Technology + Sequence + BrainRegion + Stage +
    DevelopmentStage + AgeIntervalID + AgeInterval + AgeRange + Age + Sex,
  metadata,
  length
)
names(sample_metadata) <- c(
  "source_dataset", "sample_original", "sample_id", "sequencing_platform",
  "sequencing_modality", "brain_region_harmonized", "legacy_stage",
  "legacy_development_stage", "age_interval_id", "age_interval", "age_range",
  "reported_age", "sex", "number_of_cells_after_filtering"
)
sample_metadata$source_accession <- source_summary$source_accession[
  match(sample_metadata$source_dataset, source_summary$source_dataset)
]
sample_metadata$source_repository <- source_summary$source_repository[
  match(sample_metadata$source_dataset, source_summary$source_dataset)
]
sample_metadata$publication_doi <- source_summary$publication_doi[
  match(sample_metadata$source_dataset, source_summary$source_dataset)
]
sample_metadata$access_level <- source_summary$access_level[
  match(sample_metadata$source_dataset, source_summary$source_dataset)
]
sample_metadata$species <- "Homo sapiens"
sample_metadata$organism_part <- "brain"
sample_metadata$donor_id <- sample_metadata$sample_id
sample_metadata$publication_reference <- sample_metadata$publication_doi
sample_metadata$brain_region_original <- "not reported"
sample_metadata$reported_age_unit <- ifelse(grepl("PCW", sample_metadata$reported_age), "PCW", "years")
sample_metadata$library_protocol <- "not reported"
sample_metadata$cell_number_before_filtering <- "missing"
sample_metadata$cell_number_after_filtering <- sample_metadata$number_of_cells_after_filtering
sample_metadata$single_cell_or_single_nucleus <- ifelse(
  grepl("nucleus|snRNA|snATAC", sample_metadata$sequencing_modality, ignore.case = TRUE),
  "single-nucleus", "single-cell"
)
sample_metadata$inclusion_status <- "included"
sample_metadata$exclusion_reason <- "not applicable"
sample_metadata[] <- lapply(sample_metadata, fill_missing)
sample_metadata <- sample_metadata[, c(
  "sample_id", "donor_id", "source_dataset", "source_accession",
  "source_repository", "publication_reference", "access_level",
  "species", "organism_part", "brain_region_original",
  "brain_region_harmonized", "reported_age", "reported_age_unit",
  "age_interval_id", "age_interval", "age_range", "sex", "sequencing_modality",
  "single_cell_or_single_nucleus", "sequencing_platform", "library_protocol",
  "cell_number_before_filtering", "cell_number_after_filtering",
  "number_of_cells_after_filtering", "inclusion_status", "exclusion_reason"
)]
sample_metadata <- sample_metadata[order(sample_metadata$source_dataset, sample_metadata$sample_id), ]
write_tsv(sample_metadata, file.path(paths$sample, "sample_metadata_harmonized.tsv"))

region_map <- data.frame(
  brain_region_harmonized = sort(unique(metadata$BrainRegion)),
  region_level_1 = NA_character_,
  region_level_2 = NA_character_,
  stringsAsFactors = FALSE
)
write_tsv(region_map, file.path(paths$sample, "brain_region_mapping.tsv"))

platform_map <- unique(metadata[, c("Technology", "Sequence")])
names(platform_map) <- c("sequencing_platform", "sequencing_modality")
platform_map$library_protocol <- NA_character_
platform_map <- platform_map[order(platform_map$sequencing_platform, platform_map$sequencing_modality), ]
write_tsv(platform_map, file.path(paths$sample, "sequencing_platform_mapping.tsv"))

cell_metadata <- data.frame(
  cell_id = metadata$Cells,
  sample_id = metadata$Sample_ID,
  donor_id = metadata$Sample_ID,
  source_dataset = metadata$Dataset,
  source_cell_id = metadata$Cells,
  brain_region_harmonized = metadata$BrainRegion,
  age_interval_id = metadata$AgeIntervalID,
  age_interval = metadata$AgeInterval,
  age_range = metadata$AgeRange,
  reported_age = metadata$Age,
  sex = metadata$Sex,
  sequencing_modality = metadata$Sequence,
  sequencing_platform = metadata$Technology,
  cell_type_raw = metadata$CellType_raw,
  cluster_id = NA_character_,
  cell_type = NA_character_,
  n_counts = "missing",
  n_genes = "missing",
  percent_mito = "missing",
  percent_ribo = "missing",
  quality_control_status = "included_after_metadata_filtering",
  stringsAsFactors = FALSE
)
cell_metadata[] <- lapply(cell_metadata, fill_missing)
gz <- gzfile(file.path(paths$cell, "cell_metadata_harmonized.tsv.gz"), "wt")
write_tsv(cell_metadata, gz)
close(gz)

cell_counts_by_dataset <- source_summary[, c(
  "source_dataset", "number_of_cells_after_filtering",
  "number_of_samples_or_donors", "number_of_brain_regions", "number_of_age_intervals"
)]
write_tsv(cell_counts_by_dataset, file.path(paths$cell, "cell_counts_by_dataset.tsv"))

cell_counts_by_age_interval <- as.data.frame(
  xtabs(~ AgeIntervalID, metadata),
  stringsAsFactors = FALSE
)
names(cell_counts_by_age_interval) <- c("age_interval_id", "cell_count")
write_tsv(cell_counts_by_age_interval, file.path(paths$cell, "cell_counts_by_age_interval.tsv"))

cell_counts_by_brain_region <- as.data.frame(
  xtabs(~ BrainRegion, metadata),
  stringsAsFactors = FALSE
)
names(cell_counts_by_brain_region) <- c("brain_region_harmonized", "cell_count")
write_tsv(cell_counts_by_brain_region, file.path(paths$cell, "cell_counts_by_brain_region.tsv"))

cell_counts_by_cell_type <- as.data.frame(
  xtabs(~ CellType_raw, metadata),
  stringsAsFactors = FALSE
)
names(cell_counts_by_cell_type) <- c("cell_type", "cell_count")
cell_counts_by_cell_type <- cell_counts_by_cell_type[cell_counts_by_cell_type$cell_count > 0, ]
write_tsv(cell_counts_by_cell_type, file.path(paths$cell, "cell_counts_by_cell_type.tsv"))
write_tsv(cell_counts_by_cell_type, file.path(paths$anno, "cell_type_counts.tsv"))

cell_counts_by_age_region_celltype <- as.data.frame(
  xtabs(~ AgeIntervalID + BrainRegion + CellType_raw, metadata),
  stringsAsFactors = FALSE
)
names(cell_counts_by_age_region_celltype) <- c(
  "age_interval_id", "brain_region_harmonized", "cell_type_raw", "cell_count"
)
cell_counts_by_age_region_celltype <- cell_counts_by_age_region_celltype[
  cell_counts_by_age_region_celltype$cell_count > 0,
]
write_tsv(
  cell_counts_by_age_region_celltype,
  file.path(paths$cell, "cell_counts_by_age_region_celltype.tsv")
)

cell_counts_by_sex_technology <- as.data.frame(
  xtabs(~ Sex + Technology + Sequence, metadata),
  stringsAsFactors = FALSE
)
names(cell_counts_by_sex_technology) <- c(
  "sex", "sequencing_platform", "sequencing_modality", "cell_count"
)
cell_counts_by_sex_technology <- cell_counts_by_sex_technology[
  cell_counts_by_sex_technology$cell_count > 0,
]
write_tsv(cell_counts_by_sex_technology, file.path(paths$cell, "cell_counts_by_sex_technology.tsv"))

field_names <- c(
  "reported_age", "sex", "brain_region_harmonized", "donor_id", "sample_id",
  "sequencing_modality", "sequencing_platform", "library_protocol",
  "source_accession", "access_level", "publication_reference"
)
missingness <- do.call(rbind, lapply(split(sample_metadata, sample_metadata$source_dataset), function(df) {
  data.frame(
    source_dataset = df$source_dataset[[1]],
    field = field_names,
    complete_fraction = vapply(field_names, function(field) {
      vals <- df[[field]]
      mean(!(is.na(vals) | vals == "" | vals == "NA" | vals == "Unknown" | vals == "missing" | vals == "not reported"))
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
}))
missingness$status <- ifelse(
  missingness$complete_fraction == 1, "complete",
  ifelse(missingness$complete_fraction == 0, "missing", "partial")
)
write_tsv(missingness, file.path(paths$sample, "metadata_missingness_summary.tsv"))
write_tsv(missingness, file.path(paths$qc, "metadata_missingness_summary.tsv"))

qc_by_sample <- sample_metadata[, c(
  "sample_id", "source_dataset", "cell_number_before_filtering", "cell_number_after_filtering",
  "brain_region_harmonized", "age_interval_id", "age_interval", "age_range",
  "sex", "sequencing_modality",
  "sequencing_platform"
)]
write_tsv(qc_by_sample, file.path(paths$qc, "qc_metrics_by_sample.tsv"))
write_tsv(cell_counts_by_dataset, file.path(paths$qc, "qc_metrics_by_dataset.tsv"))

removed_gene_patterns <- data.frame(
  category = c(
    "ERCC spike-ins", "ribosomal protein genes", "mitochondrial genes",
    "uncharacterized LOC/LINC/RP/AC/AL/AP/CT/CH/FAM genes", "MALAT1",
    "hemoglobin genes", "open reading frame shorthand"
  ),
  pattern = c(
    "^ERCC", "^RPLP|^RPSL", "^MT-|^mt-",
    "^LOC|^LINC|^RP[0-9]|^AC[0-9]|^AL[0-9]|^AP[0-9]|^CT[0-9]|^CH[0-9]|^FAM[0-9]",
    "MALAT1", "^HB[^(P)]", "orf"
  ),
  source_script = "HARNexus/code/datasets/datasets_integration_02.R",
  removed_gene_count = "computed during full integration rerun; see objects_filtered.rds feature count",
  stringsAsFactors = FALSE
)
write_tsv(removed_gene_patterns, file.path(paths$qc, "removed_gene_categories.tsv"))

celltype_cluster_mapping <- data.frame(
  cluster_id = NA_character_,
  cell_type = NA_character_,
  note = "Final cluster-to-cell-type mapping is exported by export_sciencedb_integrated_object.R after loading objects_celltypes.rds.",
  stringsAsFactors = FALSE
)
write_tsv(celltype_cluster_mapping, file.path(paths$cell, "celltype_cluster_mapping.tsv"))
write_tsv(celltype_cluster_mapping, file.path(paths$anno, "celltype_cluster_mapping.tsv"))

marker_genes <- data.frame(
  cell_type = c(
    rep("Radial glia", 3), rep("Endothelial cells", 4),
    rep("Inhibitory neurons", 3), rep("Oligodendrocyte progenitor cells", 5),
    rep("Microglia", 3), "Neuroblasts", rep("Excitatory neurons", 3),
    rep("Astrocytes", 5), rep("Oligodendrocytes", 3)
  ),
  marker_gene = c(
    "PAX6", "VIM", "GLI3", "CLDN5", "PECAM1", "VWF", "FLT1",
    "GAD1", "GAD2", "SLC6A1", "PDGFRA", "CSPG4", "OLIG1", "OLIG2",
    "SOX10", "CX3CR1", "P2RY12", "CSF1R", "STMN2", "SLC17A7",
    "CAMK2A", "SATB2", "GFAP", "AQP4", "ALDH1L1", "FGFR3", "GJA1",
    "MOG", "MAG", "CLDN11"
  ),
  validation_use = "marker-based manual cell-type annotation",
  stringsAsFactors = FALSE
)
write_tsv(marker_genes, file.path(paths$anno, "marker_gene_list.tsv"))
write_tsv(
  data.frame(
    rule_id = seq_len(nrow(marker_genes)),
    cell_type = marker_genes$cell_type,
    marker_gene = marker_genes$marker_gene,
    expected_pattern = "enriched expression or accessibility in the annotated cell type",
    stringsAsFactors = FALSE
  ),
  file.path(paths$anno, "celltype_annotation_rules.tsv")
)
write_tsv(
  data.frame(
    cell_type = unique(marker_genes$cell_type),
    marker_gene_count = as.integer(table(marker_genes$cell_type)[unique(marker_genes$cell_type)]),
    marker_expression_summary_status = "computed by export_sciencedb_integrated_object.R when the integrated Seurat object is loaded",
    stringsAsFactors = FALSE
  ),
  file.path(paths$anno, "marker_expression_summary.tsv")
)
write_tsv(
  data.frame(
    cluster_id = "pending_heavy_export",
    cell_type = "pending_heavy_export",
    marker_support = "cluster-level marker validation is written by export_sciencedb_integrated_object.R",
    stringsAsFactors = FALSE
  ),
  file.path(paths$anno, "cluster_marker_validation.tsv")
)

controlled <- source_summary[source_summary$access_level != "public", c(
  "source_dataset", "source_repository", "source_accession",
  "publication_doi", "access_level", "science_db_redistribution_plan"
)]
controlled$access_url_or_database <- controlled$source_repository
controlled$original_publication <- controlled$publication_doi
controlled$access_instructions <- "Request raw controlled-access human data from the original repository or data access committee under the original data-use terms."
controlled$redistribution_status <- "restricted raw data are not redistributed in ScienceDB; only permitted harmonized metadata and derived outputs are staged."
write_tsv(controlled, file.path(paths$source, "controlled_access_handling.tsv"))
writeLines(
  c(
    "Controlled-access data handling",
    "",
    "Restricted or controlled-access raw human data are not redistributed in this ScienceDB package.",
    "The package records original accessions, repositories, publications and access instructions.",
    "Users must request raw controlled-access data from the original repository or data access committee and follow the original data-use conditions.",
    "",
    paste(capture.output(write_tsv(controlled, stdout())), collapse = "\n")
  ),
  file.path(paths$readme, "controlled_access_notes.txt")
)

workflow_config <- c(
  "resource_name: Human brain single-cell and single-nucleus transcriptomic resource across age intervals",
  paste0("integration_dir: ", integration_dir),
  paste0("source_code_repository_local_path: ", repo_dir),
  "normalization: Seurat NormalizeData",
  "highly_variable_genes: 3000",
  "principal_components: 50",
  "integration: Seurat v5 RPCAIntegration",
  "clustering: FindNeighbors and FindClusters resolution 1",
  "celltype_annotation: manual marker-based annotation"
)
writeLines(workflow_config, file.path(paths$workflow, "workflow_config.yaml"))
writeLines(capture.output(sessionInfo()), file.path(paths$workflow, "software_versions.txt"))
writeLines(
  c(
    "name: BrainOmicsData-sciencedb",
    "channels:",
    "  - conda-forge",
    "  - bioconda",
    "dependencies:",
    "  - r-base",
    "  - r-seurat",
    "  - r-ggplot2",
    "  - r-data.table",
    "  - python"
  ),
  file.path(paths$workflow, "conda_environment.yml")
)
writeLines(
  c(
    "# Example commands",
    "",
    "Run from the BrainOmicsData repository root:",
    "",
    "```sh",
    "bash 01_datasets_download.sh",
    "bash 02_datasets_preprocessing.sh",
    "bash 03_datasets_integration.sh",
    "bash 04_datasets_plotting.sh",
    "bash 05_sciencedb.sh",
    "```"
  ),
  file.path(paths$workflow, "example_commands.md")
)

readme <- c(
  "# Human brain single-cell and single-nucleus transcriptomic resource across age intervals",
  "",
  "This package describes an integrated human brain sc/snRNA-seq resource prepared for a Neuroscience Bulletin Data Paper submission.",
  "",
  "The current export contains harmonized source-dataset, sample-level and cell-level metadata derived from the integrated BrainOmicsData/HARNexus workflow. Controlled-access source datasets are cited by their original accessions and are not redistributed as restricted raw data.",
  "",
  "Core folders:",
  "- 01_source_datasets: source accessions, publication DOIs, access levels and redistribution notes.",
  "- 02_sample_metadata: harmonized sample metadata, age interval mapping, region mapping and missingness summary.",
  "- 03_cell_metadata: harmonized cell metadata and count summaries.",
  "- 04_processed_objects: integrated h5ad/h5Seurat/RDS outputs when generated by the heavy export step.",
  "- 05_embeddings_and_clusters: PCA/UMAP coordinates and cluster assignments when generated by the heavy export step.",
  "- 06_quality_control: QC, missingness, integration and filtering summaries.",
  "- 07_celltype_annotation: marker genes and annotation rules.",
  "- 08_reusable_workflow: scripts, configuration and software versions.",
  "- 10_figures_and_tables: source data for Data Paper figures and tables.",
  "",
  "Missing values are encoded as NA. Files ending in .tsv.gz are gzip-compressed tab-separated tables."
)
writeLines(readme, file.path(paths$readme, "README.md"))
writeLines(
  c(
    "Reuse notes",
    "",
    "The ScienceDB record should assign final DOI and CSTR identifiers before manuscript submission.",
    "Restricted or controlled-access raw human data should be requested from the original repositories or data access committees.",
    "Derived metadata and processed outputs should not be used to attempt donor re-identification."
  ),
  file.path(paths$readme, "license_and_reuse_notes.txt")
)

dictionary <- data.frame(
  file_name = c(
    "01_source_datasets/source_dataset_summary.tsv",
    "02_sample_metadata/sample_metadata_harmonized.tsv",
    "03_cell_metadata/cell_metadata_harmonized.tsv.gz",
    "03_cell_metadata/cell_counts_by_age_region_celltype.tsv",
    "05_embeddings_and_clusters/umap_coordinates.tsv.gz",
    "05_embeddings_and_clusters/pca_coordinates.tsv.gz",
    "05_embeddings_and_clusters/cluster_assignments.tsv.gz",
    "06_quality_control/qc_metrics_by_cell.tsv.gz",
    "06_quality_control/metadata_missingness_summary.tsv",
    "07_celltype_annotation/marker_gene_list.tsv"
  ),
  column_name = c(
    "source_dataset; source_accession; publication_doi; access_level",
    "sample_id; brain_region_harmonized; age_interval_id; age_interval; age_range; sex; sequencing_modality",
    "cell_id; sample_id; source_dataset; brain_region_harmonized; age_interval_id; age_interval; age_range; cell_type_raw",
    "age_interval_id; brain_region_harmonized; cell_type_raw; cell_count",
    "cell_id; umap_1; umap_2",
    "cell_id; PC_1 ... PC_n",
    "cell_id; source_dataset; cluster_id; cell_type",
    "cell_id; source_dataset; n_counts; n_genes; percent_mito",
    "source_dataset; field; complete_fraction; status",
    "cell_type; marker_gene; validation_use"
  ),
  description = c(
    "Dataset provenance, repository identifiers, DOI status and access route.",
    "One row per harmonized dataset-sample-region-age-sex-technology combination.",
    "One row per retained cell or nucleus after metadata filtering.",
    "Sparse coverage counts for age-by-region-by-raw-cell-type combinations.",
    "UMAP coordinates exported from the Seurat plotting object.",
    "PCA coordinates exported from the Seurat plotting object.",
    "Cluster and major cell-type assignments exported from the Seurat plotting object.",
    "Cell-level QC metrics exported from the Seurat plotting object.",
    "Dataset-level completeness assessment for key metadata fields.",
    "Representative markers used for manual cell-type annotation validation."
  ),
  file = c(
    "01_source_datasets/source_dataset_summary.tsv",
    "02_sample_metadata/sample_metadata_harmonized.tsv",
    "03_cell_metadata/cell_metadata_harmonized.tsv.gz",
    "03_cell_metadata/cell_counts_by_age_region_celltype.tsv",
    "05_embeddings_and_clusters/umap_coordinates.tsv.gz",
    "05_embeddings_and_clusters/pca_coordinates.tsv.gz",
    "05_embeddings_and_clusters/cluster_assignments.tsv.gz",
    "06_quality_control/qc_metrics_by_cell.tsv.gz",
    "06_quality_control/metadata_missingness_summary.tsv",
    "07_celltype_annotation/marker_gene_list.tsv"
  ),
  field = c(
    "source_dataset; source_accession; publication_doi; access_level",
    "sample_id; brain_region_harmonized; age_interval_id; age_interval; age_range; sex; sequencing_modality",
    "cell_id; sample_id; source_dataset; brain_region_harmonized; age_interval_id; age_interval; age_range; cell_type_raw",
    "age_interval_id; brain_region_harmonized; cell_type_raw; cell_count",
    "cell_id; umap_1; umap_2",
    "cell_id; PC_1 ... PC_n",
    "cell_id; source_dataset; cluster_id; cell_type",
    "cell_id; source_dataset; n_counts; n_genes; percent_mito",
    "source_dataset; field; complete_fraction; status",
    "cell_type; marker_gene; validation_use"
  ),
  definition = c(
    "Dataset provenance, repository identifiers, DOI status and access route.",
    "One row per harmonized dataset-sample-region-age-sex-technology combination.",
    "One row per retained cell or nucleus after metadata filtering.",
    "Sparse coverage counts for age-by-region-by-raw-cell-type combinations.",
    "UMAP coordinates exported from the Seurat plotting object.",
    "PCA coordinates exported from the Seurat plotting object.",
    "Cluster and major cell-type assignments exported from the Seurat plotting object.",
    "Cell-level QC metrics exported from the Seurat plotting object.",
    "Dataset-level completeness assessment for key metadata fields.",
    "Representative markers used for manual cell-type annotation validation."
  ),
  stringsAsFactors = FALSE
)
write_tsv(dictionary, file.path(paths$readme, "data_dictionary.tsv"))

message("Metadata ScienceDB export complete: ", out_dir)
