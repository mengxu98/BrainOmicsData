#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

source("functions/data_paths.R")
source("functions/processed_object.R")
source("functions/integration.R")

integration_dir <- brainomics_data_path("integration_25")
summary_dir <- file.path(integration_dir, "evaluation", "reference_summary")
output_dir <- file.path(integration_dir, "evaluation", "resource_comparison")
summary_file <- file.path(summary_dir, "reference_summary.tsv")
dataset_file <- file.path(summary_dir, "dataset_summary.tsv")
if (any(!file.exists(c(summary_file, dataset_file)))) {
  stop("The frozen reference summary is required before resource comparison")
}
if (dir.exists(output_dir) && length(list.files(output_dir, all.files = FALSE))) {
  stop("Resource comparison refuses to overwrite: ", output_dir)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

summary <- fread(summary_file, sep = "\t")
dataset <- fread(dataset_file, sep = "\t")
summary_value <- function(field) {
  value <- summary[Field == field, Value]
  if (length(value) != 1L) stop("Reference summary field changed: ", field)
  value
}
reference_datasets <- sort(as.character(dataset$Dataset))
if (length(reference_datasets) != 22L || anyDuplicated(reference_datasets)) {
  stop("Resource comparison requires the frozen 22-dataset reference")
}

overlap_definitions <- list(
  `Kim et al. 2024` = list(
    evidence = "exact source-study overlap",
    datasets = c("AllenM1", "EGAS00001006537", "GSE168408")
  ),
  `STAB2` = list(
    evidence = "exact source-study overlap",
    datasets = c("GSE168408", "GSE186538", "GSE97942", "Li_et_al_2018")
  ),
  `ZEBRA` = list(
    evidence = "exact source-study overlap",
    datasets = "GSE168408"
  ),
  `Jorstad et al. 2023` = list(
    evidence = paste(
      "cohort-level relation to the Allen human M1 resource; exact cell-level",
      "overlap is not asserted without a barcode crosswalk"
    ),
    datasets = "AllenM1"
  )
)
overlap <- rbindlist(lapply(names(overlap_definitions), function(resource) {
  definition <- overlap_definitions[[resource]]
  absent <- setdiff(definition$datasets, reference_datasets)
  if (length(absent) > 0L) {
    stop(
      resource, " overlap list contains non-reference datasets: ",
      paste(absent, collapse = ",")
    )
  }
  data.table(
    Comparator = resource,
    Evidence_Level = definition$evidence,
    Reference_Dataset = definition$datasets
  )
}))

overlap_summary <- overlap[, .(
  Source_Overlap_Count = .N,
  Source_Overlap = paste(sort(Reference_Dataset), collapse = ";")
), by = Comparator]

resources <- data.table(
  Resource = c(
    "BrainOmicsData revised reference", "Kim et al. 2024", "STAB2",
    "ZEBRA", "Jorstad et al. 2023"
  ),
  DOI = c(
    NA_character_, "10.1038/s12276-024-01328-6",
    "10.1093/nar/gkad955", "10.1093/nar/gkad990",
    "10.1126/science.adf6812"
  ),
  Source_Studies = c(22L, 8L, 19L, 33L, 1L),
  Human_Cells_or_Nuclei = c(
    summary_value("Cells_or_Nuclei"), "393060", "1504591",
    "2743355", "1155822"
  ),
  Donors = c(
    summary_value("Known_Donors"), "80",
    "not reported as one deduplicated total",
    paste(
      "donor mappings provided; at least 196 controls and 88 Alzheimer",
      "disease donors, with additional conditions"
    ),
    "6 donor IDs in the original 10x analysis release"
  ),
  Age_Coverage = c(
    paste0("5 post-conception weeks to 104 years; ", summary_value("Age_Intervals"),
           " harmonized age intervals"),
    "7 gestational weeks to 90 years",
    "4 post-conception weeks to 95 years in 1239606 public-metadata cells; published total 1504591",
    "exact overall observed range not reliably resolved; public metadata 2253535 cells, 2158103 numeric ages and 95432 blanks; raw 0-93 year codes include zero-coded adult donors and rounded infant ages",
    "29 to 60 years in the original 10x analysis release"
  ),
  Region_Coverage = c(
    paste0(summary_value("Brain_Regions"), " harmonized brain-region labels"),
    "multiple source regions", "63 human brain subregions",
    "39 brain regions across the human and mouse resource", "8 cortical areas"
  ),
  Modalities = c(
    "scRNA-seq and snRNA-seq; RNA arm retained from multiome sources",
    "scRNA-seq and snRNA-seq", "scRNA-seq and snRNA-seq",
    "human and mouse scRNA-seq and snRNA-seq", "snRNA-seq"
  ),
  Annotation_Depth = c(
    paste0(summary_value("Main_Cell_Types"), " adopted major cell classes"),
    "10 major cell types and 22 subtypes", "18 cell types and 71 subtypes",
    "hierarchical global, cortical and non-cortical annotations",
    "24 subclasses and 153 cross-area consensus types"
  ),
  Access_Model = c(
    "versioned derived package, machine-readable sidecars and public workflow",
    "public source datasets and processed atlas", "web database and downloads",
    "web database and Zenodo downloads", "public Allen/BICAN data resources"
  ),
  Validation_Focus = c(
    paste(
      "Raw/scVI/Harmony/RPCA latent-space evaluation, donor-level uncertainty,",
      "source-label concordance, stability and independent EGAD mapping"
    ),
    "scVI integration and developmental lineage analyses",
    "marker-supported cell types across region and developmental stage",
    "hierarchical integration across species, region and disease",
    "donor-integrated neuronal neighborhoods and cross-area consensus taxonomy"
  ),
  Distinct_Reuse_Boundary = c(
    paste(
      "auditable donor-sample-library hierarchy and age-region harmonization;",
      "donor totals reconcile reported source identities across studies without genotype verification;",
      "broad RPCA neighbor mapping with support scores; calibrated rejection remains incomplete"
    ),
    "smaller developmental atlas with scVI covariate modeling",
    "web exploration of dense spatiotemporal coverage",
    "cross-species and neurodegenerative-disease exploration",
    "high-resolution adult cortical taxonomy; numeric row restricted to original 10x release with 6 donors, distinct from the Methods report of 5 newly profiled donors; AllenM1 is a cohort relation, not verified cell-level overlap"
  )
)
resources <- merge(
  resources,
  overlap_summary,
  by.x = "Resource", by.y = "Comparator", all.x = TRUE, sort = FALSE
)
resources[Resource == "BrainOmicsData revised reference", `:=`(
  Source_Overlap_Count = NA_integer_,
  Source_Overlap = NA_character_
)]
resources[Resource == "Jorstad et al. 2023", Source_Overlap_Count := NA_integer_]
setcolorder(resources, c(
  "Resource", "DOI", "Source_Studies", "Human_Cells_or_Nuclei", "Donors",
  "Age_Coverage", "Region_Coverage", "Modalities", "Annotation_Depth",
  "Access_Model", "Validation_Focus", "Source_Overlap_Count",
  "Source_Overlap", "Distinct_Reuse_Boundary"
))

resource_file <- file.path(output_dir, "resource_comparison.tsv")
overlap_file <- file.path(output_dir, "source_overlap.tsv")
input_file <- file.path(output_dir, "input_manifest.tsv")
output_file <- file.path(output_dir, "output_manifest.tsv")
fwrite(resources, resource_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(overlap, overlap_file, sep = "\t", quote = FALSE)
inputs <- c(summary_file, dataset_file)
fwrite(data.table(
  File = inputs,
  Size_Bytes = file.info(inputs)$size,
  SHA256 = vapply(inputs, processed_file_sha256, character(1L))
), input_file, sep = "\t", quote = FALSE)
outputs <- c(resource_file, overlap_file, input_file)
fwrite(data.table(
  File = basename(outputs),
  Size_Bytes = file.info(outputs)$size,
  SHA256 = vapply(outputs, processed_file_sha256, character(1L))
), output_file, sep = "\t", quote = FALSE)
message("[resource-comparison] completed five-resource comparison")
