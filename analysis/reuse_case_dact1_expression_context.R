#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SeuratObject)
})

source("functions/data_paths.R")
source("functions/processed_object.R")

source("functions/utils.R")
data.table::setDTthreads(2L)
integration_dir <- brainomics_data_path("integration_25")
output_dir <- file.path(
  integration_dir, "evaluation", "reuse_case_dact1_s15_pfc"
)
if (dir.exists(output_dir)) {
  stop("Reuse-case output already exists: ", output_dir)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(
  c(
    paste0("PID=", Sys.getpid()),
    paste0("Started_UTC=", format(Sys.time(), tz = "UTC", usetz = TRUE))
  ),
  file.path(output_dir, "_INCOMPLETE")
)

seed <- 20260730L
bootstrap_replicates <- 10000L
minimum_cells <- 20L
minimum_total_counts <- 1000
target_gene <- "DACT1"
target_age <- "S15"
target_region <- "Prefrontal cortex"
target_celltype <- "Oligodendrocytes"

assignment_file <- file.path(integration_dir, "annotation", "celltype_assignments.rds")
assignment <- as.data.table(read_celltype_assignments())
metadata <- readRDS(file.path(integration_dir, "metadata_filtered.rds"))
setDT(metadata)
fields <- c("Cells", "Dataset", "Global_Donor_ID", "Sample_ID", "AgeIntervalID", "BrainRegion")
stopifnot(all(fields %in% names(metadata)))
metadata <- metadata[, ..fields]
metadata <- metadata[match(assignment$Cells, Cells)]
stopifnot(identical(metadata$Cells, assignment$Cells), !anyDuplicated(metadata$Cells))
gc()
metadata[, Main_CellType := assignment$CellType]
required_metadata <- c(
  "Cells", "Dataset", "Global_Donor_ID", "Sample_ID",
  "AgeIntervalID", "BrainRegion", "Main_CellType"
)
if (!all(required_metadata %in% names(metadata))) {
  stop("Reuse-case metadata fields are incomplete")
}
cohort <- metadata[
  AgeIntervalID == target_age & BrainRegion == target_region &
    !is.na(Global_Donor_ID) & nzchar(Global_Donor_ID),
  ..required_metadata
]
if (nrow(cohort) == 0L || anyDuplicated(cohort$Cells)) {
  stop("Prespecified S15 prefrontal-cortex cohort is empty or duplicated")
}

# Read only the source datasets contributing to the prespecified cohort.
cell_level <- list()
input_files <- character()
for (dataset in sort(unique(cohort$Dataset))) {
  dataset_dir <- file.path("../../data/BrainOmicsData/processed", dataset)
  format <- fread(file.path(dataset_dir, "processed_object_format.tsv"))
  object_file <- file.path(dataset_dir, format$Full_Object[[1L]])
  object <- load_processed_object(object_file)
  canonical <- fread(file.path(dataset_dir, format$Canonical_Metadata[[1L]]),
    select = c("Cells", "Original_Cell_ID")
  )
  original <- canonical$Original_Cell_ID[match(colnames(object), canonical$Cells)]
  stopifnot(!anyNA(original), !anyDuplicated(original))
  qualified <- setNames(paste(dataset, original, sep = "::"), colnames(object))
  input_files <- c(input_files, object_file)
  count_layers <- SeuratObject::Layers(object, assay = "RNA", search = "^counts")
  for (layer in count_layers) {
    layer_cells <- SeuratObject::Cells(object[["RNA"]], layer = layer)
    selected <- layer_cells[qualified[layer_cells] %in% cohort$Cells]
    if (!length(selected)) next
    counts <- SeuratObject::LayerData(object, assay = "RNA", layer = layer, cells = selected)
    if (!target_gene %in% rownames(counts)) stop(dataset, " did not measure ", target_gene)
    context <- cohort[match(unname(qualified[colnames(counts)]), Cells)]
    stopifnot(!anyNA(context$Cells))
    cell_level[[paste(dataset, layer)]] <- data.table(
      Dataset = context$Dataset,
      Donor_ID = context$Global_Donor_ID, Sample_ID = context$Sample_ID,
      CellType = context$Main_CellType, Gene_Count = as.numeric(counts[target_gene, ]),
      Total_Count = as.numeric(Matrix::colSums(counts))
    )
    rm(counts)
    gc()
  }
  rm(object, canonical)
  gc()
}
cell_level <- rbindlist(cell_level, use.names = TRUE)
if (nrow(cell_level) != nrow(cohort)) {
  stop("Reuse-case count layers do not cover every prespecified cohort cell")
}

donor_pseudobulk <- cell_level[, .(
  Cells = .N,
  Samples = uniqueN(Sample_ID),
  Gene_Count = sum(Gene_Count),
  Total_Count = sum(Total_Count)
), by = .(Dataset, Donor_ID, CellType)]
donor_pseudobulk[, Included :=
  Cells >= minimum_cells & Total_Count >= minimum_total_counts]
donor_pseudobulk[, Gene_CPM :=
  fifelse(Total_Count > 0, 1e6 * Gene_Count / Total_Count, NA_real_)]
donor_pseudobulk[, Log1p_Gene_CPM := log1p(Gene_CPM)]

target_all <- donor_pseudobulk[CellType == target_celltype]
target_audit <- target_all[Included == TRUE]
if (sum(target_all$Cells) != sum(cohort$Main_CellType == target_celltype) ||
  uniqueN(target_all$Donor_ID) < 2L ||
  !setequal(target_all$Dataset, cohort[Main_CellType == target_celltype, Dataset])) {
  stop(
    "Prespecified unfiltered target coverage changed: cells=",
    sum(target_all$Cells), "; donors=", uniqueN(target_all$Donor_ID),
    "; datasets=", uniqueN(target_all$Dataset)
  )
}
if (sum(target_audit$Cells) > sum(target_all$Cells) ||
  uniqueN(target_audit$Donor_ID) < 2L ||
  any(!target_audit$Dataset %in% target_all$Dataset)) {
  stop(
    "Prespecified eligible target coverage changed: cells=",
    sum(target_audit$Cells),
    "; donors=", uniqueN(target_audit$Donor_ID),
    "; datasets=", uniqueN(target_audit$Dataset)
  )
}

context_summary <- donor_pseudobulk[Included == TRUE, .(
  Cells = sum(Cells),
  Donors = uniqueN(Donor_ID),
  Datasets = uniqueN(Dataset),
  Median_Log1p_CPM = stats::median(Log1p_Gene_CPM),
  Mean_Log1p_CPM = mean(Log1p_Gene_CPM),
  Donor_Fraction_Detected = mean(Gene_Count > 0)
), by = CellType][order(-Mean_Log1p_CPM)]

target_coverage <- data.table(
  Cohort = c("all target donor-cell-type groups", "analysis-eligible target"),
  Cells = c(sum(target_all$Cells), sum(target_audit$Cells)),
  Donors = c(uniqueN(target_all$Donor_ID), uniqueN(target_audit$Donor_ID)),
  Datasets = c(uniqueN(target_all$Dataset), uniqueN(target_audit$Dataset)),
  Inclusion_Rule = c(
    "S15; Prefrontal cortex; Oligodendrocytes; known donor",
    paste0(
      "preceding cohort plus Cells >= ", minimum_cells,
      " and Total_Count >= ", minimum_total_counts
    )
  )
)

bootstrap_dataset_equal <- function(paired, replicates, seed) {
  set.seed(seed)
  datasets <- sort(unique(paired$Dataset))
  estimates <- numeric(replicates)
  for (iteration in seq_len(replicates)) {
    means <- vapply(datasets, function(dataset) {
      values <- paired[Dataset == dataset, Difference]
      mean(sample(values, length(values), replace = TRUE))
    }, numeric(1L))
    estimates[[iteration]] <- mean(means)
  }
  stats::quantile(estimates, c(0.025, 0.975), names = FALSE)
}

comparators <- setdiff(
  sort(unique(donor_pseudobulk$CellType[donor_pseudobulk$Included])),
  target_celltype
)
paired_effects <- rbindlist(lapply(seq_along(comparators), function(index) {
  comparator <- comparators[[index]]
  target <- donor_pseudobulk[
    CellType == target_celltype & Included == TRUE,
    .(Dataset, Donor_ID, Target = Log1p_Gene_CPM)
  ]
  other <- donor_pseudobulk[
    CellType == comparator & Included == TRUE,
    .(Dataset, Donor_ID, Comparator = Log1p_Gene_CPM)
  ]
  paired <- merge(target, other, by = c("Dataset", "Donor_ID"))
  paired[, Difference := Target - Comparator]
  if (nrow(paired) < 3L || uniqueN(paired$Dataset) < 2L) {
    return(NULL)
  }
  dataset_means <- paired[, .(Mean = mean(Difference)), by = Dataset]
  confidence <- bootstrap_dataset_equal(
    paired, bootstrap_replicates, seed + index
  )
  data.table(
    Target_CellType = target_celltype,
    Comparator_CellType = comparator,
    Paired_Donors = nrow(paired),
    Datasets = uniqueN(paired$Dataset),
    Dataset_Equal_Mean_Log1p_CPM_Difference = mean(dataset_means$Mean),
    CI95_Lower = confidence[[1L]],
    CI95_Upper = confidence[[2L]]
  )
}), use.names = TRUE)

contract <- data.table(
  Field = c(
    "Gene", "Age_Interval", "Brain_Region", "Target_CellType",
    "Independent_Unit", "Minimum_Cells", "Minimum_Total_Counts",
    "Bootstrap_Replicates", "Dataset_Weighting", "Expression_Scale",
    "Cell_Use", "Interpretation_Boundary"
  ),
  Value = c(
    target_gene, target_age, target_region, target_celltype, "donor",
    minimum_cells, minimum_total_counts, bootstrap_replicates,
    "equal weight to each represented dataset",
    "log1p(CPM) from donor-by-cell-type raw-count pseudobulk",
    "all eligible cohort cells; no cell reduction",
    paste(
      "RNA expression context only; does not test regulatory-element",
      "activity, mediation or causality"
    )
  )
)

outputs <- c(
  "reuse_case_contract.tsv", "donor_celltype_pseudobulk.tsv.gz",
  "celltype_context_summary.tsv", "paired_donor_effects.tsv",
  "target_coverage.tsv"
)
fwrite(contract, file.path(output_dir, outputs[[1L]]), sep = "\t")
fwrite(
  donor_pseudobulk,
  file.path(output_dir, outputs[[2L]]),
  sep = "\t", compress = "gzip"
)
fwrite(context_summary, file.path(output_dir, outputs[[3L]]), sep = "\t")
fwrite(paired_effects, file.path(output_dir, outputs[[4L]]), sep = "\t")
fwrite(target_coverage, file.path(output_dir, outputs[[5L]]), sep = "\t")
manifest <- data.table(
  File = outputs,
  Size_Bytes = file.info(file.path(output_dir, outputs))$size
)
fwrite(manifest, file.path(output_dir, "manifest.tsv"), sep = "\t")
writeLines(
  c(
    "Status=complete",
    paste0("Cohort_Cells=", nrow(cohort)),
    paste0("Target_Cells=", sum(target_audit$Cells)),
    paste0("Target_Donors=", uniqueN(target_audit$Donor_ID)),
    paste0("Target_Datasets=", uniqueN(target_audit$Dataset))
  ),
  file.path(output_dir, "_SUCCESS")
)
unlink(file.path(output_dir, "_INCOMPLETE"))
message("DACT1 donor-level expression-context reuse case completed")

write_annotation_input(output_dir, setdiff(
  list.files(output_dir, pattern = "[.]tsv([.]gz)?$"), "annotation_input.tsv"
))
