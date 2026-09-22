suppressPackageStartupMessages({
  library(data.table)
})

source("functions/data_paths.R")
source("functions/dataset_metadata.R")
source("functions/processed_object.R")
source("functions/integration.R")


integration_dir <- brainomics_data_path("integration_25")
annotation_dir <- file.path(integration_dir, "annotation")
checkpoint_dir <- file.path(annotation_dir, "lisi_checkpoints")
label_file <- file.path(annotation_dir, "source_labels_harmonized.rds")
label_contract_file <- file.path(
  annotation_dir,
  "source_label_cell_contract.tsv"
)
metrics <- c("Source_CellClass", "Source_CellType")
checkpoint_files <- stats::setNames(
  file.path(
    checkpoint_dir,
    paste0(tolower(method_levels), "_source_label_lisi.rds")
  ),
  method_levels
)
required_files <- c(label_file, label_contract_file, checkpoint_files)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop(
    "Source-label cLISI aggregation inputs are missing: ",
    paste(missing_files, collapse = ", ")
  )
}

labels <- readRDS(label_file)
eligible <- labels$Source_Label_Evaluation_Eligible %in% TRUE &
  labels$Mapping_Status == "mapped_for_common_comparison" &
  !is.na(labels$Harmonized_CellClass) &
  !is.na(labels$Harmonized_CellType)
evaluation_labels <- as.data.table(labels[eligible, c(
  "Cells", "Dataset", "Donor_ID", "Harmonized_CellClass",
  "Harmonized_CellType"
)])
evaluation_labels[, Donor_ID := normalize_missing_metadata(Donor_ID)]
if (anyDuplicated(evaluation_labels$Cells)) {
  stop("Common source-label evaluation cells are duplicated")
}

label_contract <- read.delim(
  label_contract_file,
  sep = "\t",
  quote = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
contract_values <- stats::setNames(
  as.character(label_contract$Value),
  label_contract$Field
)
expected_cells <- as.integer(contract_values[["Mapped_Evaluation_Cells"]])
if (nrow(evaluation_labels) != expected_cells) {
  stop("Common source-label cohort differs from its cell contract")
}
expected_labels <- data.frame(Source_CellClass = evaluation_labels$Harmonized_CellClass,
  Source_CellType = evaluation_labels$Harmonized_CellType, row.names = evaluation_labels$Cells)

donor_rows <- vector("list", length(method_levels))
class_rows <- vector("list", length(method_levels))
type_rows <- vector("list", length(method_levels))
numeric_rows <- vector("list", length(method_levels))
for (method_index in seq_along(method_levels)) {
  method <- method_levels[[method_index]]
  thisutils::log_message("[source-label-lisi-aggregate] ", paste("Aggregating", method, "source-label cLISI"))
  score <- readRDS(checkpoint_files[[method]])
  if (!identical(names(score), metrics) ||
    nrow(score) != nrow(evaluation_labels) ||
    !identical(rownames(score), evaluation_labels$Cells) ||
    !identical(attr(score, "method", exact = TRUE), method) ||
    !identical(attr(score, "source_labels"), expected_labels) || any(!is.finite(as.matrix(score))) || any(as.matrix(score) < 1)) {
    stop(method, " source-label cLISI checkpoint failed aggregation checks")
  }
  score_dt <- cbind(
    evaluation_labels,
    as.data.table(score)
  )
  known <- !is.na(score_dt$Donor_ID)
  donor_wide <- score_dt[known, .(
    Cells = .N,
    Source_CellClass = mean(Source_CellClass),
    Source_CellType = mean(Source_CellType)
  ), by = .(Dataset, Donor_ID)]
  donor_long <- melt(
    donor_wide,
    id.vars = c("Dataset", "Donor_ID", "Cells"),
    measure.vars = metrics,
    variable.name = "Metric",
    value.name = "Donor_Mean_cLISI"
  )
  donor_long[, Method := method]
  donor_rows[[method_index]] <- donor_long

  class_rows[[method_index]] <- score_dt[known, .(
    Cells = .N,
    Donor_Mean_cLISI = mean(Source_CellClass)
  ), by = .(
    Dataset, Donor_ID,
    Harmonized_Label = Harmonized_CellClass
  )][, `:=`(Method = method, Metric = "Source_CellClass")]
  type_rows[[method_index]] <- score_dt[known, .(
    Cells = .N,
    Donor_Mean_cLISI = mean(Source_CellType)
  ), by = .(
    Dataset, Donor_ID,
    Harmonized_Label = Harmonized_CellType
  )][, `:=`(Method = method, Metric = "Source_CellType")]

  numeric_rows[[method_index]] <- data.frame(
    Method = method,
    Numeric_Precision = attr(score, "numeric_precision", exact = TRUE),
    R_Version = R.version.string,
    RANN_Version = as.character(packageVersion("RANN")),
    Perplexity = attr(score, "perplexity", exact = TRUE),
    Nearest_Neighbor_Error_Bound = attr(score, "nn_eps", exact = TRUE),
    Neighbor_Method = attr(score, "neighbor_method", exact = TRUE),
    Neighbor_Parameters = attr(score, "neighbor_parameters", exact = TRUE),
    Neighbor_Package_Version = attr(score, "neighbor_package_version", exact = TRUE),
    attr(score, "neighbor_accuracy", exact = TRUE),
    LISI_Source_Commit = attr(score, "lisi_source_commit", exact = TRUE),
    stringsAsFactors = FALSE
  )
  rm(score, score_dt, donor_wide, donor_long)
  gc()
}

donor_scores <- rbindlist(donor_rows, use.names = TRUE)
setcolorder(
  donor_scores,
  c("Dataset", "Donor_ID", "Method", "Metric", "Cells", "Donor_Mean_cLISI")
)
dataset_scores <- donor_scores[, .(
  Known_Donors = .N,
  Cells = sum(Cells),
  Donor_Equal_Mean_cLISI = mean(Donor_Mean_cLISI),
  Donor_Equal_Median_cLISI = stats::median(Donor_Mean_cLISI)
), by = .(Dataset, Method, Metric)]
dataset_scores[, Method_Order := match(Method, method_levels)]
setorder(dataset_scores, Metric, Dataset, Method_Order)
dataset_scores[, Method_Order := NULL]

effect_rows <- list()
effect_index <- 0L
for (metric in metrics) {
  raw <- dataset_scores[
    Method == "Raw" & Metric == metric,
    .(Dataset, Raw_cLISI = Donor_Equal_Mean_cLISI)
  ]
  for (method in method_levels) {
    effect_index <- effect_index + 1L
    comparison <- merge(
      raw,
      dataset_scores[
        Method == method & Metric == metric,
        .(Dataset, Method_cLISI = Donor_Equal_Mean_cLISI)
      ],
      by = "Dataset"
    )
    difference <- comparison$Method_cLISI - comparison$Raw_cLISI
    n <- length(difference)
    estimate <- mean(difference)
    standard_error <- if (n > 1L) stats::sd(difference) / sqrt(n) else NA_real_
    critical <- if (n > 1L) stats::qt(0.975, df = n - 1L) else NA_real_
    effect_rows[[effect_index]] <- data.frame(
      Metric = metric,
      Method = method,
      Reference_Method = "Raw",
      Independent_Unit = "source dataset",
      Datasets = n,
      Mean_Paired_Difference = estimate,
      CI95_Lower = estimate - critical * standard_error,
      CI95_Upper = estimate + critical * standard_error,
      Median_Paired_Difference = stats::median(difference),
      IQR_Paired_Difference = stats::IQR(difference),
      Direction = "negative difference indicates better label preservation",
      stringsAsFactors = FALSE
    )
  }
}
effects <- rbindlist(effect_rows, use.names = TRUE)

label_donor_scores <- rbindlist(
  c(class_rows, type_rows),
  use.names = TRUE,
  fill = TRUE
)
label_summary <- label_donor_scores[, .(
  Cells = sum(Cells),
  Known_Donors = .N,
  Datasets = uniqueN(Dataset),
  Donor_Equal_Mean_cLISI = mean(Donor_Mean_cLISI),
  Donor_Equal_Median_cLISI = stats::median(Donor_Mean_cLISI)
), by = .(Method, Metric, Harmonized_Label)]
label_summary[, Method_Order := match(Method, method_levels)]
setorder(label_summary, Metric, Harmonized_Label, Method_Order)
label_summary[, Method_Order := NULL]

output_files <- c(
  file.path(annotation_dir, "source_lisi_by_donor.tsv.gz"),
  file.path(annotation_dir, "source_lisi_by_dataset.tsv"),
  file.path(annotation_dir, "source_lisi_effects.tsv"),
  file.path(annotation_dir, "source_lisi_by_label_and_donor.tsv.gz"),
  file.path(annotation_dir, "source_lisi_by_label.tsv"),
  file.path(annotation_dir, "source_lisi_numeric_audit.tsv"),
  file.path(annotation_dir, "source_lisi_contract.tsv")
)
manifest_file <- file.path(annotation_dir, "source_lisi_manifest.tsv")
incomplete_file <- file.path(annotation_dir, "source_lisi_INCOMPLETE")
success_file <- file.path(annotation_dir, "source_lisi_SUCCESS")
if (any(file.exists(c(output_files, manifest_file, success_file)))) {
  stop("Source-label cLISI aggregation refuses to overwrite formal outputs")
}
writeLines(
  c(
    paste0("PID=", Sys.getpid()),
    paste0("Started_UTC=", format(Sys.time(), tz = "UTC", usetz = TRUE))
  ),
  incomplete_file,
  useBytes = TRUE
)

write_tsv(
  as.data.frame(donor_scores),
  output_files[[1L]]
)
write_tsv(
  as.data.frame(dataset_scores),
  output_files[[2L]]
)
write_tsv(
  as.data.frame(effects),
  output_files[[3L]]
)
write_tsv(
  as.data.frame(label_donor_scores),
  output_files[[4L]]
)
write_tsv(
  as.data.frame(label_summary),
  output_files[[5L]]
)
write_tsv(
  do.call(rbind, numeric_rows),
  output_files[[6L]]
)
write_tsv(
  data.frame(
    Metric = "original-study-label cLISI",
    Space = "method-specific 50-dimensional latent representation",
    Cells = nrow(evaluation_labels),
    Source_Label_Datasets = uniqueN(evaluation_labels$Dataset),
    Known_Donors = uniqueN(evaluation_labels$Donor_ID[!is.na(
      evaluation_labels$Donor_ID
    )]),
    Source_Label_Levels = paste(metrics, collapse = ";"),
    Interpretation = "lower cLISI indicates less local mixing of the fixed source labels",
    Inference_Unit = "source dataset",
    Cell_Level_Use = "descriptive only; no cell-level hypothesis test",
    Cohort = paste(
      "all cells with reviewed original-study common labels;",
      "source unknown, mixed, low-quality and unavailable labels excluded;",
      "no random selection"
    ),
    stringsAsFactors = FALSE
  ),
  output_files[[7L]]
)

if (any(!file.exists(output_files)) || any(file.info(output_files)$size <= 0)) {
  stop("Source-label cLISI aggregation did not produce every formal output")
}
manifest <- data.frame(
  File = basename(output_files),
  Size_Bytes = as.numeric(file.info(output_files)$size),
  stringsAsFactors = FALSE
)
write_tsv(manifest, manifest_file)
writeLines(
  c(
    "Status=complete",
    paste0("Completed_UTC=", format(Sys.time(), tz = "UTC", usetz = TRUE))
  ),
  success_file,
  useBytes = TRUE
)
unlink(incomplete_file)

thisutils::log_message("[source-label-lisi-aggregate] ", paste(
  "Completed source-label cLISI aggregation for",
  format(nrow(evaluation_labels), big.mark = ","),
  "cells and",
  uniqueN(evaluation_labels$Dataset),
  "source-labelled datasets"
))
