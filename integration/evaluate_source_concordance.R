#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
source("functions/utils.R")
source("functions/data_paths.R")
setDTthreads(2L)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% c(0L, 3L))
if (length(args) == 3L) {
  input <- normalizePath(args[1], mustWork = TRUE)
  map_file <- normalizePath(args[2], mustWork = TRUE)
  out <- args[3]
  if (dir.exists(out)) stop("Output directory already exists: ", out)
  d <- fread(input)
  class_map <- fread(map_file, na.strings = "")
  stopifnot(nrow(d) == 1661798L, !anyDuplicated(d$CellID), uniqueN(d$Donor) == 26L,
    all(d$Mapping_Status == "mapped"), !anyDuplicated(class_map$Source_CellClass),
    all(d$Source_CellClass %in% class_map$Source_CellClass))
  d[, Source_Comparable := class_map$Source_Comparable[match(Source_CellClass, class_map$Source_CellClass)]]
  d[Source_Comparable == "", Source_Comparable := NA_character_]
  restricted <- restrict_reference_labels(d$Working_CellType, !is.na(d$Source_Comparable))
  d[, `:=`(Restricted_Label = restricted$Prediction, Restricted_Status = restricted$Mapping_Status)]
  d[, Candidate_Broad := restrict_reference_labels(Working_CellType, rep(TRUE, .N))$Prediction]
  d[, `:=`(Comparable = !is.na(Source_Comparable), Assigned = !is.na(Restricted_Label),
    Known_OutOfReference = Source_CellClass == "Erythrocyte")]
  d[, Candidate_Error := Comparable & (is.na(Candidate_Broad) | Candidate_Broad != Source_Comparable)]
  d[, Restricted_Error := Assigned & Restricted_Label != Source_Comparable]
  stopifnot(sum(d$Comparable) == 605701L, sum(!d$Comparable) == 1056097L,
    all(is.na(d$Restricted_Label[!d$Comparable])),
    !anyNA(d$Candidate_Error), !anyNA(d$Restricted_Error))
  summarize <- function(fields) d[, .(
    Cells = .N, Comparable_Cells = sum(Comparable), Assigned_Cells = sum(Assigned),
    Unsupported_Cells = sum(!Assigned), Unsupported_Rate = mean(!Assigned),
    Candidate_Error_Cells = if (any(Comparable)) sum(Candidate_Error[Comparable]) else NA_integer_,
    Candidate_Error_Rate = if (any(Comparable)) mean(Candidate_Error[Comparable]) else NA_real_,
    Assigned_Error_Cells = if (any(Assigned)) sum(Restricted_Error[Assigned]) else NA_integer_,
    Assigned_Error_Rate = if (any(Assigned)) mean(Restricted_Error[Assigned]) else NA_real_,
    Correct_Assigned_Per_All_Comparable = if (any(Comparable)) sum(Assigned & !Restricted_Error) / sum(Comparable) else NA_real_,
    Known_OutOfReference_Cells = sum(Known_OutOfReference),
    Forced_OutOfReference_Candidate_Errors = if (any(Known_OutOfReference)) sum(Known_OutOfReference) else NA_integer_,
    Known_OutOfReference_Candidate_Error_Rate = if (any(Known_OutOfReference)) 1 else NA_real_,
    Candidate_Inhibitory_Cells = sum(Working_CellType == "Inhibitory neurons"),
    Candidate_Inhibitory_Support_GE_0_7 = sum(Working_CellType == "Inhibitory neurons" & Neighbor_Vote_Support >= 0.7)
  ), by = fields]
  reports <- list(
    restricted_by_source = summarize("Source_CellClass"),
    restricted_by_donor_source = summarize(c("Donor", "Source_CellClass")),
    restricted_by_donor = summarize("Donor")
  )
  for (report in reports) stopifnot(
    sum(report$Cells) == nrow(d),
    sum(report$Comparable_Cells) == sum(d$Comparable),
    sum(report$Assigned_Cells) == sum(d$Assigned),
    sum(report$Unsupported_Cells) == sum(!d$Assigned),
    sum(report$Known_OutOfReference_Cells) == sum(d$Known_OutOfReference),
    sum(report$Forced_OutOfReference_Candidate_Errors, na.rm = TRUE) == sum(d$Known_OutOfReference),
    sum(report$Candidate_Inhibitory_Cells) == sum(d$Working_CellType == "Inhibitory neurons"),
    sum(report$Candidate_Inhibitory_Support_GE_0_7) == sum(d$Working_CellType == "Inhibitory neurons" & d$Neighbor_Vote_Support >= 0.7)
  )
  dir.create(out, recursive = TRUE)
  for (name in names(reports)) fwrite(reports[[name]], file.path(out, paste0(name, ".tsv")), sep = "\t", na = "NA")
  fwrite(d[, .(Cells = .N), by = .(Source_CellClass, Working_CellType, Restricted_Label, Restricted_Status)],
    file.path(out, "restricted_source_prediction.tsv"), sep = "\t", na = "NA")
  fwrite(d[, .(CellID, Donor, Source_CellClass, Candidate_Label = Working_CellType,
    Neighbor_Vote_Support, Restricted_Label, Restricted_Status, UMAP_1, UMAP_2)],
    file.path(out, "restricted_query_predictions.tsv.gz"), sep = "\t", na = "NA")
  quit(save = "no")
}

root <- file.path(brainomics_results_dir(), "heldout_mapping")
mapping <- fread(file.path(root, "rpca_mapping/query_rpca_mapping.tsv.gz"))
stopifnot(
  nrow(mapping) == 1655074L, !anyDuplicated(mapping$Cells),
  uniqueN(mapping$Global_Donor_ID) == 26L, !anyNA(mapping$Majority_Supported),
  all(mapping$Vote_Support >= 0 & mapping$Vote_Support <= 1),
  all(mapping$Vote_Margin >= 0 & mapping$Vote_Margin <= 1)
)
output_dir <- file.path(root, "source_concordance")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
donors <- mapping[, .(
  Cells = .N,
  Comparable_Cells = sum(Comparable),
  Comparable_Concordance = if (any(Comparable)) mean(Common_Class_Concordant[Comparable]) else NA_real_,
  Majority_Support_Rate = mean(Majority_Supported),
  Source_OutOfReference_Fraction = mean(!Comparable)
), by = Global_Donor_ID]
set.seed(20260730L)
summary <- rbindlist(lapply(c(
  "Comparable_Concordance", "Majority_Support_Rate",
  "Source_OutOfReference_Fraction"
), function(metric) {
  values <- donors[[metric]]
  values <- values[is.finite(values)]
  boot <- replicate(10000L, mean(sample(values, length(values), replace = TRUE)))
  interval <- quantile(boot, c(0.025, 0.975), names = FALSE)
  data.table(
    Metric = metric, Donors = length(values), Estimate = mean(values),
    CI_Lower = interval[1], CI_Upper = interval[2]
  )
}))
status <- mapping[, .(Cells = .N), by = .(
  Mapping_Support =
    fifelse(Majority_Supported, "majority_supported", "low_support")
)]
confusion <- mapping[, .(Cells = .N), by = .(
  Source_CellType_Display,
  Source_Common_Class, Comparison_Status, Predicted_CellType
)]
source_summary <- mapping[, .(
  Cells = .N, Donors = uniqueN(Global_Donor_ID),
  Mean_Vote_Support = mean(Vote_Support), Majority_Support_Rate = mean(Majority_Supported),
  Comparable = all(Comparable), Common_Class_Concordance =
    if (any(Comparable)) mean(Common_Class_Concordant[Comparable]) else NA_real_
),
by = Source_CellType_Display
]
prediction_summary <- mapping[, .(
  Cells = .N, Donors = uniqueN(Global_Donor_ID),
  Mean_Vote_Support = mean(Vote_Support), Median_Vote_Support = median(Vote_Support),
  Mean_Vote_Margin = mean(Vote_Margin), Majority_Support_Rate = mean(Majority_Supported)
),
by = Predicted_CellType
]
fwrite(donors, file.path(output_dir, "source_concordance_by_donor.tsv"), sep = "\t")
fwrite(summary, file.path(output_dir, "donor_equal_metric_summary.tsv"), sep = "\t")
fwrite(status, file.path(output_dir, "mapping_support_summary.tsv"), sep = "\t")
fwrite(confusion, file.path(output_dir, "source_reference_confusion.tsv"), sep = "\t")
fwrite(source_summary, file.path(output_dir, "source_mapping_summary.tsv"), sep = "\t")
fwrite(prediction_summary, file.path(output_dir, "predicted_celltype_summary.tsv"), sep = "\t")

# Expression evidence is summarized after mapping, never used to overwrite labels.
query_dir <- file.path(root, "egad_query_input")
manifest <- fread(file.path(query_dir, "manifest.tsv"))
normalization <- fread(file.path(query_dir, "normalization_manifest.tsv"))
stopifnot(identical(manifest$Shard_ID, normalization$Shard_ID))
marker_rows <- list()
for (i in seq_len(nrow(manifest))) {
  directory <- dirname(file.path(query_dir, manifest$Cells_File[i]))
  marker_file <- file.path(directory, "marker_counts.rds")
  total_file <- file.path(directory, "full_rna_library_sizes.rds")
  stopifnot(
    digest::digest(file = marker_file, algo = "sha256", serialize = FALSE) == normalization$Marker_Counts_SHA256[i],
    digest::digest(file = total_file, algo = "sha256", serialize = FALSE) == normalization$Library_Sizes_SHA256[i]
  )
  counts <- readRDS(marker_file)
  totals <- readRDS(total_file)
  index <- match(colnames(counts), mapping$Cells)
  stopifnot(
    !anyNA(index), identical(names(totals), colnames(counts)),
    length(index) == manifest$Analysis_Cells[i], uniqueN(mapping$Global_Donor_ID[index]) == 1L
  )
  donor <- mapping$Global_Donor_ID[index[1L]]
  counts@x <- log1p(counts@x * rep.int(10000 / totals, diff(counts@p)))
  label_rows <- list()
  for (label_system in c("Predicted_CellType", "Source_CellType_Display")) {
    labels <- mapping[[label_system]][index]
    for (label in unique(labels)) {
      selected <- labels == label
      values <- counts[, selected, drop = FALSE]
      label_rows[[length(label_rows) + 1L]] <- data.table(
        Global_Donor_ID = donor,
        Label_System = label_system, CellType = label, Gene = rownames(values),
        Measured_Cells = sum(selected), SumLog = Matrix::rowSums(values),
        Positive_Cells = Matrix::rowSums(values > 0)
      )
    }
  }
  marker_rows[[i]] <- rbindlist(label_rows)
  rm(counts, totals, values, label_rows)
  gc()
}
marker_donor <- rbindlist(marker_rows)
marker_donor[, `:=`(
  MeanLog = SumLog / Measured_Cells,
  PctPositive = Positive_Cells / Measured_Cells
)]
marker_summary <- marker_donor[, .(
  Measured_Donors = uniqueN(Global_Donor_ID),
  Measured_Cells = sum(Measured_Cells), MeanLog = sum(SumLog) / sum(Measured_Cells),
  PctPositive = sum(Positive_Cells) / sum(Measured_Cells),
  Donor_Equal_MeanLog = mean(MeanLog), Donor_Equal_PctPositive = mean(PctPositive)
),
by = .(Label_System, CellType, Gene)
]
stopifnot(
  all(is.finite(marker_summary$MeanLog)),
  all(marker_summary$PctPositive >= 0 & marker_summary$PctPositive <= 1)
)
fwrite(marker_donor, file.path(output_dir, "marker_expression_by_donor.tsv"), sep = "\t")
fwrite(marker_summary, file.path(output_dir, "marker_expression_by_celltype.tsv"), sep = "\t")
fwrite(normalization, file.path(output_dir, "marker_measured_gene_coverage.tsv"), sep = "\t")
write_tsv(
  data.frame(
    Field = c(
      "Mapping_Method", "Query_Cells", "Reference_Annotation_SHA256",
      "Inference_Unit", "Confidence_Interval", "Support_Rule", "OutOfReference_Rule"
    ),
    Value = c(
      "RPCA KNN, k=30, cosine, Annoy", nrow(mapping),
      digest::digest(
        file = "results/annotation/celltype_assignments.rds",
        algo = "sha256", serialize = FALSE
      ), "held-out donor", "10000 donor bootstrap resamples, percentile 95% CI",
      "Majority vote > 0.5 and positive margin; descriptive support, not calibrated probability",
      "Source-defined absence from reference; excluded only from direct concordance, not from mapping or plots"
    )
  ),
  file.path(output_dir, "source_concordance_parameters.tsv")
)
write_result_manifest(output_dir, setdiff(list.files(output_dir, pattern = "[.]tsv$"), "result_manifest.tsv"))
message("Source-label concordance results written")
