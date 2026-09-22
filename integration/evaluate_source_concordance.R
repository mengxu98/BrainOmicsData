#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
source("functions/utils.R")
setDTthreads(2L)
root <- "../../data/BrainOmicsData/integration_25/heldout_mapping"
validate_annotation_statistics(file.path(root, "rpca_mapping"))
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
        file = "../../data/BrainOmicsData/integration_25/annotation/celltype_assignments.rds",
        algo = "sha256", serialize = FALSE
      ), "held-out donor", "10000 donor bootstrap resamples, percentile 95% CI",
      "Majority vote > 0.5 and positive margin; descriptive support, not calibrated probability",
      "Source-defined absence from reference; excluded only from direct concordance, not from mapping or plots"
    )
  ),
  file.path(output_dir, "source_concordance_audit.tsv")
)
write_annotation_input(output_dir, setdiff(list.files(output_dir, pattern = "[.]tsv$"), "annotation_input.tsv"))
message("PASS: current RPCA donor-level mapping evaluation; no scVI-derived mapping metrics reused")
