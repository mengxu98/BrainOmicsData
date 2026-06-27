#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1]]
}

lisi_results_file <- value_after("--lisi-results-file", "results/lisi/lisi_results.rds")
out_dir <- value_after(
  "--out-dir",
  "../../data/ScienceDB/HumanBrain_sc_snRNAseq_age_resource"
)

qc_dir <- file.path(out_dir, "06_quality_control")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, file) {
  utils::write.table(
    x, file = file, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE, na = "missing"
  )
}

if (file.exists(lisi_results_file)) {
  lisi <- readRDS(lisi_results_file)
  lisi <- as.data.frame(lisi)
  summary <- do.call(rbind, lapply(names(lisi), function(method) {
    values <- as.numeric(lisi[[method]])
    data.frame(
      method = method,
      median = stats::median(values, na.rm = TRUE),
      q25 = as.numeric(stats::quantile(values, 0.25, na.rm = TRUE)),
      q75 = as.numeric(stats::quantile(values, 0.75, na.rm = TRUE)),
      mean = mean(values, na.rm = TRUE),
      max = max(values, na.rm = TRUE),
      n_cells_or_points = sum(!is.na(values)),
      status = "computed",
      source_file = lisi_results_file,
      stringsAsFactors = FALSE
    )
  }))
} else {
  summary <- data.frame(
    method = c("Raw", "Harmony", "RPCA"),
    median = "missing",
    q25 = "missing",
    q75 = "missing",
    mean = "missing",
    max = "missing",
    n_cells_or_points = "missing",
    status = "missing_lisi_results_rds; run bash 04_datasets_plotting.sh first",
    source_file = lisi_results_file,
    stringsAsFactors = FALSE
  )
}

write_tsv(summary, file.path(qc_dir, "lisi_before_after_integration.tsv"))
write_tsv(summary, file.path(qc_dir, "lisi_before_after_rpca.tsv"))

batch_summary <- data.frame(
  metric = "dataset_label_lisi",
  comparison = "Raw_vs_Harmony_vs_RPCA",
  summary_file = "06_quality_control/lisi_before_after_integration.tsv",
  interpretation = "Higher LISI indicates stronger mixing of source-dataset labels in the embedding.",
  status = if (all(summary$status == "computed")) "computed" else "missing_lisi_results",
  stringsAsFactors = FALSE
)
write_tsv(batch_summary, file.path(qc_dir, "batch_integration_summary.tsv"))

message("LISI summary written to: ", qc_dir)
