#!/usr/bin/env Rscript
# Bounded replay of the recorded age and donor-distance calculations.
# The centroid input is retained output of evaluate_current_age_and_donor_structure.R.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: replay_age_from_centroids.R INPUT.rds NEW_OUTPUT_DIR")
if (dir.exists(args[2])) stop("Output directory must be new; existing results are never overwritten")
dir.create(args[2], recursive = TRUE)
out <- normalizePath(args[2])
input <- readRDS(args[1])
d <- as.data.table(input$donors)
cent <- input$centroids
ad <- as.data.table(input$all_donors)
allcent <- input$all_centroids
stopifnot(nrow(d) == 30L, length(cent) == 4L, length(allcent) == 4L)
pred <- folds <- list()
for (study in unique(d$Dataset)) {
  testall <- which(d$Dataset == study)
  train <- which(d$Dataset != study &
    !d$Canonical_Donor_ID %in% d$Canonical_Donor_ID[testall])
  test <- if (length(train)) testall[
    d$Age[testall] >= min(d$Age[train]) &
      d$Age[testall] <= max(d$Age[train])] else integer()
  ok <- length(train) >= 10L && length(test) >= 5L && uniqueN(d$Age[test]) >= 3L
  folds[[study]] <- data.table(Dataset = study, Train_Donors = length(train),
    Test_Donors = length(test), Outside_Training_Range = length(testall) - length(test),
    Included = ok)
  if (!ok) next
  for (method in names(cent)) {
    e <- cent[[method]]
    nr <- sqrt(rowSums(e^2))
    for (metric in c("euclidean", "cosine")) {
      dm <- if (metric == "euclidean") {
        outer(nr[test]^2, nr[train]^2, "+") -
          2 * tcrossprod(e[test, , drop = FALSE], e[train, , drop = FALSE])
      } else {
        1 - tcrossprod((e/nr)[test, , drop = FALSE], (e/nr)[train, , drop = FALSE])
      }
      pr <- apply(dm, 1, function(v) median(d$Age[train][order(v)[1:5]]))
      pred[[paste(study, method, metric)]] <- data.table(Dataset = study,
        Donor = d$Canonical_Donor_ID[test], Method = method, Distance = metric,
        Age = d$Age[test], Predicted_Age = pr, Median_Baseline = median(d$Age[train]))
    }
  }
}
fwrite(rbindlist(folds), file.path(out, "heldout_fold_coverage.tsv"), sep = "\t")
p <- rbindlist(pred)
p[, `:=`(Absolute_Error = abs(Age - Predicted_Age), Baseline_Error = abs(Age - Median_Baseline))]
fwrite(p, file.path(out, "age_predictions.tsv"), sep = "\t")
ds <- p[, .(Donors = .N, MAE = mean(Absolute_Error), Baseline_MAE = mean(Baseline_Error),
  Spearman = suppressWarnings(cor(Age, Predicted_Age, method = "spearman"))),
  by = .(Dataset, Method, Distance)]
fwrite(ds, file.path(out, "age_metrics_by_dataset.tsv"), sep = "\t")
fwrite(ds[, .(Datasets = .N, Donors = sum(Donors), Dataset_Equal_MAE = mean(MAE),
  Dataset_Equal_Baseline_MAE = mean(Baseline_MAE)), by = .(Method, Distance)],
  file.path(out, "age_metrics_summary.tsv"), sep = "\t")
struc <- list()
for (dataset in unique(ad$Dataset)) {
  ii <- which(ad$Dataset == dataset)
  if (length(ii) < 4L) next
  raw <- as.vector(dist(allcent$Raw[ii, , drop = FALSE]))
  for (method in names(allcent)) struc[[paste(dataset, method)]] <-
    data.table(Dataset = dataset, Donors = length(ii), Method = method,
      Spearman_Donor_Distance_to_Raw = suppressWarnings(cor(raw,
        as.vector(dist(allcent[[method]][ii, , drop = FALSE])), method = "spearman")))
}
fwrite(rbindlist(struc), file.path(out, "donor_structure_by_dataset.tsv"), sep = "\t")
message("Age predictions, study summaries and donor-distance correlations replayed from retained centroids")
