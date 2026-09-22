#!/usr/bin/env Rscript
# Donor-level age validation of the existing, age-unsupervised atlas embeddings.
suppressPackageStartupMessages(library(data.table))
source("functions/data_paths.R")
source("functions/metadata_schema.R")
source("functions/dataset_metadata.R")
setDTthreads(8L)
root <- brainomics_data_path("integration_25")
out <- file.path(root, "evaluation", "age_signal_preservation")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
args <- commandArgs(trailingOnly = TRUE)

if ("--expression-prepare" %in% args) {
  # Reuse the fixed age cohort; do not choose populations from expression results.
  suppressPackageStartupMessages(library(SeuratObject))
  source("functions/processed_object.R")
  source("functions/integration.R")
  d <- fread(file.path(out, "donor_metadata.tsv"))[
    CellType == "Excitatory neurons" & BrainRegion == "Frontal cortex" & Unit == "Years"]
  stopifnot(nrow(d) > 0L, !anyDuplicated(d$Global_Donor_ID))
  a <- as.data.table(readRDS(file.path(root, "annotation", "celltype_assignments.rds")))
  selected_cells <- a[CellType == "Excitatory neurons", Cells]
  rm(a)
  counts_by_study <- list()
  covariates <- list()
  for (study in unique(d$Dataset)) {
    message("Preparing donor pseudobulk: ", study)
    folder <- brainomics_data_path("processed", study)
    format <- fread(file.path(folder, "processed_object_format.tsv"))
    m <- fread(file.path(folder, format$Canonical_Metadata[[1L]]), select = c(
      "Cells", "Original_Cell_ID", "Global_Donor_ID", "BrainRegion", "Age_num", "Unit",
      "Sex", "Assay_Type", "Sequencing_Platform", "Library_Chemistry", "source_cell_type_label"))
    m[, Qualified_Cell := paste(study, Original_Cell_ID, sep = "::")]
    m[, BrainRegion := standardize_brain_region(BrainRegion)]
    ds <- d[Dataset == study]
    m <- m[Qualified_Cell %in% selected_cells & Global_Donor_ID %in% ds$Global_Donor_ID &
      BrainRegion == "Frontal cortex" & Unit == "Years"]
    stopifnot(!anyDuplicated(m$Cells), !anyDuplicated(m$Qualified_Cell),
      all(abs(m$Age_num - ds$Age[match(m$Global_Donor_ID, ds$Global_Donor_ID)]) < 1e-8))
    observed <- m[, .N, by = Global_Donor_ID]
    stopifnot(identical(observed$N[match(ds$Global_Donor_ID, observed$Global_Donor_ID)], ds$Cells))
    covariates[[study]] <- m[, .(Cells = .N, Sex = paste(sort(unique(Sex)), collapse = ";"),
      Assay = paste(sort(unique(Assay_Type)), collapse = ";"),
      Platform = paste(sort(unique(Sequencing_Platform)), collapse = ";"),
      Chemistry = paste(sort(unique(Library_Chemistry)), collapse = ";"),
      Source_Labels = uniqueN(source_cell_type_label)), by = .(Global_Donor_ID)]
    object <- load_processed_object(file.path(folder, format$Full_Object[[1L]]))
    sums <- matrix(0, nrow(object[["RNA"]]), nrow(ds),
      dimnames = list(rownames(object[["RNA"]]), ds$Global_Donor_ID))
    seen <- character()
    for (layer in Layers(object, assay = "RNA", search = "^counts")) {
      cells <- intersect(Cells(object[["RNA"]], layer = layer), m$Cells)
      if (!length(cells)) next
      x <- LayerData(object, assay = "RNA", layer = layer, cells = cells)
      donor <- m$Global_Donor_ID[match(colnames(x), m$Cells)]
      design <- Matrix::sparseMatrix(i = seq_along(donor), j = match(donor, ds$Global_Donor_ID),
        x = 1, dims = c(length(donor), nrow(ds)))
      sums[rownames(x), ] <- sums[rownames(x), , drop = FALSE] + as.matrix(x %*% design)
      seen <- c(seen, colnames(x))
    }
    stopifnot(!anyDuplicated(seen), setequal(seen, m$Cells), all(colSums(sums) > 0))
    crosswalk <- feature_crosswalk(read_feature_metadata(file.path(folder, format$Feature_Metadata[[1L]])),
      rownames(sums), study)
    rownames(sums) <- crosswalk$Canonical_Feature_ID
    counts_by_study[[study]] <- sums
    rm(object, x, sums, m)
    gc()
  }
  genes <- Reduce(intersect, lapply(counts_by_study, rownames))
  stopifnot(length(genes) > 0L)
  # Library sizes use all measured genes, before intersecting source feature sets.
  sizes <- unlist(lapply(counts_by_study, colSums), use.names = FALSE)
  counts <- do.call(cbind, lapply(counts_by_study, function(x) x[genes, , drop = FALSE]))
  cpm <- sweep(counts, 2, sizes, "/") * 1e6
  d <- d[match(colnames(counts), Global_Donor_ID)]
  stopifnot(identical(d$Global_Donor_ID, colnames(counts)), all(is.finite(cpm)))
  saveRDS(list(donors = d, counts = counts, cpm = cpm, library_sizes = sizes),
    file.path(out, "donor_expression.rds"))
  fwrite(merge(d, rbindlist(covariates), by = "Global_Donor_ID", suffixes = c("", "_observed")),
    file.path(out, "expression_cohort.tsv"), sep = "\t")
  message("Expression preparation complete: ", nrow(d), " donors; ", length(genes), " common genes")
}

if ("--expression-check" %in% args) {
  e <- readRDS(file.path(out, "donor_expression.rds"))
  d <- e$donors
  cov <- fread(file.path(out, "expression_cohort.tsv"))
  d[, Sex := cov$Sex[match(Global_Donor_ID, cov$Global_Donor_ID)]]
  # Fixed from coverage before expression inspection: largest discovery study,
  # pediatric replication study, and their shared pediatric age support.
  studies <- c("Velmeshev_2023", "GSE204683")
  indices <- lapply(studies, function(s) which(d$Dataset == s & between(d$Age, 1, 14)))
  stopifnot(length(indices[[1L]]) >= 6L, length(indices[[2L]]) >= 5L)
  expressed <- rowSums(e$cpm[, indices[[1L]], drop = FALSE] >= 1) >= ceiling(length(indices[[1L]]) / 2)
  y <- log2(e$cpm[expressed, , drop = FALSE] + 1)
  results <- list()
  designs <- character()
  for (i in seq_along(studies)) {
    ii <- indices[[i]]
    dd <- as.data.frame(d[ii])
    dd$Age_log <- log1p(dd$Age)
    dd$Sex <- normalize_missing_metadata(dd$Sex)
    known_sex <- !anyNA(dd$Sex) && all(tolower(dd$Sex) %in% c("female", "male", "f", "m"))
    formula <- if (known_sex && length(unique(dd$Sex)) == 2L) ~ Age_log + Sex else ~ Age_log
    design <- model.matrix(formula, dd)
    stopifnot(qr(design)$rank == ncol(design), nrow(design) > ncol(design) + 1L)
    fit <- lm.fit(design, t(y[, ii, drop = FALSE]))
    df <- nrow(design) - fit$rank
    variance <- colSums(fit$residuals^2) / df
    se <- sqrt(variance * solve(crossprod(design))["Age_log", "Age_log"])
    beta <- fit$coefficients["Age_log", ]
    p <- 2 * pt(-abs(beta / se), df)
    p[!is.finite(p)] <- 1
    results[[i]] <- data.table(Gene = rownames(y), Study = studies[[i]], Donors = length(ii),
      Beta = beta, SE = se, CI_Lower = beta - qt(0.975, df) * se,
      CI_Upper = beta + qt(0.975, df) * se, P = p, FDR = p.adjust(p, "BH"))
    designs <- c(designs, paste(studies[[i]], paste(deparse(formula), collapse = " "),
      "n =", length(ii), "; assay/platform/chemistry recorded in expression_cohort.tsv"))
  }
  discovery <- results[[1L]][FDR < 0.05]
  replication <- results[[2L]][match(discovery$Gene, Gene)]
  replication[, Replication_FDR := p.adjust(P, "BH")]
  replicated <- sum(replication$Replication_FDR < 0.05 & sign(replication$Beta) == sign(discovery$Beta))
  fwrite(rbindlist(results), file.path(out, "expression_age_associations.tsv.gz"), sep = "\t")
  report <- c("Exploratory donor pseudobulk check; not proof of embedding preservation.",
    "Fixed cohort: frontal excitatory neurons, ages 1-14 years; log2(CPM+1) ~ log1p(age) (+ sex when known).",
    "Genes: CPM >= 1 in at least half of discovery donors; discovery BH FDR < 0.05.",
    "Replication: two-sided age test, BH within discovery candidates, concordant direction required.",
    designs, paste("Tested genes:", nrow(y)), paste("Discovery genes:", nrow(discovery)),
    if (nrow(discovery)) paste("Replicated candidates:", replicated)
    else paste("No candidates proceeded to formal cross-study replication assessment because none passed the discovery threshold.",
      "Coefficients for both studies are retained as exploratory source data."),
    "Low donor count and broad-cell-type composition limit inference; no proof of absent age biology if this check fails.",
    if (replicated == 0L) "STOP: this check does not establish a reproducible signal for embedding-retention claims."
    else "Candidate associations require composition/confounding review before an embedding-retention analysis.")
  writeLines(report, file.path(out, "expression_check.txt"))
  message(paste(report, collapse = "\n"))
}

if ("--prepare" %in% args) {
  message("Reading full metadata (not the age-empty plotting cache)")
  m <- as.data.table(readRDS(file.path(root, "metadata_filtered.rds")))
  fields <- c("Cells", "Dataset", "Global_Donor_ID", "Donor_ID_Verification_Status",
    "BrainRegion", "Age", "Age_num", "Age_Lower", "Age_Upper", "Unit",
    "Age_Representative_Method", "Continuous_Age_Eligibility")
  m <- m[, ..fields]
  gc()
  m[, BrainRegion := standardize_brain_region(BrainRegion)]
  a <- as.data.table(readRDS(file.path(root, "annotation", "celltype_assignments.rds")))
  stopifnot(!anyDuplicated(m$Cells), !anyDuplicated(a$Cells),
    nrow(m) == nrow(a), setequal(m$Cells, a$Cells))
  m <- m[match(a$Cells, Cells)]
  m[, CellType := a$CellType]
  rm(a)
  m[, eligible := !is.na(normalize_missing_metadata(Global_Donor_ID)) &
    !grepl("unknown", Donor_ID_Verification_Status, ignore.case = TRUE) &
    is.finite(Age_num) & is.finite(Age_Lower) & is.finite(Age_Upper) &
    Age_Lower == Age_Upper & Age_num == Age_Lower &
    !grepl("not eligible|mean|midpoint|range", paste(Continuous_Age_Eligibility,
      Age_Representative_Method, Age), ignore.case = TRUE) & Dataset != "GSE144136"]
  m[is.na(eligible), eligible := FALSE]
  fwrite(m[, .(Cells = .N), by = .(Dataset, Unit, eligible)], file.path(out, "age_eligibility.tsv"), sep = "\t")
  donor_check <- m[eligible == TRUE, .(Ages = uniqueN(paste(Age_num, Unit))), by = Global_Donor_ID]
  m[Global_Donor_ID %in% donor_check[Ages != 1L, Global_Donor_ID], eligible := FALSE]
  fwrite(donor_check[Ages != 1L], file.path(out, "donor_age_conflicts.tsv"), sep = "\t")
  m[, group := .GRP, by = .(Dataset, Global_Donor_ID, BrainRegion, CellType)]
  donors <- m[eligible == TRUE, .(Cells = .N, Age = Age_num[[1]], Unit = Unit[[1]]),
    by = .(group, Dataset, Global_Donor_ID, BrainRegion, CellType)]
  donors <- donors[Cells >= 20L]
  setorder(donors, group)
  fwrite(donors, file.path(out, "donor_metadata.tsv"), sep = "\t")
  coverage <- donors[, .(Donors = .N, Cells = sum(Cells), Ages = uniqueN(Age),
    Age_Min = min(Age), Age_Max = max(Age)), by = .(CellType, BrainRegion, Unit, Dataset)]
  fwrite(coverage, file.path(out, "coverage.tsv"), sep = "\t")
  selected <- which(m$eligible & m$group %in% donors$group)
  group <- match(m$group[selected], donors$group)
  cell_ids <- m$Cells
  rm(m)
  gc()
  centroids <- list()
  for (method in c("Raw", "RPCA", "Harmony", "scVI")) {
    message("Aggregating ", method, " donor centroids")
    emb <- readRDS(file.path(root, "annotation", "reductions", paste0(tolower(method), "_latent.rds")))
    stopifnot(identical(rownames(emb), cell_ids), ncol(emb) == 50L)
    sums <- rowsum(emb[selected, , drop = FALSE], group, reorder = TRUE)
    stopifnot(identical(as.integer(rownames(sums)), seq_len(nrow(donors))), all(is.finite(sums)))
    centroids[[method]] <- sums / donors$Cells
    rm(emb, sums)
    gc()
  }
  saveRDS(list(donors = donors, centroids = centroids), file.path(out, "donor_centroids.rds"))
  message("Preparation complete: ", out)
}

if ("--evaluate" %in% args) {
  # Select the population from coverage alone, before examining predictions.
  at <- match("--evaluate", args)
  celltype <- args[at + 1L]
  region <- args[at + 2L]
  unit <- args[at + 3L]
  inputs <- readRDS(file.path(out, "donor_centroids.rds"))
  d <- inputs$donors
  # The two GSE178175 source matrices have unresolved sample identity; do not
  # treat them as independent age-labelled donors. Resource cells are retained.
  indices <- which(d$CellType == celltype & d$BrainRegion == region & d$Unit == unit &
    d$Dataset != "GSE178175")
  d <- copy(d[indices])
  stopifnot(uniqueN(d$Dataset) >= 2L, !anyDuplicated(d$Global_Donor_ID))
  k <- 5L
  repeats <- 1000L
  set.seed(2026)
  predictions <- list()
  folds <- list()
  for (study in unique(d$Dataset)) {
    train <- which(d$Dataset != study)
    test_all <- which(d$Dataset == study)
    # Evaluate interpolation only; report discarded extrapolation separately.
    test <- test_all[d$Age[test_all] >= min(d$Age[train]) &
      d$Age[test_all] <= max(d$Age[train])]
    folds[[study]] <- data.table(Dataset = study, Training_Donors = length(train),
      Test_Donors = length(test), Outside_Training_Range = length(test_all) - length(test),
      Training_Age_Min = min(d$Age[train]), Training_Age_Max = max(d$Age[train]),
      Included_Test = length(train) >= 10L && length(test) >= 5L && uniqueN(d$Age[test]) >= 3L)
    if (length(train) < 10L || length(test) < 5L || uniqueN(d$Age[test]) < 3L) next
    # Negative control preserves each training study's age distribution.
    # It is a predictive control, not a formal permutation p-value.
    permuted <- replicate(repeats, {
      y <- d$Age[train]
      for (s in unique(d$Dataset[train])) {
        ii <- which(d$Dataset[train] == s)
        y[ii] <- sample(y[ii], length(ii), replace = FALSE)
      }
      y
    })
    for (method in names(inputs$centroids)) {
      x <- inputs$centroids[[method]][indices, , drop = FALSE]
      norm <- sqrt(rowSums(x^2))
      stopifnot(all(is.finite(norm)), all(norm > 0))
      distances <- list(
        euclidean = outer(norm[test]^2, norm[train]^2, "+") -
          2 * tcrossprod(x[test, , drop = FALSE], x[train, , drop = FALSE]),
        cosine = 1 - tcrossprod((x / norm)[test, , drop = FALSE], (x / norm)[train, , drop = FALSE])
      )
      # Euclidean matches the recorded reference graph metric. Retain the
      # initial cosine calculation as sensitivity, not a selectively discarded run.
      for (geometry in names(distances)) {
        distance <- distances[[geometry]]
        nn <- t(apply(distance, 1L, function(z) order(z)[seq_len(k)]))
        predicted <- apply(nn, 1L, function(ii) median(d$Age[train][ii]))
        shuffled <- t(apply(nn, 1L, function(ii) apply(permuted[ii, , drop = FALSE], 2L, median)))
        key <- paste(study, method, geometry, sep = ":")
        predictions[[key]] <- data.table(Dataset = study, Global_Donor_ID = d$Global_Donor_ID[test],
          Method = method, Distance = geometry, Age = d$Age[test], Predicted_Age = predicted,
          Median_Baseline = median(d$Age[train]),
          Shuffled_Training_MAE = rowMeans(abs(shuffled - d$Age[test])))
      }
    }
  }
  p <- rbindlist(predictions)
  stopifnot(nrow(p) > 0L, uniqueN(p$Dataset) >= 2L)
  p[, `:=`(Absolute_Error = abs(Age - Predicted_Age),
    Baseline_Error = abs(Age - Median_Baseline))]
  raw <- p[Method == "Raw", .(Global_Donor_ID, Distance, Absolute_Error)]
  p[, Raw_Error := raw$Absolute_Error[match(paste(Global_Donor_ID, Distance), paste(raw$Global_Donor_ID, raw$Distance))]]
  fold_metrics <- p[, .(Donors = .N, Age_Min = min(Age), Age_Max = max(Age),
    MAE = mean(Absolute_Error), Median_Baseline_MAE = mean(Baseline_Error),
    Shuffled_Training_MAE = mean(Shuffled_Training_MAE),
    Spearman = suppressWarnings(cor(Age, Predicted_Age, method = "spearman"))),
    by = .(Dataset, Method, Distance)]
  summaries <- list()
  comparisons <- unique(p[, .(Method, Distance)])
  for (i in seq_len(nrow(comparisons))) {
    method <- comparisons$Method[[i]]
    geometry <- comparisons$Distance[[i]]
    rows <- p[Method == method & Distance == geometry]
    effect <- function(z) c(MAE = mean(z$Absolute_Error),
      Improvement_vs_Raw = mean(z$Raw_Error - z$Absolute_Error),
      Improvement_vs_median = mean(z$Baseline_Error - z$Absolute_Error),
      Improvement_vs_shuffled = mean(z$Shuffled_Training_MAE - z$Absolute_Error))
    groups <- split(seq_len(nrow(rows)), rows$Dataset)
    estimate <- rowMeans(sapply(groups, function(ii) effect(rows[ii])))
    draws <- replicate(repeats, rowMeans(sapply(groups, function(ii) {
      effect(rows[sample(ii, length(ii), replace = TRUE)])
    })))
    summaries[[i]] <- data.table(Method = method, Distance = geometry, Metric = names(estimate),
      Estimate = as.numeric(estimate), CI_Lower = apply(draws, 1L, quantile, 0.025),
      CI_Upper = apply(draws, 1L, quantile, 0.975),
      Donors = nrow(rows), Datasets = length(groups), CellType = celltype,
      BrainRegion = region, Age_Unit = unit)
  }
  fwrite(rbindlist(folds), file.path(out, "heldout_fold_coverage.tsv"), sep = "\t")
  fwrite(p, file.path(out, "age_predictions.tsv"), sep = "\t")
  fwrite(fold_metrics, file.path(out, "age_metrics_by_dataset.tsv"), sep = "\t")
  fwrite(rbindlist(summaries), file.path(out, "age_metrics_summary.tsv"), sep = "\t")
  writeLines(c(
    paste("Population:", celltype, "|", region, "|", unit),
    "Inputs: current formal cell labels and full metadata; >=20 cells per donor-region-class.",
    "Only exact point ages and uniquely aged, known donors; small studies can contribute training donors but are not scored as test cohorts.",
    "GSE178175 is not used in training or testing because source sample identity remains unresolved; no resource cells removed.",
    "Whole-study held-out predictions; >=10 training donors, >=5 interpolation test donors and >=3 test ages.",
    "Fixed 5-nearest donor centroids, median neighbor age, 50 dimensions; no outcome-based tuning.",
    "Euclidean is primary (annotation/seed_assignments/clustering_contract.tsv: NN_Metric); initial cosine retained as sensitivity.",
    "All methods use the same donor groups and cells. No new integration or cell subsampling.",
    "Negative control: 1000 within-training-study donor-age permutations; test ages unchanged; no permutation p-value.",
    "Primary estimates are dataset-equal means of donor errors; CIs use 1000 within-study donor bootstraps, seed 2026.",
    "CIs condition on fitted predictions and the included studies; they do not estimate between-study generalization uncertainty.",
    "This is transductive validation of fixed, age-unsupervised embeddings, not a newly built independent atlas.",
    "Supports age structure only in the tested population/range; not specific developmental gene programs or causality."
  ), file.path(out, "analysis_description.txt"))
  print(rbindlist(summaries))
}
