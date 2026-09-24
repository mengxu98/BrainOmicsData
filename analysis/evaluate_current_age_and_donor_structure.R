#!/usr/bin/env Rscript
# Held-out age evaluation and donor-geometry comparison.
suppressPackageStartupMessages({library(data.table);library(Matrix);library(SeuratObject)})
setDTthreads(4L)

args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
analysis_root <- normalizePath(args[1])
figure_data_dir <- normalizePath(args[2])
output_dir <- file.path(figure_data_dir, "tables", "age_signal")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

find_unique_file <- function(root, filename) {
  candidates <- list.files(root, recursive = TRUE, full.names = TRUE)
  matches <- candidates[basename(candidates) == filename]
  if (length(matches) != 1L) {
    stop("Expected one ", filename, " under ", root, "; found ", length(matches))
  }
  matches[[1L]]
}

metadata_file <- Sys.getenv("BRAINOMICS_METADATA_FILE", unset = "")
if (!nzchar(metadata_file)) {
  metadata_file <- find_unique_file(analysis_root, "metadata_working.rds")
}
metadata <- as.data.table(readRDS(metadata_file))

annotation_file <- Sys.getenv(
  "BRAINOMICS_ANNOTATION_TABLE",
  unset = file.path("results", "annotation", "cluster_annotation.tsv")
)
if (!file.exists(annotation_file)) {
  stop("Set BRAINOMICS_ANNOTATION_TABLE to the adopted cluster annotation")
}
annotation <- fread(annotation_file)
annotation_type <- if ("CellType" %in% names(annotation)) {
  "CellType"
} else if ("Working_CellType" %in% names(annotation)) {
  "Working_CellType"
} else {
  stop("The adopted annotation lacks a cell-type column")
}
stopifnot(
  all(c("Cluster", "Cells") %in% names(annotation)),
  nrow(annotation) == 75L,
  !anyDuplicated(annotation$Cluster),
  sum(annotation$Cells) == nrow(metadata),
  "Cluster" %in% names(metadata)
)
adopted_type <- annotation[[annotation_type]][
  match(as.character(metadata$Cluster), as.character(annotation$Cluster))
]
if (anyNA(adopted_type)) {
  stop("The adopted annotation does not cover every metadata cluster")
}
metadata[, Working_CellType := adopted_type]
eligible <-
  !is.na(metadata$Global_Donor_ID) & nzchar(metadata$Global_Donor_ID) &
  !grepl("unknown", metadata$Donor_ID_Verification_Status, ignore.case = TRUE) &
  is.finite(metadata$Age_num) & is.finite(metadata$Age_Lower) &
  is.finite(metadata$Age_Upper) & metadata$Age_Lower == metadata$Age_Upper &
  metadata$Age_num == metadata$Age_Lower &
  !grepl("mean|midpoint|range",
         paste(metadata$Age_Representative_Method, metadata$Age), ignore.case = TRUE) &
  !metadata$Dataset %in% c("GSE178175", "GSE144136")
if ("Continuous_Age_Eligibility" %in% names(metadata)) {
  eligible <- eligible & !grepl(
    "not eligible|mean|midpoint|range",
    metadata$Continuous_Age_Eligibility, ignore.case = TRUE
  )
}
eligible[is.na(eligible)] <- FALSE
metadata[, eligible := eligible]

conflicts <- metadata[eligible == TRUE,
  .(Ages = uniqueN(paste(Age_num, Unit))), by = Canonical_Donor_ID]
metadata[Canonical_Donor_ID %in% conflicts[Ages != 1L, Canonical_Donor_ID],
         eligible := FALSE]
fwrite(conflicts[Ages != 1L], file.path(output_dir, "donor_age_conflicts.tsv"), sep = "\t")
fwrite(metadata[, .(Cells = .N), by = .(Dataset, Unit, eligible)],
       file.path(output_dir, "age_eligibility.tsv"), sep = "\t")

# Restrict the age evaluation to frontal-cortex excitatory neurons with exact ages.
selected <- which(
  metadata$eligible & metadata$Working_CellType == "Excitatory neurons" &
  metadata$BrainRegion == "Frontal cortex" & metadata$Unit == "Years"
)
age_cells <- metadata[selected]
age_cells[, group := paste(Dataset, Canonical_Donor_ID, sep = "|")]
donors <- age_cells[, .(
  Dataset = Dataset[1L], Canonical_Donor_ID = Canonical_Donor_ID[1L],
  Cells = .N, Age = Age_num[1L]
), by = group][Cells >= 20L]
age_cells <- age_cells[group %in% donors$group]
selected <- match(age_cells$Cells, metadata$Cells)
stopifnot(!anyNA(selected))
age_membership <- sparseMatrix(
  i = seq_along(selected), j = match(age_cells$group, donors$group),
  x = 1 / donors$Cells[match(age_cells$group, donors$group)],
  dims = c(length(selected), nrow(donors))
)
fwrite(donors, file.path(output_dir, "donor_metadata.tsv"), sep = "\t")

metadata[, All_Donor_Group := paste(Dataset, Canonical_Donor_ID, sep = "|")]
all_donors <- metadata[, .(
  Dataset = Dataset[1L], Canonical_Donor_ID = Canonical_Donor_ID[1L], Cells = .N
), by = All_Donor_Group]
all_membership <- sparseMatrix(
  i = seq_len(nrow(metadata)), j = match(metadata$All_Donor_Group, all_donors$All_Donor_Group),
  x = 1 / all_donors$Cells[match(metadata$All_Donor_Group, all_donors$All_Donor_Group)],
  dims = c(nrow(metadata), nrow(all_donors))
)

reductions <- c(
  Raw = "pca", RPCA = "integrated.rpca",
  Harmony = "integrated.harmony", scVI = "integrated.scvi"
)
age_centroids <- list()
all_centroids <- list()
for (method in names(reductions)) {
  embedding_file <- find_unique_file(
    analysis_root, paste0("embedding_", reductions[[method]], ".rds")
  )
  embedding <- readRDS(embedding_file)
  stopifnot(identical(metadata$Cells, rownames(embedding)), ncol(embedding) == 50L)
  age_centroids[[method]] <- as.matrix(crossprod(
    age_membership, embedding[selected, , drop = FALSE]
  ))
  all_centroids[[method]] <- as.matrix(crossprod(all_membership, embedding))
  rm(embedding)
  gc(FALSE)
}
rm(age_membership, all_membership)
gc(FALSE)
saveRDS(
  list(donors = donors, centroids = age_centroids,
       all_donors = all_donors, all_centroids = all_centroids),
  file.path(output_dir, "donor_centroids.rds")
)

prediction_rows <- list()
fold_rows <- list()
for (study in unique(donors$Dataset)) {
  all_test <- which(donors$Dataset == study)
  train <- which(
    donors$Dataset != study &
    !donors$Canonical_Donor_ID %in% donors$Canonical_Donor_ID[all_test]
  )
  test <- if (length(train)) {
    all_test[
      donors$Age[all_test] >= min(donors$Age[train]) &
      donors$Age[all_test] <= max(donors$Age[train])
    ]
  } else integer()
  included <- length(train) >= 10L && length(test) >= 5L &&
    uniqueN(donors$Age[test]) >= 3L
  fold_rows[[study]] <- data.table(
    Dataset = study, Train_Donors = length(train), Test_Donors = length(test),
    Outside_Training_Range = length(all_test) - length(test), Included = included
  )
  if (!included) next

  for (method in names(age_centroids)) {
    embedding <- age_centroids[[method]]
    norms <- sqrt(rowSums(embedding^2))
    for (distance in c("euclidean", "cosine")) {
      distance_matrix <- if (distance == "euclidean") {
        outer(norms[test]^2, norms[train]^2, "+") -
          2 * tcrossprod(embedding[test, , drop = FALSE],
                         embedding[train, , drop = FALSE])
      } else {
        1 - tcrossprod(
          (embedding / norms)[test, , drop = FALSE],
          (embedding / norms)[train, , drop = FALSE]
        )
      }
      predicted <- apply(distance_matrix, 1L, function(value) {
        median(donors$Age[train][order(value)[1:5]])
      })
      prediction_rows[[paste(study, method, distance)]] <- data.table(
        Dataset = study, Donor = donors$Canonical_Donor_ID[test],
        Method = method, Distance = distance, Age = donors$Age[test],
        Predicted_Age = predicted, Median_Baseline = median(donors$Age[train])
      )
    }
  }
}
fwrite(rbindlist(fold_rows), file.path(output_dir, "heldout_fold_coverage.tsv"), sep = "\t")
predictions <- rbindlist(prediction_rows)
if (nrow(predictions)) {
  predictions[, `:=`(
    Absolute_Error = abs(Age - Predicted_Age),
    Baseline_Error = abs(Age - Median_Baseline)
  )]
  fwrite(predictions, file.path(output_dir, "age_predictions.tsv"), sep = "\t")
  by_dataset <- predictions[, .(
    Donors = .N, MAE = mean(Absolute_Error),
    Baseline_MAE = mean(Baseline_Error),
    Spearman = suppressWarnings(cor(Age, Predicted_Age, method = "spearman"))
  ), by = .(Dataset, Method, Distance)]
  fwrite(by_dataset, file.path(output_dir, "age_metrics_by_dataset.tsv"), sep = "\t")
  fwrite(by_dataset[, .(
    Datasets = .N, Donors = sum(Donors), Dataset_Equal_MAE = mean(MAE),
    Dataset_Equal_Baseline_MAE = mean(Baseline_MAE)
  ), by = .(Method, Distance)],
  file.path(output_dir, "age_metrics_summary.tsv"), sep = "\t")
}

structure_rows <- list()
for (study in unique(all_donors$Dataset)) {
  rows <- which(all_donors$Dataset == study)
  if (length(rows) < 4L) next
  raw_distances <- as.vector(dist(all_centroids$Raw[rows, , drop = FALSE]))
  for (method in names(all_centroids)) {
    structure_rows[[paste(study, method)]] <- data.table(
      Dataset = study, Donors = length(rows), Method = method,
      Spearman_Donor_Distance_to_Raw = suppressWarnings(cor(
        raw_distances,
        as.vector(dist(all_centroids[[method]][rows, , drop = FALSE])),
        method = "spearman"
      ))
    )
  }
}
fwrite(rbindlist(structure_rows),
       file.path(output_dir, "donor_structure_by_dataset.tsv"), sep = "\t")
message("Age evaluation and donor-geometry summaries written")
