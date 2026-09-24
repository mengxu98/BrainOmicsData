#!/usr/bin/env Rscript
# Serialize completed per-cell LISI scores without recomputing neighbors.
suppressPackageStartupMessages(library(data.table))
setDTthreads(4L)

args <- commandArgs(TRUE)
stopifnot(length(args) %in% c(3L, 4L))
analysis_root <- normalizePath(args[1])
package_dir <- args[2]

existing <- file.path(package_dir, "validation", "lisi.tsv.gz")
if (file.exists(existing)) {
  header <- fread(cmd = paste("gzip -cd", shQuote(existing)), sep = "\t", nrows = 2L)
  if (ncol(header) == 9L) {
    rows <- as.integer(system2(
      "bash", c("-c", shQuote(paste("gzip -cd", shQuote(existing), "| wc -l"))),
      stdout = TRUE
    ))
    if (!is.na(rows) && rows - 1L == 2602031L) {
      message("Existing per-cell LISI table retained")
      quit(save = "no", status = 0L)
    }
  }
  stop("Existing validation/lisi.tsv.gz does not match the expected 2,602,031 x 9 layout")
}

results_dir <- normalizePath(
  Sys.getenv("BRAINOMICS_FROZEN_DIR", unset = file.path(dirname(analysis_root), "frozen_run")),
  mustWork = TRUE
)
clusters <- readRDS(file.path(results_dir, "cluster_assignments.rds"))
stopifnot(nrow(clusters) == 2602031L, !anyDuplicated(clusters$Cell))

dataset_score_file <- function(directory, method, reduction) {
  candidates <- c(
    file.path(directory, method, paste0(reduction, "_dataset_lisi.rds")),
    file.path(directory, paste0(tolower(method), "_dataset_lisi.rds")),
    file.path(directory, paste0(reduction, "_dataset_lisi.rds"))
  )
  matches <- candidates[file.exists(candidates)]
  if (length(matches) != 1L) {
    stop("Expected one dataset-LISI result for ", method, "; found ", length(matches))
  }
  matches[[1L]]
}
score_values <- function(path) {
  score <- readRDS(path)
  if (is.list(score) && all(c("Cell", "Score") %in% names(score))) {
    return(list(Cell = score$Cell, Score = score$Score))
  }
  if (is.numeric(score) && !is.null(names(score))) {
    return(list(Cell = names(score), Score = unname(score)))
  }
  stop("Unsupported score object: ", path)
}

full_lisi_dir <- Sys.getenv("BRAINOMICS_FULL_LISI_DIR", unset = "")
if (!nzchar(full_lisi_dir)) {
  directories <- list.dirs(analysis_root, recursive = TRUE, full.names = TRUE)
  matches <- directories[file.exists(file.path(
    directories, "tables", "full_lisi", "latent50", "Raw", "cell_scores.rds"
  ))]
  if (length(matches) != 1L) {
    stop("Set BRAINOMICS_FULL_LISI_DIR; found ", length(matches), " candidate directories")
  }
  full_lisi_dir <- file.path(matches[[1L]], "tables", "full_lisi")
}

methods <- c(
  Raw = "pca", scVI = "integrated.scvi",
  Harmony = "integrated.harmony", RPCA = "integrated.rpca"
)
output <- data.table(Cells = sprintf("Cell%07d", seq_len(nrow(clusters))))
if (length(args) != 4L) stop("Supply the dataset-LISI result directory as argument 4")
dataset_lisi_dir <- normalizePath(args[4], mustWork = TRUE)

for (method in names(methods)) {
  dataset_scores <- score_values(dataset_score_file(
    dataset_lisi_dir, method, methods[[method]]
  ))
  index <- match(clusters$Cell, as.character(dataset_scores$Cell))
  if (anyNA(index) || anyDuplicated(dataset_scores$Cell)) {
    stop("Dataset-LISI cells differ from the released cohort for ", method)
  }
  dataset_scores$Score <- dataset_scores$Score[index]

  source_scores <- readRDS(file.path(
    full_lisi_dir, "latent50", method, "cell_scores.rds"
  ))
  stopifnot(
    identical(as.character(source_scores$Cells), as.character(clusters$Cell)),
    length(dataset_scores$Score) == nrow(clusters),
    nrow(source_scores$cLISI) == nrow(clusters),
    all(is.finite(dataset_scores$Score)),
    all(is.finite(source_scores$cLISI[, "Source_Full"]))
  )
  name <- if (method == "Raw") "Raw_PCA" else method
  output[, (paste0(name, "_iLISI")) := dataset_scores$Score]
  output[, (paste0(name, "_cLISI")) := source_scores$cLISI[, "Source_Full"]]
}

setcolorder(
  output,
  c("Cells", paste0(c("Raw_PCA", "scVI", "Harmony", "RPCA"), "_iLISI"),
    paste0(c("Raw_PCA", "scVI", "Harmony", "RPCA"), "_cLISI"))
)
dir.create(file.path(package_dir, "validation"), recursive = TRUE, showWarnings = FALSE)
temporary <- file.path(package_dir, "validation", "lisi.partial")
fwrite(output, temporary, sep = "\t", quote = FALSE, na = "", nThread = 4L)
stopifnot(
  system2("pigz", c("-n", "-f", "-p", "4", shQuote(temporary))) == 0L,
  file.rename(paste0(temporary, ".gz"), existing)
)
message("Per-cell LISI table written")
