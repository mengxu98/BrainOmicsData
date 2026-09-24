#!/usr/bin/env Rscript
# Full-cohort LISI summaries in a 50-dimensional or two-dimensional space.
suppressPackageStartupMessages({
  library(data.table)
  library(RcppHNSW)
  library(lisi)
  library(Rcpp)
  library(digest)
})
source("functions/data_paths.R")

args <- commandArgs(TRUE)
stopifnot(length(args) == 4L)
analysis_root <- normalizePath(args[1])
figure_data_dir <- normalizePath(args[2])
method <- args[3]
space <- args[4]
output_dir <- file.path(figure_data_dir, "tables", "full_lisi", space, method)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setDTthreads(4L)

find_unique_file <- function(root, filename) {
  candidates <- list.files(root, recursive = TRUE, full.names = TRUE)
  matches <- candidates[basename(candidates) == filename]
  if (length(matches) != 1L) {
    stop("Expected one ", filename, " under ", root, "; found ", length(matches))
  }
  matches[[1L]]
}
configured_file <- function(variable, root, filename) {
  value <- Sys.getenv(variable, unset = "")
  if (nzchar(value)) normalizePath(value, mustWork = TRUE) else find_unique_file(root, filename)
}

pin <- readLines(file.path(find.package("lisi"), "BRAINOMICS_PINNED_SOURCE"))
stopifnot(
  "numeric_patch=double_precision_probability_normalization_v1" %in% pin,
  "numeric_patch_sha256=6c8d1e25cf62aaa385ef9fff90eeedfb31519768a9853c5d5a116a85837171dd" %in% pin
)
pin_fields <- strsplit(pin, "=", fixed = TRUE)
pin_values <- setNames(
  vapply(pin_fields, function(value) paste(value[-1L], collapse = "="), character(1)),
  vapply(pin_fields, `[[`, character(1), 1L)
)

metadata_file <- configured_file(
  "BRAINOMICS_METADATA_FILE", analysis_root, "metadata_working.rds"
)
label_file <- configured_file(
  "BRAINOMICS_SOURCE_LABEL_FILE", analysis_root, "independent_labels_by_cell.rds"
)
annotation_file <- Sys.getenv(
  "BRAINOMICS_ANNOTATION_TABLE",
  unset = file.path(brainomics_repo_root(), "results", "annotation", "cluster_annotation.tsv")
)
if (!file.exists(annotation_file)) {
  stop("Set BRAINOMICS_ANNOTATION_TABLE to the adopted cluster annotation")
}

metadata <- as.data.table(readRDS(metadata_file))
source_labels <- as.data.table(readRDS(label_file))
index <- match(metadata$Cells, source_labels$Cells)
stopifnot(!anyNA(index), !anyDuplicated(metadata$Cells), nrow(metadata) == 2602031L)
source_labels <- source_labels[index]
metadata[, Source_Full := as.character(source_labels$Source_Coarse)]

# Preserve the author-provided TAC cycling state as its own source category.
tac <- which(metadata$Dataset == "GSE217511" & grepl("^TAC", metadata$CellType_raw))
stopifnot(length(tac) == 1209L, all(metadata$Source_Full[tac] == "Unassigned"))
metadata[tac, Source_Full := "TAC_cycling_state"]
metadata[is.na(Source_Full) | Source_Full == "Unassigned",
         Source_Full := "Unresolved_source_label"]
stopifnot(sum(metadata$Source_Full == "Unresolved_source_label") == 31L)

annotation <- fread(annotation_file)
if (!"CellType" %in% names(annotation) && "Working_CellType" %in% names(annotation)) {
  annotation[, CellType := Working_CellType]
}
stopifnot(all(c("Cluster", "CellType") %in% names(annotation)))
metadata[, Formal_Type := annotation$CellType[match(Cluster, annotation$Cluster)]]
stopifnot(!anyNA(metadata$Formal_Type))

schemes <- c("Source_Full", "Formal_Type", "Dataset")
labels <- lapply(schemes, function(name) as.integer(factor(metadata[[name]])) - 1L)
names(labels) <- schemes
category_counts <- vapply(labels, function(value) uniqueN(value), integer(1))
rm(source_labels, annotation)
gc(FALSE)

cppFunction('
List drop_self(IntegerMatrix ix, NumericMatrix ds, int start) {
  int n = ix.nrow();
  IntegerMatrix oi(n, 89);
  NumericMatrix od(n, 89);
  for (int i = 0; i < n; i++) {
    int q = 0;
    for (int j = 0; j < ix.ncol() && q < 89; j++) {
      if (ix(i, j) == start + i) continue;
      oi(i, q) = ix(i, j);
      od(i, q) = ds(i, j);
      q++;
    }
    if (q != 89) stop("insufficient nonself neighbors");
  }
  return List::create(Named("idx") = oi, Named("dist") = od);
}')

reductions <- c(
  Raw = "pca", RPCA = "integrated.rpca",
  Harmony = "integrated.harmony", scVI = "integrated.scvi"
)
stopifnot(method %in% names(reductions), space %in% c("latent50", "umap2"))
if (space == "umap2") {
  reductions <- c(
    Raw = "umap.unintegrated", RPCA = "umap.rpca",
    Harmony = "umap.harmony", scVI = "umap.scvi"
  )
}
embedding_file <- find_unique_file(
  analysis_root, paste0("embedding_", reductions[[method]], ".rds")
)
embedding <- readRDS(embedding_file)
expected_dimensions <- if (space == "latent50") 50L else 2L
stopifnot(
  identical(rownames(embedding), metadata$Cells),
  ncol(embedding) == expected_dimensions,
  all(is.finite(embedding))
)

message("Building full-cohort neighbors for ", method, " (", space, ")")
neighbor_index <- hnsw_build(
  embedding, distance = "euclidean", M = 24, ef = 200,
  n_threads = 4, random_seed = 2026
)
scores <- matrix(
  NA_real_, nrow(metadata), length(schemes),
  dimnames = list(NULL, schemes)
)
purity <- scores
for (start in seq.int(1L, nrow(metadata), by = 50000L)) {
  end <- min(start + 49999L, nrow(metadata))
  rows <- start:end
  neighbors <- hnsw_search(
    embedding[rows, , drop = FALSE], neighbor_index,
    k = 90L, ef = 400, n_threads = 4
  )
  neighbors <- drop_self(neighbors$idx, neighbors$dist, start)
  neighbor_ids <- t(neighbors$idx) - 1L
  neighbor_distances <- t(neighbors$dist)
  for (scheme in schemes) {
    value <- 1 / lisi:::compute_simpson_index(
      neighbor_distances, neighbor_ids, labels[[scheme]],
      category_counts[[scheme]], 30
    )
    stopifnot(
      all(is.finite(value)), all(value >= 1 - 1e-12),
      all(value <= category_counts[[scheme]] + 1e-12)
    )
    scores[rows, scheme] <- pmax(value, 1)
    purity[rows, scheme] <- rowMeans(
      matrix(labels[[scheme]][neighbors$idx[, 1:30, drop = FALSE]],
             nrow = length(rows)) == labels[[scheme]][rows]
    )
  }
}
stopifnot(!anyNA(scores), nrow(scores) == 2602031L)
saveRDS(
  list(Cells = metadata$Cells, cLISI = scores, Neighbor30_Purity = purity,
       method = method, sampling = FALSE),
  file.path(output_dir, "cell_scores.rds")
)

summaries <- list()
for (scheme in schemes) {
  metadata[, `:=`(cLISI = scores[, scheme], Neighbor30_Purity = purity[, scheme])]
  donor <- metadata[, .(
    Cells = .N, cLISI = mean(cLISI), Neighbor30_Purity = mean(Neighbor30_Purity)
  ), by = .(Dataset, Canonical_Donor_ID)]
  fwrite(donor, file.path(output_dir, paste0(scheme, "_donor_summary.tsv")), sep = "\t")
  fwrite(donor[, .(
    Donors = .N, Cells = sum(Cells), cLISI = mean(cLISI),
    Neighbor30_Purity = mean(Neighbor30_Purity)
  ), by = Dataset], file.path(output_dir, paste0(scheme, "_dataset_summary.tsv")), sep = "\t")
  fwrite(metadata[, .(
    Cells = .N, cLISI = mean(cLISI), Neighbor30_Purity = mean(Neighbor30_Purity)
  ), by = c("Cluster", scheme)],
  file.path(output_dir, paste0(scheme, "_cluster_summary.tsv")), sep = "\t")
  summaries[[scheme]] <- data.table(
    Method = method, Label_Scheme = scheme, Cells = nrow(metadata),
    Label_Categories = category_counts[[scheme]], Mean = mean(scores[, scheme]),
    Median = median(scores[, scheme]), Q25 = unname(quantile(scores[, scheme], .25)),
    Q75 = unname(quantile(scores[, scheme], .75))
  )
}
fwrite(rbindlist(summaries), file.path(output_dir, "summary.tsv"), sep = "\t")
parameters <- data.table(
  Method = method,
  Space = space,
  Cells = nrow(metadata),
  Dimensions = ncol(embedding),
  Sampling = FALSE,
  Perplexity = 30,
  Neighbors_Excluding_Self = 89L,
  HNSW_M = 24L,
  HNSW_Construction_Ef = 200L,
  HNSW_Search_Ef = 400L,
  Threads = 4L,
  Seed = 2026L,
  Metadata_SHA256 = digest(metadata_file, file = TRUE, algo = "sha256"),
  Source_Labels_SHA256 = digest(label_file, file = TRUE, algo = "sha256"),
  Annotation_SHA256 = digest(annotation_file, file = TRUE, algo = "sha256"),
  Embedding_SHA256 = digest(embedding_file, file = TRUE, algo = "sha256"),
  LISI_Source_Commit = unname(pin_values[["source_commit"]]),
  LISI_Numeric_Patch_SHA256 = unname(pin_values[["numeric_patch_sha256"]])
)
fwrite(parameters, file.path(output_dir, "parameters.tsv"), sep = "\t")
message("Full-cohort LISI summaries written for ", method, " (", space, ")")
