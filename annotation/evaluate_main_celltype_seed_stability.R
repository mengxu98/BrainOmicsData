#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

# Retain minority contributions; majority labels alone conceal split/merge events.
stability_metrics <- function(reference, clusters) {
  stopifnot(
    length(reference) == length(clusters), length(reference) > 1L,
    !anyNA(reference), !anyNA(clusters)
  )
  counts <- data.table(Type = as.character(reference), Cluster = as.character(clusters))[
    , .(Overlap = .N),
    by = .(Type, Cluster)
  ]
  counts[, Cells := sum(Overlap), by = Type]
  counts[, Cluster_Cells := sum(Overlap), by = Cluster]
  counts[, `:=`(
    Recall = Overlap / Cells, Purity = Overlap / Cluster_Cells,
    Jaccard = Overlap / (Cells + Cluster_Cells - Overlap)
  )]
  setorder(counts, Cluster, -Overlap, Type)
  mapping <- counts[, .(
    Mapped_CellType = Type[1L], Majority_Cells = Overlap[1L],
    Cluster_Cells = Cluster_Cells[1L]
  ), by = Cluster]
  counts[, Mapped_CellType := mapping$Mapped_CellType[match(Cluster, mapping$Cluster)]]
  predicted <- counts[, .(Predicted_Cells = sum(Overlap)), by = .(Type = Mapped_CellType)]
  summary <- counts[, .(
    Cells = Cells[1L],
    Concordant_Cells = sum(Overlap[Mapped_CellType == Type]),
    Best_Cluster_Jaccard = max(Jaccard), Largest_Fragment_Fraction = max(Recall),
    Occupied_Clusters = .N
  ), by = Type]
  summary[, Predicted_Cells := predicted$Predicted_Cells[match(Type, predicted$Type)]]
  summary[is.na(Predicted_Cells), Predicted_Cells := 0L]
  summary[, `:=`(
    Cell_Concordance = Concordant_Cells / Cells,
    Precision = fifelse(Predicted_Cells > 0, Concordant_Cells / Predicted_Cells, NA_real_),
    F1 = 2 * Concordant_Cells / (Cells + Predicted_Cells)
  )]
  list(counts = counts, mapping = mapping, summary = summary)
}

stability_ari <- function(left, right) {
  stopifnot(length(left) == length(right), length(left) > 1L, !anyNA(left), !anyNA(right))
  choose2 <- function(x) x * (x - 1) / 2
  tab <- table(left, right)
  observed <- sum(choose2(tab))
  a <- sum(choose2(rowSums(tab)))
  b <- sum(choose2(colSums(tab)))
  expected <- a * b / choose2(sum(tab))
  maximum <- (a + b) / 2
  if (maximum == expected) 1 else (observed - expected) / (maximum - expected)
}

run_stability_audit <- function() {
  source("functions/data_paths.R")
  source("functions/processed_object.R")
  source("functions/utils.R")
  setDTthreads(2L)
  started <- Sys.time()
  annotation_dir <- brainomics_data_path("integration_25", "annotation")
  input_dir <- file.path(annotation_dir, "seed_assignments")
  output_dir <- file.path(annotation_dir, "celltype_seed_stability")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  unlink(file.path(output_dir, "_SUCCESS"))
  assignment_file <- file.path(annotation_dir, "celltype_assignments.rds")
  main <- as.data.table(readRDS(assignment_file))
  stopifnot(
    nrow(main) == 2712452L, !anyDuplicated(main$Cells),
    !anyNA(main[, .(Cells, Cluster, CellType)]), uniqueN(main$CellType) == 16L
  )
  primary <- "RPCA_Seed20260730_Resolution2_0"
  directories <- file.path(input_dir, c(".", "resolution_0_8_to_1_2", "resolution_1_5"))
  variants <- data.table(Cells = main$Cells)
  provenance <- list()
  graph_hash <- NULL
  for (directory in directories) {
    path <- file.path(directory, "cluster_assignments.rds")
    contract_path <- file.path(directory, "clustering_contract.tsv")
    contract <- fread(contract_path)
    field <- function(key) {
      value <- contract[Field == key, Value]
      stopifnot(length(value) == 1L, !is.na(value))
      value
    }
    hash <- processed_file_sha256(path)
    stopifnot(hash == field("Assignment_File_SHA256"))
    if (is.null(graph_hash)) graph_hash <- field("Graph_Checkpoint_SHA256")
    stopifnot(graph_hash == field("Graph_Checkpoint_SHA256"))
    x <- as.data.table(readRDS(path))
    stopifnot(identical(x$Cells, main$Cells), !anyDuplicated(x$Cells))
    columns <- grep("^RPCA_Seed[0-9]+_Resolution[0-9]+_[0-9]+$", names(x), value = TRUE)
    stopifnot(length(columns) > 0L)
    for (column in columns) {
      stopifnot(!anyNA(x[[column]]))
      if (column %in% names(variants)) {
        stopifnot(identical(variants[[column]], x[[column]]))
      } else {
        variants[, (column) := x[[column]]]
      }
    }
    provenance[[directory]] <- data.table(
      Input = file.path(basename(directory), basename(path)),
      SHA256 = hash, Contract_SHA256 = processed_file_sha256(contract_path),
      Graph_SHA256 = graph_hash, Assignments = paste(columns, collapse = ";")
    )
  }
  stopifnot(
    primary %in% names(variants),
    identical(paste0("C", sprintf("%02d", as.integer(variants[[primary]]))), as.character(main$Cluster))
  )
  columns <- setdiff(names(variants), "Cells")
  design <- data.table(
    Assignment = columns,
    Seed = as.integer(sub("^RPCA_Seed([0-9]+)_.*$", "\\1", columns)),
    Resolution = as.numeric(sub("_", ".", sub(".*_Resolution", "", columns))),
    Is_Baseline = columns == primary
  )
  design[, Seed_ID := LETTERS[match(Seed, sort(unique(Seed)))]]
  setorder(design, Resolution, Seed)
  design[, Display := paste0(Seed_ID, " / ", Resolution)]
  summaries <- maps <- flows <- cluster_matches <- list()
  for (column in design$Assignment) {
    result <- stability_metrics(main$CellType, variants[[column]])
    result$summary[, Assignment := column]
    result$mapping[, Assignment := column]
    result$counts[, Assignment := column]
    summaries[[column]] <- result$summary
    maps[[column]] <- result$mapping
    flows[[column]] <- result$counts
    overlap <- stability_metrics(main$Cluster, variants[[column]])$counts
    setorder(overlap, Type, -Jaccard, -Overlap, Cluster)
    best <- overlap[, .SD[1L], by = Type]
    best[, Assignment := column]
    cluster_matches[[column]] <- best
  }
  summary <- rbindlist(summaries)
  setnames(summary, "Type", "Main_CellType")
  cluster_overlap <- rbindlist(cluster_matches)
  reference_labels <- unique(main[, .(Cluster, CellType)])
  stopifnot(!anyDuplicated(reference_labels$Cluster))
  cluster_overlap[, Main_CellType := reference_labels$CellType[match(Type, reference_labels$Cluster)]]
  class_overlap <- cluster_overlap[, .(
    Cluster_Weighted_Best_Jaccard =
      weighted.mean(Jaccard, Cells)
  ), by = .(Assignment, Main_CellType)]
  summary <- merge(summary, class_overlap, by = c("Assignment", "Main_CellType"), sort = FALSE)
  overall <- summary[, .(
    Cells = sum(Cells),
    Cell_Weighted_Broad_Class_Concordance = sum(Concordant_Cells) / sum(Cells),
    Macro_Recall = mean(Cell_Concordance), Macro_F1 = mean(F1),
    Minimum_Class_Recall = min(Cell_Concordance)
  ), by = Assignment]
  pairs <- combn(design$Assignment, 2L, simplify = FALSE)
  ari <- rbindlist(lapply(pairs, function(pair) {
    a <- design[Assignment == pair[1L]]
    b <- design[Assignment == pair[2L]]
    data.table(
      Assignment_A = pair[1L], Assignment_B = pair[2L],
      Comparison_Type = if (a$Resolution == b$Resolution) "seed" else if (a$Seed == b$Seed) "resolution" else "seed_and_resolution",
      Adjusted_Rand_Index = stability_ari(variants[[pair[1L]]], variants[[pair[2L]]]), Cells = nrow(main)
    )
  }))
  files <- list(
    "seed_main_celltype_concordance.tsv" = summary,
    "seed_main_celltype_stability_summary.tsv" = overall,
    "seed_cluster_main_celltype_mapping.tsv" = rbindlist(maps),
    "cluster_celltype_overlap.tsv" = rbindlist(flows),
    "primary_cluster_best_overlap.tsv" = cluster_overlap,
    "cluster_pairwise_ari.tsv" = ari, "stability_design.tsv" = design,
    "stability_input_provenance.tsv" = rbindlist(provenance)
  )
  for (name in names(files)) fwrite(files[[name]], file.path(output_dir, name), sep = "\t", quote = FALSE)
  manifest <- data.table(
    File = names(files),
    SHA256 = vapply(file.path(output_dir, names(files)), processed_file_sha256, character(1L)),
    Annotation_SHA256 = processed_file_sha256(assignment_file), Cells = nrow(main)
  )
  fwrite(manifest, file.path(output_dir, "seed_main_celltype_stability_manifest.tsv"), sep = "\t")
  write_annotation_input(output_dir, c(names(files), "seed_main_celltype_stability_manifest.tsv"))
  writeLines(
    c(
      "Status=complete", paste0("Cells=", nrow(main)), paste0("Assignments=", nrow(design)),
      paste0("Elapsed_Seconds=", round(as.numeric(difftime(Sys.time(), started, units = "secs")), 2))
    ),
    file.path(output_dir, "_SUCCESS")
  )
  message("Completed seed/resolution audit with ", nrow(design), " assignments")
}

if (sys.nframe() == 0L) run_stability_audit()
