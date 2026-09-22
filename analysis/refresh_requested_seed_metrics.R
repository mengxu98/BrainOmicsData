#!/usr/bin/env Rscript
# Reaggregate saved partitions against current labels; preserve original run provenance.
source("functions/utils.R")
source("functions/data_paths.R")
source("analysis/evaluate_main_celltype_seed_stability.R")
setDTthreads(2L)
base <- brainomics_data_path("integration_25", "annotation")
main_file <- file.path(base, "celltype_assignments.rds")
main <- as.data.table(read_celltype_assignments())
focus <- main$Cluster == "C64"
stopifnot(sum(focus) == 5251L, all(main$CellType[focus] == "CGE-derived inhibitory neurons"))
old_file <- brainomics_data_path("annotation_history", "c64_excitatory_formal_20260908", "data", "integration_25", "annotation", "celltype_assignments.rds")
old <- readRDS(old_file)
stopifnot(
  identical(main$Cells, old$Cells), identical(main$Cluster, old$Cluster),
  sum(main$CellType != old$CellType) == 21759L
)
out <- file.path(base, "requested_seed_stability")
dir.create(out, showWarnings = FALSE)
rows <- classes <- overlaps <- list()
for (seed in c(2026L, 42L, 11L)) {
  input <- file.path(base, "seed_sensitivity_20260908", paste0("run_", seed))
  path <- file.path(input, paste0("seed_", seed, ".rds"))
  result <- readRDS(path)
  metrics <- stability_metrics(main$CellType, result$Labels)
  pred <- metrics$mapping$Mapped_CellType[match(result$Labels, metrics$mapping$Cluster)]
  fm <- stability_metrics(ifelse(focus, "C64", "Other"), result$Labels)
  stopifnot(!anyNA(pred))
  rows[[as.character(seed)]] <- data.table(
    Seed = seed, Cells = nrow(main), Cell_Types = uniqueN(main$CellType),
    C64_Cells = sum(focus), C64_Label_Retention = mean(pred[focus] == main$CellType[focus]),
    C64_Assigned_Excitatory = sum(pred[focus] == "Excitatory neurons"),
    C64_Assigned_CGE = sum(pred[focus] == "CGE-derived inhibitory neurons"),
    C64_Assigned_Radial_Glia = sum(pred[focus] == "Radial glia"),
    C64_Best_Cluster_Jaccard = fm$summary[Type == "C64", Best_Cluster_Jaccard],
    Whole_Atlas_Label_Retention = mean(pred == main$CellType),
    ARI_To_Annotation_Baseline = stability_ari(main$Cluster, result$Labels)
  )
  classes[[as.character(seed)]] <- copy(metrics$summary)[, Seed := seed]
  focus_overlap <- copy(fm$counts[Type == "C64"])
  setnames(focus_overlap, "Mapped_CellType", "C64_Majority_Class")
  focus_overlap[, Mapped_CellType := metrics$mapping$Mapped_CellType[
    match(Cluster, metrics$mapping$Cluster)
  ]]
  overlaps[[as.character(seed)]] <- focus_overlap[, Seed := seed]
}
fwrite(rbindlist(rows), file.path(out, "summary.tsv"), sep = "\t")
fwrite(rbindlist(classes), file.path(out, "celltype_metrics.tsv"), sep = "\t")
fwrite(rbindlist(overlaps), file.path(out, "c64_cluster_overlap.tsv"), sep = "\t")
write_annotation_input(out, c("summary.tsv", "celltype_metrics.tsv", "c64_cluster_overlap.tsv"))
writeLines("Saved partitions reaggregated; original run files unchanged.", file.path(out, "_SUCCESS"))
print(rbindlist(rows)[, 1:10])
