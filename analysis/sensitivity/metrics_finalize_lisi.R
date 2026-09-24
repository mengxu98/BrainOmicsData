#!/usr/bin/env Rscript
# Aggregate full-cohort LISI scores by source and compare each method with Raw PCA.
suppressPackageStartupMessages(library(data.table))

figure_data_dir <- commandArgs(TRUE)[1]
input_dir <- file.path(figure_data_dir, "tables", "full_lisi")
methods <- c("Raw", "RPCA", "Harmony", "scVI")
spaces <- c("latent50", "umap2")
schemes <- c("Dataset", "Source_Full", "Formal_Type")
setDTthreads(4L)

dataset_rows <- list()
effect_rows <- list()
cell_rows <- list()
for (space in spaces) {
  for (scheme in schemes) {
    donor_tables <- lapply(methods, function(method) {
      directory <- file.path(input_dir, space, method)
      summary_table <- fread(file.path(directory, "summary.tsv"))
      stopifnot(
        nrow(summary_table) == 3L,
        all(summary_table$Method == method),
        all(summary_table$Cells == 2602031L)
      )
      parameter_file <- file.path(directory, "parameters.tsv")
      if (file.exists(parameter_file)) {
        parameters <- fread(parameter_file)
        stopifnot(
          nrow(parameters) == 1L,
          parameters$Method == method,
          parameters$Space == space,
          parameters$Cells == 2602031L
        )
      }
      table <- fread(file.path(directory, paste0(scheme, "_donor_summary.tsv")))
      table[, Method := method]
      table
    })
    names(donor_tables) <- methods

    for (method in methods) {
      values <- merge(
        donor_tables[[method]],
        donor_tables$Raw[, .(Dataset, Canonical_Donor_ID, Raw_Score = cLISI)],
        by = c("Dataset", "Canonical_Donor_ID")
      )
      values[, Change_From_Raw := cLISI - Raw_Score]
      by_source <- values[, .(
        Cells = sum(Cells), Donors = .N, Mean = mean(cLISI),
        Change_From_Raw = mean(Change_From_Raw),
        Neighbor30_Purity = mean(Neighbor30_Purity)
      ), by = Dataset]
      by_source[, `:=`(Method = method, Space = space, Label_Scheme = scheme)]
      dataset_rows[[paste(space, scheme, method)]] <- by_source

      sensitivity_sets <- list(
        all_studies = by_source,
        exclude_Ma = by_source[Dataset != "Ma_et_al_2022"],
        exclude_three_linked = by_source[
          !Dataset %in% c("Li_et_al_2018", "Ma_et_al_2022", "GSE186538")
        ]
      )
      for (sensitivity in names(sensitivity_sets)) {
        selected <- sensitivity_sets[[sensitivity]]
        margin <- qt(.975, nrow(selected) - 1L) *
          sd(selected$Change_From_Raw) / sqrt(nrow(selected))
        effect_rows[[paste(space, scheme, method, sensitivity)]] <- data.table(
          Space = space, Label_Scheme = scheme, Method = method,
          Sensitivity = sensitivity, Studies = nrow(selected),
          Study_Equal_Mean = mean(selected$Mean),
          Change_From_Raw = mean(selected$Change_From_Raw),
          CI95_Lower = mean(selected$Change_From_Raw) - margin,
          CI95_Upper = mean(selected$Change_From_Raw) + margin
        )
      }
    }
  }
  for (method in methods) {
    summary <- fread(file.path(input_dir, space, method, "summary.tsv"))
    summary[, Space := space]
    cell_rows[[paste(space, method)]] <- summary
  }
}

fwrite(rbindlist(dataset_rows),
       file.path(figure_data_dir, "tables", "full_lisi_dataset_summary.tsv"), sep = "\t")
fwrite(rbindlist(effect_rows),
       file.path(figure_data_dir, "tables", "full_lisi_study_effects.tsv"), sep = "\t")
fwrite(rbindlist(cell_rows),
       file.path(figure_data_dir, "tables", "full_lisi_cell_summary.tsv"), sep = "\t")
message("Full-cohort LISI aggregate tables written")
