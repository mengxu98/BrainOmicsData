#!/usr/bin/env Rscript
# Prepare Figure 5 source tables from the frozen study-equal summaries.
suppressPackageStartupMessages(library(data.table))
source("functions/config.R")

input <- file.path(doc, "tables/gene_reuse")
output <- "figures"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
genes <- c("PPP4R2", "GXYLT2", "KCNJ3", "SMCHD1", "CTSD", "MRPL23",
           "FHIT", "SLC10A7", "KLHDC4", "DACT1", "DAAM1")

expression <- fread(file.path(input, "expression_study_equal.tsv"))
setnames(expression, "Symbol", "Gene")
stopifnot(setequal(expression$Gene, genes),
          !anyDuplicated(expression[, .(CellType, Gene)]))
fwrite(expression, file.path(output, "fig5_expression_heatmap_source.tsv"), sep = "\t")

study <- fread(file.path(input, "effects_by_study.tsv"))
pooled <- fread(file.path(input, "effects_study_equal.tsv"))
setnames(study, c("Symbol", "Effect"), c("Gene", "Estimate"))
setnames(pooled, c("Symbol", "Effect"), c("Gene", "Estimate"))
pooled[, Dataset := "Study-equal"]
effects <- rbindlist(list(study, pooled), fill = TRUE)
stopifnot(setequal(unique(effects$Gene), genes),
          !anyDuplicated(effects[, .(Dataset, Gene)]))
fwrite(effects, file.path(output, "fig5_oligo_microglia_effect_source.tsv"), sep = "\t")
