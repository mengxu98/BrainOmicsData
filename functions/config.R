# Shared configuration for the manuscript figure scripts.
#
# BRAINOMICS_ANALYSIS_DIR selects the analysis run; the default is the run that
# produced the released figures.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
setDTthreads(2L)
root <- Sys.getenv(
  "BRAINOMICS_ANALYSIS_DIR",
  unset = file.path(
    "results", "analysis_run"
  )
)
run <- file.path(root, "07_downstream", "revision_20260918")
doc <- file.path(run, "25_formal_manuscript_review_20260918")
source("functions/utils.R")

a <- fread(file.path(run, "18_final_annotation_20260918", "cluster_annotation.tsv"))
brainomics_celltype_colors <- setNames(
  a$Colour[!duplicated(a$Working_CellType)],
  a$Working_CellType[!duplicated(a$Working_CellType)]
)
type_levels <- c(
  "Astrocytes", "Endothelial cells", "Mural cells", "Fibroblasts",
  "Excitatory neurons", "Inhibitory neurons", "Microglia", "Lymphocytes",
  "Neural progenitors", "Oligodendrocyte progenitor cells",
  "Differentiating oligodendrocytes", "Oligodendrocytes"
)
brainomics_celltype_colors <- brainomics_celltype_colors[type_levels]
stopifnot(
  nrow(a) == 75L,
  sum(a$Cells) == 2602031L,
  !anyNA(brainomics_celltype_colors)
)
