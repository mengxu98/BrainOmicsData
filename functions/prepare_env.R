packages <- c(
  "Seurat", "patchwork", "scop", "grid", "ggplot2", "thisplot",
  "thisutils"
)
missing_packages <- packages[!vapply(
  packages,
  requireNamespace,
  logical(1L),
  quietly = TRUE
)]
if (length(missing_packages) > 0L) {
  stop(
    "The frozen BrainOmicsData R environment is incomplete; missing: ",
    paste(missing_packages, collapse = ", "),
    ". Install the version-locked environment before running the workflow."
  )
}

library(grid)
library(ggplot2)
library(patchwork)
library(Seurat)
library(thisutils)
library(scop)

source("functions/utils.R")

check_dir <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    thisutils::log_message(
      "{.path {dir_path}} does not exist. Creating it"
    )
    dir.create(dir_path, recursive = TRUE)
  }
  return(dir_path)
}

retain_source_metadata <- function(metadata, canonical_columns) {
  if (!is.data.frame(metadata)) {
    stop("metadata must be a data.frame")
  }
  missing <- setdiff(canonical_columns, names(metadata))
  if (length(missing) > 0L) {
    stop(
      "canonical metadata columns are missing: ",
      paste(missing, collapse = ", ")
    )
  }
  if (anyNA(metadata$Cells) || any(metadata$Cells == "") ||
    anyDuplicated(metadata$Cells)) {
    stop("metadata Cells must be complete and unique")
  }
  metadata[
    ,
    unique(c(canonical_columns, names(metadata))),
    drop = FALSE
  ]
}

standardize_source_sex_value <- function(value) {
  raw <- trimws(as.character(value))
  normalized <- tolower(raw)
  result <- rep(NA_character_, length(normalized))
  result[!is.na(normalized) & normalized %in% c("f", "female")] <-
    "Female"
  result[!is.na(normalized) & normalized %in% c("m", "male")] <-
    "Male"
  result
}

color_sets <- attr(thisplot::chinese_colors, "color_sets", exact = TRUE)

color_celltypes <- c(
  "Radial glia" = "#8076A3",
  "Glial progenitor cells" = "#A88FBF",
  "Neuroblasts" = "#ED5736",
  "Excitatory neurons" = "#0AA344",
  "Inhibitory neurons" = "#2177B8",
  "Astrocytes" = "#D70440",
  "Oligodendrocyte progenitor cells" = "#F9BD10",
  "Oligodendrocytes" = "#B14B28",
  "Microglia" = "#006D87",
  "Endothelial cells" = "#5E7987",
  "Histaminergic neurons" = "#A24D70",
  "Lymphocytes" = "#7E57C2",
  "Mural cells" = "#8C6D31",
  "Ependymal cells" = "#4DBBD5"
)


color_stages <- c(
  colorRampPalette(
    c("#0AA344", "#006D87")
  )(7),
  colorRampPalette(
    c("#2B73AF", "#003D74")
  )(8)
)
names(color_stages) <- paste0("S", 1:15)
