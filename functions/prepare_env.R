if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak")
}

if (!requireNamespace("thisutils", quietly = TRUE)) {
  pak::pak("thisutils")
}

packages <- c(
  "Seurat", "patchwork", "mengxu98/scop", "grid", "ggplot2"
)
thisutils::check_r(packages)

log_message <- thisutils::log_message

library(grid)
library(ggplot2)
library(patchwork)
library(Seurat)
library(scop)

check_dir <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    log_message(
      "{.path {dir_path}} does not exist. Creating it"
    )
    dir.create(dir_path, recursive = TRUE)
  }
  return(dir_path)
}

color_sets <- attr(thisplot::chinese_colors, "color_sets", exact = TRUE)

colors32 <- color_sets$ChineseSet32
colors128 <- color_sets$ChineseSet128

color_celltypes <- c(
  "Radial glia" = "#8076A3",
  "Neuroblasts" = "#ED5736",
  "Excitatory neurons" = "#0AA344",
  "Inhibitory neurons" = "#2177B8",
  "Astrocytes" = "#D70440",
  "Oligodendrocyte progenitor cells" = "#F9BD10",
  "Oligodendrocytes" = "#B14B28",
  "Microglia" = "#006D87",
  "Endothelial cells" = "#5E7987"
)

color_stages1 <- colorRampPalette(
  c("#0AA344", "#006D87")
)(7)
color_stages2 <- colorRampPalette(
  c("#2B73AF", "#003D74")
)(8)
color_stages <- c(color_stages1, color_stages2)
names(color_stages) <- paste0("S", 1:15)
