#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(scop)
  library(SeuratObject)
  library(ggplot2)
})
source("functions/config.R")
setDTthreads(2)
grDevices::pdfFonts(Arial = grDevices::pdfFonts("ArialMT")[[1]])
figure_dir <- "figures"
output_dir <- "results/annotation"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

formal <- as.data.table(read_celltype_assignments())
mapping <- unique(formal[, .(Cluster, CellType)])
type_counts <- formal[, .(Cells = .N), by = .(Group = CellType)]
adopted <- a[, .(Cluster, CellType = Working_CellType, Cells)]
stopifnot(
  nrow(mapping) == 75L,
  !anyDuplicated(mapping$Cluster),
  nrow(type_counts) == 12L,
  sum(type_counts$Cells) == 2602031L,
  type_counts[Group == "Excitatory neurons", Cells] == 725557L,
  identical(
    mapping$CellType[match(adopted$Cluster, mapping$Cluster)],
    adopted$CellType
  )
)
adopted_counts <- adopted[, .(Cells = sum(Cells)), by = CellType]
stopifnot(identical(
  type_counts$Cells[match(adopted_counts$CellType, type_counts$Group)],
  adopted_counts$Cells
))
type_levels <- names(brainomics_celltype_colors)
stopifnot(identical(sort(type_levels), sort(type_counts$Group)))
rm(formal)
gc()

spec <- fread("results/annotation/display_markers.tsv")
stopifnot(
  nrow(spec) == 33L,
  !anyDuplicated(spec$Gene),
  setequal(spec$Marker_Group, type_levels)
)

evidence <- fread("results/annotation/full_cluster_marker_evidence.tsv")
evidence <- evidence[Gene %in% spec$Gene]
stopifnot(
  nrow(evidence) == 75L * nrow(spec),
  all(evidence$Available_Cells > 0),
  !anyDuplicated(evidence[, .(Cluster, Gene)])
)
evidence[, Group := mapping$CellType[match(Cluster, mapping$Cluster)]]
summary <- evidence[,
  .(
    Available_Cells = sum(Available_Cells),
    MeanLog = sum(MeanLog * Available_Cells) / sum(Available_Cells),
    PctPositive = sum(PctPositive * Available_Cells) / sum(Available_Cells)
  ),
  by = .(Group, Gene)
]
summary[, Cells := type_counts$Cells[match(Group, type_counts$Group)]]
stopifnot(
  nrow(summary) == 12L * nrow(spec),
  all(summary$Available_Cells <= summary$Cells)
)
fwrite(summary, file.path(output_dir, "celltype_marker_source.tsv"), sep = "\t")

save_tight_pdf <- function(plot, path, width, height) {
  temporary <- tempfile(fileext = ".pdf")
  on.exit(unlink(temporary), add = TRUE)
  ggsave(temporary, plot, device = grDevices::cairo_pdf, width = width, height = height,
    units = "mm", family = "Arial", bg = "white")
  status <- system2("pdfcrop", c("--margins", "4", shQuote(temporary), shQuote(path)))
  if (!identical(status, 0L)) stop("pdfcrop failed for ", path)
}

# Reserve canvas width for full cell-type row labels.
grDevices::pdf(file = NULL, family = "Arial")
p_flip <- celltype_summary_group_heatmap(as.data.frame(summary), as.data.frame(spec),
  type_levels, flip = TRUE)
grDevices::dev.off()
save_tight_pdf(p_flip$plot, file.path(figure_dir, "fig3b.pdf"), 280, 140)
message("Cell-type marker heatmap completed")
