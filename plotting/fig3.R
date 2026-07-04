source("functions/prepare_env.R")
source("functions/metadata_schema.R")
source("functions/sample_schema.R")

data_dir <- "../../data/BrainOmicsData/integration/"
fig_dir <- check_dir("figures/")

objects_plot <- readRDS(
  file.path(data_dir, "objects_celltype_plot.rds")
)
objects_plot@meta.data <- add_metadata_schema(objects_plot@meta.data)
objects_plot@meta.data <- add_sample_schema(objects_plot@meta.data)

development_stages <- c(
  "Embryonic",
  "Early fetal",
  "Early fetal",
  "Early mid-fetal",
  "Early mid-fetal",
  "Late mid-fetal",
  "Late fetal",
  "Neonatal and early infancy",
  "Late infancy",
  "Early childhood",
  "Middle and late childhood",
  "Adolescence",
  "Young adulthood",
  "Middle adulthood",
  "Late adulthood"
)

objects_plot$AgeIntervalID <- factor(
  objects_plot$AgeIntervalID,
  levels = paste0("S", 1:15)
)
objects_plot$AgeInterval <- factor(
  objects_plot$AgeInterval,
  levels = unique(development_stages)
)

stage_def <- data.frame(
  AgeIntervalID = paste0("S", 1:15),
  AgeInterval = development_stages,
  AgeRange = c(
    "4-8 PCW",
    "8-10 PCW",
    "10-13 PCW",
    "13-16 PCW",
    "16-19 PCW",
    "19-24 PCW",
    "24-38 PCW",
    "0-0.5 years",
    "0.5-1 years",
    "1-6 years",
    "6-12 years",
    "12-20 years",
    "20-40 years",
    "40-60 years",
    "60+ years"
  ),
  stringsAsFactors = FALSE
)
stage_counts <- as.data.frame(
  table(AgeIntervalID = objects_plot$AgeIntervalID),
  stringsAsFactors = FALSE
)
names(stage_counts)[2] <- "N_cells"
stage_sample_counts <- aggregate(
  Sample_ID ~ AgeIntervalID,
  objects_plot@meta.data,
  function(x) length(unique(x))
)
names(stage_sample_counts)[2] <- "N_samples"
stage_def <- merge(stage_def, stage_counts, by = "AgeIntervalID", all.x = TRUE)
stage_def <- merge(
  stage_def,
  stage_sample_counts,
  by = "AgeIntervalID",
  all.x = TRUE
)
stage_def$N_cells[is.na(stage_def$N_cells)] <- 0
stage_def$N_samples[is.na(stage_def$N_samples)] <- 0
stage_def$AgeIntervalID <- factor(
  stage_def$AgeIntervalID,
  levels = paste0("S", 1:15)
)
stage_def <- stage_def[order(stage_def$AgeIntervalID), ]
stage_def$AgeIntervalID <- as.character(stage_def$AgeIntervalID)

stage_labels <- as.matrix(
  stage_def[, c(
    "AgeIntervalID",
    "AgeInterval",
    "AgeRange",
    "N_cells",
    "N_samples"
  )]
)
stage_labels[, "N_cells"] <- format(
  stage_def$N_cells,
  big.mark = ",",
  trim = TRUE
)
stage_labels[, "N_samples"] <- format(
  stage_def$N_samples,
  big.mark = ",",
  trim = TRUE
)
mat_text <- matrix(0, nrow = 15, ncol = 5)
colnames(mat_text) <- c(
  "Age interval ID",
  "Age interval",
  "Age range",
  "Cell count",
  "Biological samples"
)
rownames(mat_text) <- as.character(stage_def$AgeIntervalID)

stage_table_gp <- gpar(fontsize = 9)
stage_table_padding_mm <- 1.5
stage_cell_fill <- "#EFF7FC"
stage_cell_border <- "#9fbbd3ff"
stage_text_width_mm <- function(labels) {
  label_widths <- lapply(
    labels,
    function(label) grobWidth(textGrob(label, gp = stage_table_gp))
  )
  convertWidth(max(do.call(unit.c, label_widths)), "mm", TRUE)
}
stage_column_widths_mm <- vapply(
  seq_len(ncol(stage_labels)),
  function(j) {
    stage_text_width_mm(c(colnames(mat_text)[j], stage_labels[, j])) +
      stage_table_padding_mm * 2
  },
  numeric(1)
)
stage_color_width_mm <- 5
stage_row_height_mm <- 5.5
stage_header_height_mm <- 5
stage_table_width_mm <- sum(stage_column_widths_mm)
stage_plot_width_mm <- stage_color_width_mm + stage_table_width_mm
stage_plot_height_mm <- stage_row_height_mm *
  nrow(stage_labels) +
  stage_header_height_mm
stage_pdf_width <- (stage_plot_width_mm + 4) / 25.4
stage_pdf_height <- (stage_plot_height_mm + 2) / 25.4

draw_stage_table <- function() {
  grid.newpage()
  pushViewport(
    viewport(
      x = unit(2, "mm"),
      y = unit(1, "mm"),
      width = unit(stage_plot_width_mm, "mm"),
      height = unit(stage_plot_height_mm, "mm"),
      just = c("left", "bottom")
    )
  )

  n_rows <- nrow(stage_labels)
  table_left_mm <- stage_color_width_mm
  table_bottom_mm <- stage_header_height_mm

  grid.rect(
    x = unit(0, "mm"),
    y = unit(
      table_bottom_mm + stage_row_height_mm * (n_rows - seq_len(n_rows)),
      "mm"
    ),
    width = unit(stage_color_width_mm, "mm"),
    height = unit(stage_row_height_mm, "mm"),
    gp = gpar(
      fill = color_stages[as.character(stage_def$AgeIntervalID)],
      col = NA
    ),
    just = c("left", "bottom")
  )

  col_left_mm <- table_left_mm
  for (j in seq_len(ncol(stage_labels))) {
    for (i in seq_len(n_rows)) {
      grid.rect(
        x = unit(col_left_mm, "mm"),
        y = unit(
          table_bottom_mm + stage_row_height_mm * (n_rows - i),
          "mm"
        ),
        width = unit(stage_column_widths_mm[j], "mm"),
        height = unit(stage_row_height_mm, "mm"),
        gp = gpar(fill = stage_cell_fill, col = stage_cell_border, lwd = 0.5),
        just = c("left", "bottom")
      )
      grid.text(
        stage_labels[i, j],
        x = unit(col_left_mm + stage_table_padding_mm, "mm"),
        y = unit(
          table_bottom_mm + stage_row_height_mm * (n_rows - i + 0.5),
          "mm"
        ),
        gp = stage_table_gp,
        just = "left"
      )
    }
    grid.text(
      colnames(mat_text)[j],
      x = unit(col_left_mm + stage_column_widths_mm[j] / 2, "mm"),
      y = unit(stage_header_height_mm / 2, "mm"),
      gp = stage_table_gp
    )
    col_left_mm <- col_left_mm + stage_column_widths_mm[j]
  }

  popViewport()
}

pdf(
  file.path(fig_dir, "development_stage_annotation.pdf"),
  width = stage_pdf_width,
  height = stage_pdf_height
)
draw_stage_table()
dev.off()

p_brain_region <- CellDimPlot(
  objects_plot,
  reduction = "umap.rpca",
  group.by = "BrainRegion",
  palette = "simpsons",
  label = FALSE,
  raster = TRUE,
  xlab = "UMAP_1",
  ylab = "UMAP_2",
  theme_use = "theme_blank_axis"
)

color_stages2 <- color_stages[!duplicated(development_stages)]
names(color_stages2) <- unique(development_stages)
p_age_interval <- CellDimPlot(
  objects_plot,
  reduction = "umap.rpca",
  group.by = "AgeInterval",
  palcolor = color_stages2,
  label = FALSE,
  raster = TRUE,
  xlab = "UMAP_1",
  ylab = "UMAP_2",
  theme_use = "theme_blank_axis"
)
p_brain_region_stage <- p_brain_region / p_age_interval
ggsave(
  file.path(fig_dir, "brainregion_stage.pdf"),
  p_brain_region_stage,
  width = 12,
  height = 7
)
