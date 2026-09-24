#!/usr/bin/env Rscript

# Figure 1 panels A and C: cohort composition and the analysis workflow.
# Panel B is written by plotting/fig1b_age_coverage.R; the figures are assembled
# by functions/assemble_figures.R.
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

source("functions/utils.R")

panel_dir <- "figures"
asset_dir <- "results/annotation"
summary_dir <- Sys.getenv("BRAINOMICS_REFERENCE_SUMMARY", unset = "")
if (!nzchar(summary_dir)) {
  analysis_root <- Sys.getenv("BRAINOMICS_ANALYSIS_DIR", unset = file.path("results", "analysis_run"))
  candidates <- list.files(analysis_root, recursive = TRUE, full.names = TRUE)
  matches <- dirname(candidates[basename(candidates) == "reference_summary.tsv"])
  if (length(matches) != 1L) {
    stop("Set BRAINOMICS_REFERENCE_SUMMARY; found ", length(matches), " candidate directories")
  }
  summary_dir <- matches[[1L]]
}

dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(asset_dir, recursive = TRUE, showWarnings = FALSE)
reference <- read_tsv(file.path(summary_dir, "reference_summary.tsv"))
datasets <- read_tsv(file.path(summary_dir, "dataset_summary.tsv"))
ages <- read_tsv(file.path(summary_dir, "age_interval_summary.tsv"))
reported_ages <- read_tsv(file.path(summary_dir, "reported_age_summary.tsv"))
age_sex_specimens <- read_tsv(file.path(
  summary_dir, "reported_age_sex_specimen_summary.tsv"
))
sex <- read_tsv(file.path(summary_dir, "reported_sex_summary.tsv"))
celltypes <- read_tsv(file.path(summary_dir, "celltype_summary.tsv"))
access <- read_tsv("data/source_access_summary.tsv")
reference_value <- setNames(reference$Value, reference$Field)
expected <- c(
  Reference_Datasets = 22, Cells_or_Nuclei = 2602031,
  Known_Donors = 286, Known_Specimens = 448,
  Verified_Sequence_Libraries = 3744, Age_Intervals = 15,
  Brain_Regions = 48, Main_Cell_Types = 14
)
observed <- suppressWarnings(as.numeric(reference_value[names(expected)]))
if (!identical(unname(observed), unname(as.numeric(expected))) ||
  nrow(datasets) != 22L || sum(datasets$Cells) != 2602031 ||
  sum(ages$Cells) != 2602031 || sum(sex$Cells) != 2602031 ||
  sum(celltypes$Cells) != 2602031 ||
  sum(reported_ages$Cells) != 2602031 ||
  sum(age_sex_specimens$Specimens) != 448L) {
  stop("Figure 1 inputs differ from the reference summary")
}

access <- access[access$dataset %in% datasets$Dataset, , drop = FALSE]
if (nrow(access) != 22L || length(unique(access$dataset)) != 22L) {
  stop("Figure 1 access ledger does not contain exactly the 22 reference datasets")
}

fmt_int <- function(x) {
  format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
}
fig1_theme <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#2B2B2B"),
      axis.ticks = element_line(linewidth = 0.30, colour = "#2B2B2B"),
      axis.text = element_text(colour = "#2B2B2B"),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.7),
      panel.grid.major.y = element_line(linewidth = 0.18, colour = "#E4E7EB"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = base_size + 0.8, face = "bold"),
      plot.subtitle = element_text(size = base_size - 0.2, colour = "#4A4A4A"),
      plot.tag = element_text(size = base_size + 1, face = "bold")
    )
}

count_levels <- function(x, levels) {
  out <- as.data.frame(table(factor(x, levels = levels)), stringsAsFactors = FALSE)
  names(out) <- c("group", "count")
  out[out$count > 0, , drop = FALSE]
}

access$access_group <- ifelse(
  grepl("controlled|mixed", access$access_status, ignore.case = TRUE),
  "Public + controlled", "Public source"
)
access_totals <- count_levels(
  access$access_group, c("Public source", "Public + controlled")
)

profile_columns <- c("Predominant_Technology", "Predominant_Modality")
if (any(!profile_columns %in% names(datasets))) {
  stop("Figure 1 dataset summary lacks predominant technology or modality")
}
datasets$technology_group <- datasets$Predominant_Technology
technology_totals <- count_levels(
  datasets$technology_group,
  c(
    "10x Genomics", "Fluidigm C1", "EasySci-RNA", "snDrop-seq",
    "Smart-seq2 UMI"
  )
)

datasets$modality_group <- datasets$Predominant_Modality
modality_totals <- count_levels(
  datasets$modality_group, c("scRNA-seq", "snRNA-seq")
)

sex_totals <- aggregate(Specimens ~ Reported_Sex, age_sex_specimens, sum)
names(sex_totals) <- c("group", "count")
sex_totals$group <- factor(
  sex_totals$group,
  levels = c("Female", "Male", "Not reported")
)
sex_totals <- sex_totals[order(sex_totals$group), , drop = FALSE]
sex_totals$group <- as.character(sex_totals$group)
if (sum(sex_totals$count) != 448L ||
  !identical(sex_totals$group, c("Female", "Male", "Not reported"))) {
  stop("Figure 1 reported-sex pie does not match the 448-specimen summary")
}

make_pie <- function(df, title, palette) {
  df$fraction <- df$count / sum(df$count)
  pct_label <- ifelse(
    100 * df$fraction < 0.1,
    sprintf("%.2f", 100 * df$fraction),
    sprintf("%.1f", 100 * df$fraction)
  )
  df$legend_label <- paste0(
    df$group, " (", fmt_int(df$count), "; ",
    pct_label, "%)"
  )
  ggplot(df, aes(x = 1, y = count, fill = group)) +
    geom_col(width = 1, colour = "white", linewidth = 0.25) +
    coord_polar(theta = "y") +
    scale_fill_manual(
      values = palette, breaks = df$group,
      labels = df$legend_label, name = NULL
    ) +
    labs(title = title) +
    theme_void(base_family = "Arial") +
    theme(
      plot.title = element_text(size = 7, face = "plain", hjust = 0.5),
      legend.position = "right",
      legend.text = element_text(size = 5.0, margin = margin(l = 0.5)),
      legend.key.spacing.x = grid::unit(0, "mm"),
      legend.key.size = grid::unit(2.6, "mm"),
      legend.spacing.x = grid::unit(0.5, "mm"),
      legend.spacing.y = grid::unit(0.9, "mm"),
      legend.margin = margin(0, 0, 0, -4), plot.margin = margin(1, 1, 1, 1)
    )
}

sex_palette <- c(
  Female = "#B24AA6", Male = "#28327F",
  `Not reported` = "#8A8A8A"
)
p_access <- make_pie(
  access_totals, "Access route\nN = 22 datasets",
  c(`Public source` = "#2CA25F", `Public + controlled` = "#756BB1")
)
p_technology <- make_pie(
  technology_totals, "Technology\nN = 22 datasets",
  c(
    `10x Genomics` = "#FC8D62", `Fluidigm C1` = "#8DA0CB",
    `EasySci-RNA` = "#A6D854", `snDrop-seq` = "#E78AC3",
    `Smart-seq2 UMI` = "#66C2A5"
  )
)
p_modality <- make_pie(
  modality_totals, "Modality\nN = 22 datasets",
  c(`scRNA-seq` = "#6BAED6", `snRNA-seq` = "#2171B5")
)
p_sex <- make_pie(
  sex_totals, "Reported sex\nN = 448 specimens", sex_palette
)
fig1a <- p_access | p_technology | p_modality | p_sex

reported_ages <- reported_ages[
  order(reported_ages$Age_Display_Order), ,
  drop = FALSE
]
ggplot2::ggsave(
  file.path(panel_dir, "fig1a.pdf"), fig1a,
  device = grDevices::cairo_pdf, width = 168, height = 34, units = "mm",
  family = "Arial", bg = "white"
)
write_tsv(reference, file.path(asset_dir, "numeric_summary.tsv"))
write_tsv(datasets, file.path(asset_dir, "figure1_dataset_summary.tsv"))
write_tsv(ages, file.path(asset_dir, "figure1_age_interval_summary.tsv"))
write_tsv(
  reported_ages,
  file.path(asset_dir, "figure1_reported_age_summary.tsv")
)
write_tsv(age_sex_specimens, file.path(
  asset_dir, "figure1_reported_age_sex_specimen_counts.tsv"
))
write_tsv(
  data.frame(
    Reported_Sex = sex_totals$group,
    Specimens = sex_totals$count,
    stringsAsFactors = FALSE
  ),
  file.path(asset_dir, "figure1_reported_sex_summary.tsv")
)
write_tsv(celltypes, file.path(asset_dir, "figure1_celltype_summary.tsv"))
write_tsv(access_totals, file.path(asset_dir, "figure1_access_dataset_totals.tsv"))
write_tsv(
  technology_totals,
  file.path(asset_dir, "figure1_technology_dataset_totals.tsv")
)
write_tsv(
  modality_totals,
  file.path(asset_dir, "figure1_modality_dataset_totals.tsv")
)

message("Figure 1 resource panels written")

source("plotting/fig1b_age_coverage.R", local = new.env())
source("plotting/fig1c_workflow.R", local = new.env())
Sys.setenv(BRAINOMICS_ASSEMBLE = "fig1")
source("functions/assemble_figures.R", local = new.env())
Sys.unsetenv("BRAINOMICS_ASSEMBLE")
source("functions/export_png.R")
export_pdf_png("figures/fig1.pdf", "figures/fig1.png", 6.77)
