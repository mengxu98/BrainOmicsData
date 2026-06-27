#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1]]
}

out_dir <- value_after(
  "--out-dir",
  "../../data/ScienceDB/HumanBrain_sc_snRNAseq_age_resource"
)
fig_dir <- file.path(out_dir, "10_figures_and_tables", "figure_source_data")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(file) utils::read.delim(file, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
write_tsv <- function(x, file) {
  utils::write.table(x, file = file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "NA")
}

source_summary <- read_tsv(file.path(out_dir, "01_source_datasets", "source_dataset_summary.tsv"))
sample_meta <- read_tsv(file.path(out_dir, "02_sample_metadata", "sample_metadata_harmonized.tsv"))
missingness <- read_tsv(file.path(out_dir, "02_sample_metadata", "metadata_missingness_summary.tsv"))
counts <- read_tsv(file.path(out_dir, "03_cell_metadata", "cell_counts_by_age_region_celltype.tsv"))

write_tsv(source_summary, file.path(fig_dir, "figure2_source_dataset_coverage.tsv"))

age_counts <- aggregate(
  number_of_cells_after_filtering ~ age_interval_id + age_interval + age_range,
  sample_meta, sum
)
sample_counts <- aggregate(
  sample_id ~ age_interval_id + age_interval + age_range,
  sample_meta, function(x) length(unique(x))
)
names(sample_counts)[4] <- "sample_count"
age_coverage <- merge(
  age_counts,
  sample_counts,
  by = c("age_interval_id", "age_interval", "age_range")
)
write_tsv(age_coverage, file.path(fig_dir, "figure3_age_interval_coverage.tsv"))

region_cell_counts <- aggregate(
  number_of_cells_after_filtering ~ brain_region_harmonized,
  sample_meta, sum
)
region_sample_counts <- aggregate(
  sample_id ~ brain_region_harmonized,
  sample_meta, function(x) length(unique(x))
)
names(region_sample_counts)[2] <- "sample_count"
region_coverage <- merge(region_cell_counts, region_sample_counts, by = "brain_region_harmonized")
write_tsv(region_coverage, file.path(fig_dir, "figure3_brain_region_coverage.tsv"))

age_region <- aggregate(
  number_of_cells_after_filtering ~ age_interval_id + brain_region_harmonized,
  sample_meta, sum
)
write_tsv(age_region, file.path(fig_dir, "figure3_age_region_heatmap_source.tsv"))

write_tsv(missingness, file.path(fig_dir, "figure4_metadata_missingness_source.tsv"))

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  message("ggplot2 not installed; source data tables were written but PDF figures were skipped.")
  quit(save = "no", status = 0)
}

library(ggplot2)

theme_data_paper <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text = element_text(face = "bold"),
      legend.key.height = unit(3, "mm"),
      legend.key.width = unit(3, "mm")
    )
}

pdf(file.path(fig_dir, "figure2_dataset_cell_counts.pdf"), width = 7.2, height = 4.2, family = "sans")
print(
  ggplot(source_summary, aes(
    x = reorder(source_dataset, number_of_cells_after_filtering),
    y = number_of_cells_after_filtering,
    fill = access_level
  )) +
    geom_col(width = 0.75) +
    coord_flip() +
    scale_y_continuous(labels = function(x) format(x, big.mark = ",")) +
    labs(x = NULL, y = "Cells / nuclei after filtering", fill = "Access") +
    theme_data_paper()
)
dev.off()

age_coverage$stage_order <- as.integer(sub("^S", "", age_coverage$age_interval_id))
age_coverage <- age_coverage[order(age_coverage$stage_order), ]
age_coverage$age_interval_id <- factor(
  age_coverage$age_interval_id,
  levels = age_coverage$age_interval_id
)
pdf(file.path(fig_dir, "figure3_age_interval_counts.pdf"), width = 7.2, height = 2.6, family = "sans")
print(
  ggplot(age_coverage, aes(x = age_interval_id, y = number_of_cells_after_filtering)) +
    geom_col(fill = "#3182BD", width = 0.75) +
    geom_line(aes(y = sample_count * max(number_of_cells_after_filtering) / max(sample_count), group = 1), colour = "#D24B40") +
    geom_point(aes(y = sample_count * max(number_of_cells_after_filtering) / max(sample_count)), colour = "#D24B40", size = 1) +
    scale_y_continuous(
      labels = function(x) format(round(x), big.mark = ","),
      sec.axis = sec_axis(~ . * max(age_coverage$sample_count) / max(age_coverage$number_of_cells_after_filtering),
                          name = "Sample count")
    ) +
    labs(x = "Age interval", y = "Cells / nuclei") +
    theme_data_paper()
)
dev.off()

age_region$age_interval_id <- factor(age_region$age_interval_id, levels = paste0("S", 1:15))
pdf(file.path(fig_dir, "figure3_age_region_heatmap.pdf"), width = 8.0, height = 5.5, family = "sans")
print(
  ggplot(age_region, aes(x = age_interval_id, y = brain_region_harmonized, fill = log10(number_of_cells_after_filtering + 1))) +
    geom_tile(colour = "grey90", linewidth = 0.1) +
    scale_fill_gradient(low = "white", high = "#2166AC", name = "log10(cells + 1)") +
    labs(x = "Age interval", y = "Brain region") +
    theme_data_paper(base_size = 6)
)
dev.off()

missingness$source_dataset <- factor(missingness$source_dataset, levels = unique(source_summary$source_dataset))
missingness$field <- factor(missingness$field, levels = rev(unique(missingness$field)))
pdf(file.path(fig_dir, "figure4_metadata_missingness_heatmap.pdf"), width = 7.2, height = 4.8, family = "sans")
print(
  ggplot(missingness, aes(x = source_dataset, y = field, fill = status)) +
    geom_tile(colour = "white", linewidth = 0.15) +
    scale_fill_manual(values = c(complete = "#2C7BB6", partial = "#FDAE61", missing = "#D7191C")) +
    labs(x = NULL, y = NULL, fill = "Status") +
    theme_data_paper(base_size = 6) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
)
dev.off()

message("Resource figure source tables and draft PDFs written: ", fig_dir)
