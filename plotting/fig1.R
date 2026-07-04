#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

repo_dir <- normalizePath(".", mustWork = TRUE)
setwd(repo_dir)

out_dir <- file.path("manuscript", "revised_assets")
fig_dir <- file.path("figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

metadata_path <- "../../data/BrainOmicsData/integration/metadata_filtered.rds"
object_plot_path <- "../../data/BrainOmicsData/integration/objects_celltype_plot.rds"
lisi_path <- "results/lisi_results.rds"

read_tsv <- function(path) {
  utils::read.delim(
    path,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

write_tsv <- function(x, path) {
  utils::write.table(
    x,
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
}

fmt_int <- function(x) format(as.integer(round(x)), big.mark = ",", trim = TRUE)

theme_paper <- function(base_size = 7, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#2B2B2B"),
      axis.ticks = element_line(linewidth = 0.3, colour = "#2B2B2B"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.5, colour = "#2B2B2B"),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.7),
      strip.text = element_text(size = base_size - 0.2, face = "bold"),
      panel.grid.major.y = element_line(linewidth = 0.18, colour = "#E4E7EB"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(
        size = base_size + 0.8,
        face = "bold",
        hjust = 0
      ),
      plot.subtitle = element_text(size = base_size - 0.2, colour = "#4A4A4A"),
      plot.tag = element_text(size = base_size + 2, face = "bold")
    )
}

if (!file.exists(metadata_path)) {
  stop("Missing metadata: ", metadata_path)
}
metadata <- readRDS(metadata_path)
metadata[] <- lapply(metadata, as.character)
source("functions/metadata_schema.R")
source("functions/sample_schema.R")
metadata <- add_metadata_schema(metadata)
metadata <- add_sample_schema(metadata)

source_info <- data.frame(
  source_dataset = c(
    "AllenM1",
    "EGAD00001006049",
    "EGAS00001006537",
    "GSE103723",
    "GSE104276",
    "GSE144136",
    "GSE168408",
    "GSE178175",
    "GSE186538",
    "GSE199762",
    "GSE202210",
    "GSE204683",
    "GSE207334",
    "GSE212606",
    "GSE217511",
    "GSE261983",
    "GSE296073",
    "GSE67835",
    "GSE81475",
    "GSE97942",
    "HYPOMAP",
    "Li_et_al_2018",
    "Ma_et_al_2022",
    "Nowakowski_et_al_2017",
    "PRJCA015229",
    "ROSMAP",
    "SomaMut"
  ),
  source_repository = c(
    "Allen Brain Map",
    "EGA",
    "EGA",
    rep("GEO", 17),
    "CELLxGENE",
    "publication supplement",
    "Sestan lab / BrainSCOPE",
    "publication supplement",
    "CNGBdb",
    "ROSMAP / AD Knowledge Portal",
    "project website"
  ),
  source_accession = c(
    "AllenM1",
    "EGAD00001006049",
    "EGAS00001006537",
    "GSE103723",
    "GSE104276",
    "GSE144136",
    "GSE168408",
    "GSE178175",
    "GSE186538",
    "GSE199762; phs003509.v1.p1",
    "GSE202210",
    "GSE204683",
    "GSE207334",
    "GSE212606",
    "GSE217511",
    "GSE261983",
    "GSE296073",
    "GSE67835",
    "GSE81475",
    "GSE97942",
    "HYPOMAP",
    "Li et al. 2018",
    "Ma et al. 2022",
    "Nowakowski et al. 2017",
    "PRJCA015229",
    "ROSMAP",
    "SomaMut"
  ),
  publication_doi = c(
    "10.1038/s41586-021-03465-8",
    "10.1126/science.adf1226",
    "10.1016/j.biopsych.2022.06.033",
    "10.1038/s41422-018-0053-3",
    "10.1038/nature25980",
    "10.1038/s41593-020-0621-y",
    "10.1016/j.cell.2022.09.039",
    "10.1038/s41587-022-01231-3",
    "10.1016/j.neuron.2021.10.036",
    "10.1038/s41586-023-06981-x",
    "10.1126/scitranslmed.abo1997",
    "10.1126/sciadv.adg3754",
    "10.1126/science.abo7257",
    "10.1016/j.cell.2023.08.042",
    "10.1038/s41467-022-34975-2",
    "10.1126/science.adi5199",
    "10.1038/s41586-025-09362-8",
    "10.1073/pnas.1507125112",
    "10.1016/j.celrep.2016.08.038",
    "10.1038/nbt.4038",
    "10.1038/s41586-024-08504-8",
    "10.1126/science.aat7615",
    "10.1126/science.abo7257",
    "10.1126/science.aap8809",
    "10.1016/j.xgen.2024.100703",
    "10.1016/j.cell.2023.08.039",
    "10.1038/s41586-025-09435-8"
  ),
  access_level = c(
    "public",
    "controlled",
    "controlled",
    "public",
    "public",
    "public",
    "public",
    "public",
    "public",
    "mixed_public_controlled",
    "public",
    "public",
    "public",
    "public",
    "public",
    "public",
    "public",
    "public",
    "public",
    "public",
    "public",
    "public",
    "public",
    "public",
    "public",
    "controlled_or_restricted",
    "public"
  ),
  stringsAsFactors = FALSE
)

source_access <- read_tsv("data/source_access_summary.tsv")
source_info <- merge(
  source_info,
  source_access[, c("dataset", "source_url", "notes")],
  by.x = "source_dataset",
  by.y = "dataset",
  all.x = TRUE
)

dataset_counts <- as.data.frame(
  table(metadata$Dataset),
  stringsAsFactors = FALSE
)
names(dataset_counts) <- c("source_dataset", "number_of_cells")
donor_counts <- aggregate(Donor_ID ~ Dataset, metadata, function(x) {
  length(unique(x))
})
names(donor_counts) <- c("source_dataset", "number_of_donors")
sample_counts <- aggregate(Sample_ID ~ Dataset, metadata, function(x) {
  length(unique(x))
})
names(sample_counts) <- c("source_dataset", "number_of_biological_samples")
library_counts <- aggregate(Library_ID ~ Dataset, metadata, function(x) {
  length(unique(x))
})
names(library_counts) <- c("source_dataset", "number_of_library_records")
region_counts <- aggregate(BrainRegion ~ Dataset, metadata, function(x) {
  length(unique(x))
})
names(region_counts) <- c("source_dataset", "number_of_brain_regions")
stage_counts <- aggregate(AgeIntervalID ~ Dataset, metadata, function(x) {
  length(unique(x))
})
names(stage_counts) <- c("source_dataset", "number_of_age_intervals")

source_summary <- Reduce(
  function(x, y) merge(x, y, by = "source_dataset", all.x = TRUE),
  list(
    source_info,
    dataset_counts,
    donor_counts,
    sample_counts,
    library_counts,
    region_counts,
    stage_counts
  )
)
source_summary$access_level <- ifelse(
  source_summary$access_level == "public",
  "Accessible",
  "Restricted"
)
source_summary <- source_summary[
  order(-source_summary$number_of_cells),
]
write_tsv(source_summary, file.path(out_dir, "source_dataset_summary.tsv"))

sample_meta <- aggregate(
  Cells ~ Dataset +
    Donor_ID +
    Sample_ID +
    Technology +
    Sequence +
    BrainRegion +
    Stage +
    DevelopmentStage +
    AgeIntervalID +
    AgeInterval +
    AgeRange +
    Age +
    Sex,
  metadata,
  length
)
names(sample_meta) <- c(
  "source_dataset",
  "donor_id",
  "sample_id",
  "sequencing_platform",
  "sequencing_modality",
  "brain_region_harmonized",
  "legacy_stage",
  "legacy_development_stage",
  "age_interval_id",
  "age_interval",
  "age_range",
  "reported_age",
  "sex_raw",
  "number_of_cells"
)
sample_meta$sex_standardized <- sample_meta$sex_raw
sample_meta$sex_standardized[
  !(sample_meta$sex_standardized %in% c("Female", "Male"))
] <- "Not reported"
sample_meta$sequencing_modality_standardized <- ifelse(
  grepl("^scRNA", sample_meta$sequencing_modality),
  "scRNA-seq",
  ifelse(
    grepl("^snRNA", sample_meta$sequencing_modality),
    "snRNA-seq",
    "Other/unknown"
  )
)
collapse_one <- function(x, mixed_value = "Multiple") {
  x <- unique(as.character(x[!is.na(x) & x != ""]))
  if (length(x) == 0) {
    return("Not reported")
  }
  if (length(x) == 1) {
    return(x)
  }
  mixed_value
}
sample_level <- do.call(
  rbind,
  lapply(split(sample_meta, sample_meta$sample_id), function(x) {
    sex <- collapse_one(x$sex_standardized, mixed_value = "Not reported")
    data.frame(
      sample_id = unique(x$sample_id)[1],
      donor_id = collapse_one(x$donor_id),
      sex_standardized = sex,
      sequencing_modality_standardized = collapse_one(
        x$sequencing_modality_standardized
      ),
      sequencing_platform = collapse_one(x$sequencing_platform),
      source_dataset = collapse_one(x$source_dataset),
      stringsAsFactors = FALSE
    )
  })
)

parse_reported_age <- function(x) {
  y <- trimws(as.character(x))
  is_pcw <- grepl("PCW", y, ignore.case = TRUE)
  nums <- regmatches(y, gregexpr("[0-9]+\\.?[0-9]*", y))
  value <- vapply(
    nums,
    function(v) {
      if (length(v) == 0) {
        return(NA_real_)
      }
      mean(as.numeric(v), na.rm = TRUE)
    },
    numeric(1)
  )
  data.frame(
    reported_age = y,
    age_value = value,
    age_unit = ifelse(is_pcw, "PCW", "years"),
    age_sort = ifelse(is_pcw, value - 1000, value),
    stringsAsFactors = FALSE
  )
}

age_lookup <- unique(parse_reported_age(sample_meta$reported_age))
reported_age_cells <- aggregate(
  number_of_cells ~ reported_age,
  sample_meta,
  sum
)
reported_age_samples <- aggregate(
  sample_id ~ reported_age,
  sample_meta,
  function(x) length(unique(x))
)
names(reported_age_samples)[2] <- "sample_count"
reported_age_summary <- merge(
  reported_age_cells,
  reported_age_samples,
  by = "reported_age"
)
reported_age_summary <- merge(
  reported_age_summary,
  age_lookup,
  by = "reported_age",
  all.x = TRUE
)
reported_age_summary <- reported_age_summary[
  order(reported_age_summary$age_sort, reported_age_summary$reported_age),
]
reported_age_summary$age_display_order <- seq_len(nrow(reported_age_summary))
write_tsv(
  reported_age_summary,
  file.path(out_dir, "figure1_reported_age_summary.tsv")
)

reported_age_sex <- aggregate(
  sample_id ~ reported_age + sex_standardized,
  sample_meta,
  function(x) length(unique(x))
)
names(reported_age_sex)[3] <- "sample_count"
reported_age_sex <- merge(
  reported_age_sex,
  reported_age_summary[, c("reported_age", "age_sort", "age_display_order")],
  by = "reported_age",
  all.x = TRUE
)
reported_age_sex <- reported_age_sex[
  order(
    reported_age_sex$age_sort,
    reported_age_sex$reported_age,
    reported_age_sex$sex_standardized
  ),
]
write_tsv(
  reported_age_sex,
  file.path(out_dir, "figure1_reported_age_sex_sample_counts.tsv")
)

celltype_counts <- NULL
cluster_count <- NA_integer_
if (file.exists(object_plot_path)) {
  object_plot <- readRDS(object_plot_path)
  md <- object_plot@meta.data
  if ("CellType" %in% names(md)) {
    celltype_counts <- as.data.frame(
      table(md$CellType),
      stringsAsFactors = FALSE
    )
    names(celltype_counts) <- c("cell_type", "cell_count")
    celltype_counts <- celltype_counts[order(-celltype_counts$cell_count), ]
  }
  if ("seurat_clusters" %in% names(md)) {
    cluster_count <- length(unique(md$seurat_clusters))
  }
}

celltype_count <- if (!is.null(celltype_counts)) {
  nrow(celltype_counts)
} else {
  NA_integer_
}
sex_unknown_cells <- sum(
  metadata$Sex[!(metadata$Sex %in% c("Female", "Male"))] %in% metadata$Sex
)
numeric_summary <- data.frame(
  metric = c(
    "source_dataset_count",
    "cell_or_nucleus_count",
    "donor_count",
    "biological_sample_count",
    "library_record_count",
    "age_interval_count",
    "brain_region_count",
    "seurat_cluster_count",
    "major_cell_type_count",
    "male_cells",
    "female_cells",
    "unknown_or_nonstandard_sex_cells",
    "male_biological_samples",
    "female_biological_samples",
    "unknown_or_nonstandard_sex_biological_samples"
  ),
  value = c(
    length(unique(metadata$Dataset)),
    nrow(metadata),
    length(unique(metadata$Donor_ID)),
    length(unique(metadata$Sample_ID)),
    length(unique(metadata$Library_ID)),
    length(unique(metadata$AgeIntervalID)),
    length(unique(metadata$BrainRegion)),
    cluster_count,
    celltype_count,
    sum(metadata$Sex == "Male"),
    sum(metadata$Sex == "Female"),
    sum(!(metadata$Sex %in% c("Female", "Male"))),
    sum(sample_level$sex_standardized == "Male"),
    sum(sample_level$sex_standardized == "Female"),
    sum(sample_level$sex_standardized == "Not reported")
  ),
  stringsAsFactors = FALSE
)
write_tsv(numeric_summary, file.path(out_dir, "numeric_summary.tsv"))

verify_crossref <- function(doi) {
  if (is.na(doi) || doi == "") {
    return(list(
      status = "needs_source_publication_confirmation",
      title = NA_character_,
      year = NA_character_,
      url = NA_character_
    ))
  }
  api <- paste0(
    "https://api.crossref.org/works/",
    utils::URLencode(doi, reserved = TRUE)
  )
  out <- tryCatch(jsonlite::fromJSON(api), error = function(e) NULL)
  if (is.null(out) || is.null(out$message$title)) {
    return(list(
      status = "doi_unresolved_or_crossref_unavailable",
      title = NA_character_,
      year = NA_character_,
      url = paste0("https://doi.org/", doi)
    ))
  }
  year <- NA_character_
  if (!is.null(out$message$issued[["date-parts"]])) {
    year <- as.character(out$message$issued[["date-parts"]][[1]][1])
  }
  list(
    status = "verified_by_crossref",
    title = paste(out$message$title, collapse = " "),
    year = year,
    url = paste0("https://doi.org/", doi)
  )
}

reference_audit_path <- file.path(out_dir, "reference_audit.tsv")
if (file.exists(reference_audit_path)) {
  reference_audit <- read_tsv(reference_audit_path)
} else {
  audit_rows <- lapply(seq_len(nrow(source_info)), function(i) {
    doi <- source_info$publication_doi[[i]]
    v <- verify_crossref(doi)
    data.frame(
      source_dataset = source_info$source_dataset[[i]],
      source_accession = source_info$source_accession[[i]],
      source_repository = source_info$source_repository[[i]],
      publication_doi = ifelse(is.na(doi), "", doi),
      verification_status = v$status,
      verified_title = v$title %||% NA_character_,
      publication_year = v$year %||% NA_character_,
      source_url = source_info$source_url[[i]] %||% NA_character_,
      reference_or_access_note = source_info$notes[[i]] %||% NA_character_,
      stringsAsFactors = FALSE
    )
  })
  reference_audit <- do.call(rbind, audit_rows)
  write_tsv(reference_audit, reference_audit_path)
}

workflow_steps <- data.frame(
  step = factor(1:6),
  label = c(
    "Dataset\ncollection",
    "Metadata\ncuration",
    "Age/region\nmapping",
    "Expression\nprocessing",
    "RPCA\nintegration",
    "Cell-type\nannotation"
  ),
  detail = c(
    "27 datasets",
    "520 donors;\n603 samples",
    sprintf("15 intervals;\n%s regions", length(unique(metadata$BrainRegion))),
    "filtering;\nPCA",
    "2.24M cells;\nLISI",
    "120 clusters;\n9 cell types"
  ),
  x = seq(1, 11, by = 2),
  y = rep(1, 6),
  stringsAsFactors = FALSE
)
write_tsv(workflow_steps, file.path(out_dir, "figure1_workflow_steps.tsv"))

arrow_steps <- data.frame(
  x = workflow_steps$x[-nrow(workflow_steps)],
  y = workflow_steps$y[-nrow(workflow_steps)],
  xend = workflow_steps$x[-1],
  yend = workflow_steps$y[-1]
)
box_w <- 1.50
box_h <- 0.70

fig1c <- ggplot(workflow_steps, aes(x, y)) +
  geom_segment(
    data = arrow_steps,
    aes(x = x + box_w / 2, xend = xend - box_w / 2, y = y, yend = yend),
    linewidth = 0.45,
    arrow = arrow(length = unit(2.4, "mm"), type = "closed"),
    colour = "#56616D"
  ) +
  geom_rect(
    aes(
      xmin = x - box_w / 2,
      xmax = x + box_w / 2,
      ymin = y - box_h / 2,
      ymax = y + box_h / 2,
      fill = as.integer(step)
    ),
    colour = "#334155",
    linewidth = 0.35
  ) +
  geom_text(
    aes(label = label),
    nudge_y = 0.13,
    size = 2.05,
    fontface = "bold",
    lineheight = 0.86,
    colour = "#111827"
  ) +
  geom_text(
    aes(label = detail),
    nudge_y = -0.15,
    size = 1.8,
    lineheight = 0.86,
    colour = "#374151"
  ) +
  scale_fill_gradientn(
    colours = c("#E8F3F8", "#DDF0E7", "#F8E8D1"),
    guide = "none"
  ) +
  coord_cartesian(xlim = c(0.0, 12.0), ylim = c(0.5, 1.5), clip = "off") +
  theme_void(base_family = "Arial") +
  # theme(plot.margin = margin(6, 8, 4, 8)) +
  theme(plot.margin = margin(0, 0, 0, 0))

sex_totals <- as.data.frame(
  table(sample_level$sex_standardized),
  stringsAsFactors = FALSE
)
names(sex_totals) <- c("group", "count")
sex_totals$group <- as.character(sex_totals$group)
write_tsv(sex_totals, file.path(out_dir, "figure1_sex_sample_totals.tsv"))

modality_totals <- as.data.frame(
  table(sample_level$sequencing_modality_standardized),
  stringsAsFactors = FALSE
)
names(modality_totals) <- c("group", "count")
write_tsv(
  modality_totals,
  file.path(out_dir, "figure1_modality_sample_totals.tsv")
)

technology_totals <- as.data.frame(
  table(sample_level$sequencing_platform),
  stringsAsFactors = FALSE
)
names(technology_totals) <- c("group", "count")
technology_totals <- technology_totals[order(-technology_totals$count), ]
write_tsv(
  technology_totals,
  file.path(out_dir, "figure1_technology_sample_totals.tsv")
)

source_summary$access_group <- source_summary$access_level
access_totals <- as.data.frame(
  table(source_summary$access_group),
  stringsAsFactors = FALSE
)
names(access_totals) <- c("group", "count")
write_tsv(
  access_totals,
  file.path(out_dir, "figure1_access_dataset_totals.tsv")
)

make_pie <- function(df, title, palette) {
  df <- df[order(df$group), ]
  df$fraction <- df$count / sum(df$count)
  df$legend_label <- paste0(
    df$group,
    " (",
    df$count,
    ", ",
    sprintf("%.1f", 100 * df$fraction),
    "%)"
  )
  df$legend_label <- factor(df$legend_label, levels = df$legend_label)
  ggplot(df, aes(x = 1, y = count, fill = group)) +
    geom_col(width = 1, colour = "white", linewidth = 0.25) +
    coord_polar(theta = "y") +
    scale_fill_manual(
      values = palette,
      breaks = df$group,
      labels = as.character(df$legend_label),
      name = NULL
    ) +
    labs(title = title) +
    theme_void(base_family = "Arial") +
    theme(
      plot.title = element_text(size = 7, face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.text = element_text(size = 5.2),
      legend.key.size = unit(2.6, "mm"),
      legend.spacing.x = unit(0.6, "mm"),
      legend.margin = margin(0, 0, 0, -4),
      plot.margin = margin(1, 0, 1, 0)
    )
}

sex_palette <- c(
  "Female" = "#B24AA6",
  "Male" = "#28327F",
  "Not reported" = "#8A8A8A"
)
p_sex <- make_pie(sex_totals, "C  Reported sex", sex_palette)
p_modality <- make_pie(
  modality_totals,
  "Sequencing modality",
  c(
    "scRNA-seq" = "#6BAED6",
    "snRNA-seq" = "#2171B5",
    "Other/unknown" = "#BDBDBD"
  )
)
p_technology <- make_pie(
  technology_totals,
  "Sequencing technology",
  c(
    "STRT-seq" = "#66C2A5",
    "10X Genomics" = "#FC8D62",
    "Fluidigm C1" = "#8DA0CB",
    "SNdrop-seq" = "#E78AC3",
    "EasySci-RNA" = "#A6D854"
  )
)
p_access <- make_pie(
  access_totals,
  "Dataset access",
  c(
    "Accessible" = "#2CA25F",
    "Restricted" = "#756BB1"
  )
)
# pies_fig1 <- ((p_sex + labs(title = NULL)) |
#   (p_modality + labs(title = NULL))) /
#   ((p_technology + labs(title = NULL)) |
#     (p_access + labs(title = NULL)))
pies_fig1 <- (p_access + labs(title = NULL)) |
  (p_technology + labs(title = NULL)) |
  (p_modality + labs(title = NULL)) |
  (p_sex + labs(title = NULL))
age_bar <- reported_age_summary
age_bar$reported_age <- factor(
  age_bar$reported_age,
  levels = reported_age_summary$reported_age
)
reported_age_sex$reported_age <- factor(
  reported_age_sex$reported_age,
  levels = reported_age_summary$reported_age
)
reported_age_sex$sex_standardized <- factor(
  reported_age_sex$sex_standardized,
  levels = c("Female", "Male", "Not reported")
)
tick_idx <- seq(1, nrow(reported_age_summary), by = 6)
tick_idx <- sort(unique(c(tick_idx, nrow(reported_age_summary))))
age_tick_labels <- rep("", nrow(reported_age_summary))
age_tick_labels[tick_idx] <- as.character(reported_age_summary$reported_age[
  tick_idx
])
names(age_tick_labels) <- as.character(reported_age_summary$reported_age)

p_age_cells <- ggplot(age_bar, aes(reported_age, number_of_cells)) +
  geom_col(
    # fill = "#5AA0D8",
    fill = "black",
    width = 0.72
  ) +
  scale_y_continuous(labels = fmt_int, expand = expansion(mult = c(0, 0.06))) +
  scale_x_discrete(labels = age_tick_labels) +
  labs(x = NULL, y = "Cells/nuclei") +
  theme_paper(7) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(
      angle = 0,
      vjust = 0.5,
      margin = margin(r = 5)
    ),
    plot.margin = margin(0, 1, 0, 10)
  )

p_age_samples <- ggplot(
  reported_age_sex,
  aes(
    reported_age,
    sample_count,
    colour = sex_standardized,
    size = sample_count
  )
) +
  geom_point(
    position = position_dodge(width = 0.7),
    alpha = 0.9
  ) +
  scale_colour_manual(values = sex_palette, name = "Samples by sex") +
  scale_size_continuous(range = c(0.7, 3.0), name = "Biological samples") +
  scale_y_reverse(labels = fmt_int, expand = expansion(mult = c(0.12, 0))) +
  scale_x_discrete(labels = age_tick_labels) +
  labs(
    x = "Reported age (PCW or years)",
    y = "Samples"
  ) +
  theme_paper(7) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 5.8),
    legend.text = element_text(size = 5.4),
    legend.key.size = unit(2.8, "mm"),
    legend.margin = margin(-3, 0, 0, 0),
    legend.box.margin = margin(-6, 0, -6, 0),
    legend.box.spacing = unit(0, "mm"),
    legend.spacing.y = unit(0, "mm"),
    legend.spacing.x = unit(0.8, "mm"),
    axis.title.x = element_text(margin = margin(t = 1)),
    axis.title.y = element_text(
      angle = 0,
      vjust = 0.5,
      margin = margin(r = 5)
    ),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 4.7),
    plot.margin = margin(0, 1, -8, 10)
  )

fig1b <- p_age_cells /
  p_age_samples +
  plot_layout(heights = c(0.46, 0.54), guides = "collect") &
  theme(legend.position = "bottom")

fig1a <- pies_fig1 & theme(plot.title = element_blank())
fig1_combined <- (wrap_elements(full = fig1a) /
  wrap_elements(full = fig1b) /
  wrap_elements(full = fig1c)) +
  plot_layout(heights = c(0.25, 0.50, 0.16)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.margin = margin(1, 2, 1, 2),
    plot.tag = element_text(size = 8, family = "Arial"),
    plot.tag.position = c(0.01, 0.99)
  )

grDevices::cairo_pdf(
  file.path(fig_dir, "fig1.pdf"),
  width = 6.5,
  height = 4,
  family = "Arial"
)
print(fig1_combined)
dev.off()

grDevices::png(
  file.path(fig_dir, "fig1.png"),
  width = 6.5,
  height = 4,
  family = "Arial",
  units = "in",
  res = 600
)
print(fig1_combined)
dev.off()
