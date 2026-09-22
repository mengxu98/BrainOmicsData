#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(data.table)
  library(ggplot2)
  library(scop)
  library(Seurat)
  library(SeuratObject)
})

source("functions/processed_object.R")
source("functions/integration.R")

plot_lineage <- function(lineage_id, lineage_dir) {
  if (!dir.exists(lineage_dir) ||
    !lineage_id %chin% c("excitatory", "inhibitory")) {
    stop("A completed excitatory or inhibitory lineage directory is required")
  }

  paths <- list(
    object = file.path(lineage_dir, "objects_lineage_integrated.rds"),
    assignments = file.path(
      lineage_dir, "clustering", "lineage_cluster_assignments.tsv.gz"
    ),
    subtype_clusters = file.path(
      lineage_dir, "subtype_annotation",
      "lineage_cluster_subtype_annotation.tsv"
    ),
    subtype_cells = file.path(
      lineage_dir, "subtype_annotation",
      "lineage_cell_subtype_annotation.tsv.gz"
    ),
    subtype_contract = file.path(
      lineage_dir, "subtype_annotation",
      "lineage_subtype_annotation_contract.tsv"
    ),
    marker_contract = file.path(
      lineage_dir, "markers", "lineage_marker_contract.tsv"
    ),
    stability = file.path(
      lineage_dir, "clustering", "lineage_cluster_stability.tsv"
    )
  )
  missing <- unlist(paths)[!file.exists(unlist(paths))]
  if (length(missing) > 0L ||
    !file.exists(file.path(lineage_dir, "_VALIDATED"))) {
    stop("Validated lineage plotting inputs are missing: ", paste(missing, collapse = ", "))
  }

  final_figure_dir <- file.path(lineage_dir, "figures")
  plot_object_file <- file.path(lineage_dir, "objects_lineage_plot.rds")
  plot_manifest_file <- file.path(lineage_dir, "objects_lineage_plot_manifest.tsv")
  figure_dir <- paste0(final_figure_dir, ".tmp.", Sys.getpid())
  unlink(figure_dir, recursive = TRUE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  assignments <- fread(paths$assignments, sep = "\t")
  subtype_clusters <- fread(paths$subtype_clusters, sep = "\t")
  subtype_cells <- fread(paths$subtype_cells, sep = "\t")
  required_assignment <- c(
    "Cells", "Global_Cluster", "Main_CellType", "Lineage_Cluster"
  )
  required_subtype <- c(
    "Cells", "Lineage_Cluster", "Lineage_Subtype", "Assignment_Status"
  )
  if (any(!required_assignment %in% names(assignments)) ||
    any(!required_subtype %in% names(subtype_cells)) ||
    anyDuplicated(assignments$Cells) || anyDuplicated(subtype_cells$Cells) ||
    !identical(assignments$Cells, subtype_cells$Cells) ||
    !identical(assignments$Lineage_Cluster, subtype_cells$Lineage_Cluster)) {
    stop("Complete lineage and subtype sidecars differ in cell order or cluster")
  }

  marker_groups <- if (lineage_id == "inhibitory") {
    list(
      `PVALB interneurons` = c("PVALB", "KCNC1", "GAD1", "GAD2"),
      `SST interneurons` = c("SST", "GRM1", "GAD1", "GAD2"),
      `VIP interneurons` = c("VIP", "ADARB2", "GAD1", "GAD2"),
      `SNCG interneurons` = c("SNCG", "CCK", "RELN", "PAX6"),
      `LAMP5 interneurons` = c("LAMP5", "RELN", "CXCL14", "NDNF"),
      `MEIS2/PBX3 LGE/striatal-like inhibitory neurons` = c(
        "MEIS2", "FOXP1", "PBX3", "BCL11B"
      ),
      `Striatal inhibitory neurons` = c("PPP1R1B", "DRD1", "DRD2", "FOXP1"),
      `Immature inhibitory neurons` = c("DLX6-AS1", "SOX4", "SOX11", "DCX"),
      `Purkinje neurons` = c("PCP2", "CALB1", "GRID2", "ITPR1"),
      `Cerebellar inhibitory neurons` = c("PAX2", "TFAP2A", "TFAP2B", "SLC6A5")
    )
  } else {
    list(
      `Upper-layer IT-like neurons` = c("CUX1", "CUX2", "SATB2", "POU3F2"),
      `RORB-positive L4/5 IT-like neurons` = c("RORB", "FOXP1", "SATB2", "SCN4B"),
      `Layer 5 IT-like neurons` = c("BCL11B", "FOXP1", "SATB2", "ETV1"),
      `L5 ET-like neurons` = c("BCL11B", "FEZF2", "POU3F1", "ETV1"),
      `Layer 6 IT-like neurons` = c("FOXP1", "TLE4", "SOX5", "SATB2"),
      `Layer 6 IT CAR3-like neurons` = c("CAR3", "SULF2", "FOXP1", "TLE4"),
      `L6 CT-like neurons` = c("TLE4", "FOXP2", "TBR1", "SOX5"),
      `L5/6 NP-like neurons` = c("TSHZ2", "FOXP2", "SOX5", "TLE4"),
      `Layer 6b neurons` = c("CTGF", "CPLX3", "RELN", "TLE4"),
      `Hippocampal pyramidal neurons` = c("WFS1", "GRIK4", "CPNE4", "PCP4"),
      `Dentate granule neurons` = c("PROX1", "CALB1", "PDLIM5", "NEUROD6"),
      `Cerebellar granule neurons` = c("GABRA6", "GRM4", "ZIC1", "ZIC2"),
      `Thalamic excitatory neurons` = c("SLC17A6", "TCF7L2", "FOXP2", "GBX2"),
      `Immature excitatory state` = c("NEUROD1", "NEUROD2", "SOX11", "DCX")
    )
  }
  exclusion_groups <- list(
    `Opposing neuronal lineage` = if (lineage_id == "inhibitory") {
      c("SLC17A7", "SLC17A6", "CAMK2A")
    } else {
      c("GAD1", "GAD2", "SLC32A1")
    },
    Astroglial = c("AQP4", "SLC1A2", "SLC1A3"),
    Oligodendroglial = c("PLP1", "MOBP", "MOG"),
    Immune = c("P2RY12", "C1QA", "PTPRC")
  )
  marker_spec <- rbindlist(lapply(names(marker_groups), function(group) {
    data.table(Marker_Group = group, Gene = marker_groups[[group]])
  }))
  exclusion_spec <- rbindlist(lapply(names(exclusion_groups), function(group) {
    data.table(Marker_Group = group, Gene = exclusion_groups[[group]])
  }))
  requested_features <- unique(c(marker_spec$Gene, exclusion_spec$Gene))

  reuse_plot_object <- FALSE
  if (file.exists(plot_object_file) && file.exists(plot_manifest_file)) {
    existing_manifest <- fread(plot_manifest_file, sep = "\t")
    required_manifest_columns <- c("Requested_Features", "Available_Features")
    manifest_ok <- nrow(existing_manifest) == 1L &&
      all(required_manifest_columns %in% names(existing_manifest))
    if (manifest_ok) {
      objects_plot <- readRDS(plot_object_file)
      manifest_requested_features <- strsplit(
        existing_manifest$Requested_Features[[1L]], ";",
        fixed = TRUE
      )[[1L]]
      manifest_available_features <- strsplit(
        existing_manifest$Available_Features[[1L]], ";",
        fixed = TRUE
      )[[1L]]
      reuse_checks <- c(
        seurat = inherits(objects_plot, "Seurat"),
        cells = ncol(objects_plot) == nrow(assignments),
        cell_order = identical(colnames(objects_plot), assignments$Cells),
        requested_features = identical(
          requested_features, manifest_requested_features
        ),
        available_features = identical(
          rownames(objects_plot), manifest_available_features
        ),
        reduction = "umap.lineage.rpca" %in% Reductions(objects_plot)
      )
      reuse_plot_object <- all(reuse_checks)
      if (!reuse_plot_object) {
        rm(objects_plot)
        gc()
      }
    }
  }

  if (!reuse_plot_object) {
    objects <- readRDS(paths$object)
    if (!inherits(objects, "Seurat") ||
      !identical(colnames(objects), assignments$Cells) ||
      !"umap.lineage.rpca" %in% Reductions(objects) ||
      !identical(
        rownames(Embeddings(objects, "umap.lineage.rpca")),
        assignments$Cells
      )) {
      stop("Integrated lineage object differs from validated assignments")
    }
    available_features <- requested_features[requested_features %in% rownames(objects)]
    if (length(available_features) == 0L) {
      stop("No lineage marker feature is present in the integrated object")
    }
    objects_plot <- objects[available_features, ]
    rm(objects)
    gc()
    for (reduction in setdiff(Reductions(objects_plot), "umap.lineage.rpca")) {
      objects_plot[[reduction]] <- NULL
    }
    for (assay in setdiff(Assays(objects_plot), "RNA")) {
      objects_plot[[assay]] <- NULL
    }
    DefaultAssay(objects_plot) <- "RNA"
    objects_plot@graphs <- list()
    objects_plot@neighbors <- list()
    objects_plot@commands <- list()
    objects_plot@tools <- list()
    data_layers <- Layers(objects_plot, assay = "RNA", search = "^data")
    if (length(data_layers) > 1L) {
      objects_plot <- JoinLayers(
        objects_plot,
        assay = "RNA", layers = "data", new = "data"
      )
    }
    available_features <- rownames(objects_plot)
    manifest <- data.table(
      File = basename(plot_object_file),
      Cells = ncol(objects_plot),
      Features = nrow(objects_plot),
      Requested_Features = paste(requested_features, collapse = ";"),
      Available_Features = paste(available_features, collapse = ";")
    )
  }

  objects_plot$Lineage_Cluster <- factor(
    assignments$Lineage_Cluster,
    levels = sort(unique(assignments$Lineage_Cluster))
  )
  display_subtype <- ifelse(
    subtype_cells$Assignment_Status == "supported",
    subtype_cells$Lineage_Subtype,
    ifelse(
      subtype_cells$Assignment_Status == "excluded",
      "Excluded from subclass/state naming",
      "Unresolved"
    )
  )
  subtype_levels <- c(
    names(marker_groups)[names(marker_groups) %chin% display_subtype],
    intersect(
      c("Unresolved", "Excluded from subclass/state naming"),
      display_subtype
    )
  )
  objects_plot$Lineage_Subtype <- factor(display_subtype, levels = subtype_levels)
  Idents(objects_plot) <- "Lineage_Subtype"

  if (!reuse_plot_object) {
    temporary <- paste0(plot_object_file, ".tmp.", Sys.getpid())
    saveRDS(objects_plot, temporary, compress = FALSE)
    if (!file.rename(temporary, plot_object_file)) {
      unlink(temporary)
      stop("Could not publish complete-cell lineage plotting object")
    }
    fwrite(manifest, plot_manifest_file, sep = "\t", quote = FALSE)
  }

  available_features <- rownames(objects_plot)
  marker_spec <- marker_spec[Gene %chin% available_features]
  exclusion_spec <- exclusion_spec[Gene %chin% available_features]
  if (nrow(marker_spec) == 0L || nrow(exclusion_spec) == 0L) {
    stop("Lineage marker or exclusion-marker panel has no available features")
  }

  cluster_levels <- levels(objects_plot$Lineage_Cluster)
  cluster_colors <- setNames(
    grDevices::hcl.colors(length(cluster_levels), "Dynamic"), cluster_levels
  )
  subtype_colors <- setNames(
    grDevices::hcl.colors(
      sum(!subtype_levels %chin% c(
        "Unresolved", "Excluded from subclass/state naming"
      )),
      "Dark 3"
    ),
    subtype_levels[!subtype_levels %chin% c(
      "Unresolved", "Excluded from subclass/state naming"
    )]
  )
  subtype_colors <- c(
    subtype_colors,
    `Unresolved` = "#BDBDBD",
    `Excluded from subclass/state naming` = "#4D4D4D"
  )
  subtype_colors <- subtype_colors[subtype_levels]

  plot_warnings <- data.table(Stage = character(), Warning = character())
  capture_warnings <- function(expression, stage) {
    withCallingHandlers(
      expression,
      warning = function(condition) {
        plot_warnings <<- rbind(
          plot_warnings,
          data.table(Stage = stage, Warning = conditionMessage(condition))
        )
        invokeRestart("muffleWarning")
      }
    )
  }

  save_scop_plot <- function(plot, stem, width, height) {
    ggplot2::ggsave(
      paste0(stem, ".pdf"), plot,
      device = grDevices::cairo_pdf, width = width, height = height,
      units = "in", family = "Arial", bg = "white"
    )
  }

  save_scop_heatmap <- function(plot, stem, width, height) {
    render_heatmap <- function(plot) {
      if (inherits(plot, c("Heatmap", "HeatmapList"))) {
        ComplexHeatmap::draw(plot)
      } else {
        print(plot)
      }
    }
    grDevices::cairo_pdf(
      paste0(stem, ".pdf"),
      width = width, height = height,
      family = "Arial", onefile = FALSE
    )
    render_heatmap(plot)
    grDevices::dev.off()
  }

  # scop::GroupHeatmap provides the requested expression and percent-dot panel,
  # but its group legend duplicates every displayed lineage-cluster column. Keep
  # the color strip and column names while removing only that redundant legend.
  group_heatmap_without_group_legend <- function(...) {
    group_heatmap <- scop::GroupHeatmap
    body_parts <- as.list(body(group_heatmap))
    group_legend_expression <- which(vapply(
      body_parts,
      function(expression) {
        grepl(
          "seq_along(raw_group_by)",
          paste(deparse(expression, width.cutoff = 500L), collapse = " "),
          fixed = TRUE
        )
      },
      logical(1L)
    ))
    if (length(group_legend_expression) != 1L) {
      stop("The installed scop::GroupHeatmap group-legend contract changed")
    }
    body_parts <- append(
      body_parts, list(quote(lgd[raw_group_by] <- NULL)),
      after = group_legend_expression
    )
    body(group_heatmap) <- as.call(body_parts)
    environment(group_heatmap) <- environment(scop::GroupHeatmap)
    group_heatmap(...)
  }

  cluster_umap <- capture_warnings(scop::CellDimPlot(
    objects_plot,
    reduction = "umap.lineage.rpca",
    group.by = "Lineage_Cluster",
    palcolor = cluster_colors,
    label = TRUE,
    raster = TRUE,
    xlab = "UMAP_1",
    ylab = "UMAP_2",
    theme_use = "theme_blank"
  ), "cluster_umap")
  capture_warnings(save_scop_plot(
    cluster_umap,
    file.path(figure_dir, paste0(lineage_id, "_lineage_clusters")),
    width = 8.0, height = 6.0
  ), "cluster_umap_render")

  subtype_umap <- capture_warnings(scop::CellDimPlot(
    objects_plot,
    reduction = "umap.lineage.rpca",
    group.by = "Lineage_Subtype",
    palcolor = subtype_colors,
    label = FALSE,
    raster = TRUE,
    xlab = "UMAP_1",
    ylab = "UMAP_2",
    theme_use = "theme_blank"
  ), "subtype_umap")
  capture_warnings(save_scop_plot(
    subtype_umap,
    file.path(figure_dir, paste0(lineage_id, "_lineage_subtypes")),
    width = 8.0, height = 6.0
  ), "subtype_umap_render")

  draw_group_heatmap <- function(spec, group_by, group_colors, column_title,
                                 stem, width, height,
                                 suppress_group_legend = FALSE) {
    feature_split <- factor(spec$Marker_Group, levels = unique(spec$Marker_Group))
    split_colors <- setNames(
      grDevices::hcl.colors(nlevels(feature_split), "Dark 3"),
      levels(feature_split)
    )
    heatmap_function <- if (suppress_group_legend) {
      group_heatmap_without_group_legend
    } else {
      scop::GroupHeatmap
    }
    result <- capture_warnings(heatmap_function(
      srt = objects_plot,
      features = spec$Gene,
      layer = "data",
      group.by = group_by,
      group_palcolor = group_colors,
      cell_annotation_palcolor = group_colors,
      lib_normalize = FALSE,
      exp_method = "zscore",
      exp_legend_title = "Z-score",
      feature_split = feature_split,
      feature_split_palcolor = split_colors,
      heatmap_palette = "Spectral",
      heatmap_border = TRUE,
      cell_annotation_border = FALSE,
      feature_annotation_border = FALSE,
      add_dot = TRUE,
      dot_size = grid::unit(6, "mm"),
      nlabel = 0,
      show_row_names = TRUE,
      show_column_names = TRUE,
      column_title = column_title,
      column_names_side = "bottom",
      column_names_rot = 90,
      border = TRUE,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      use_raster = FALSE,
      legend.position = "right",
      ht_params = list(
        row_names_gp = grid::gpar(
          fontface = "italic", fontsize = 7, fontfamily = "Arial"
        )
      ),
      verbose = TRUE
    ), paste0("group_heatmap_", basename(stem)))
    capture_warnings(
      save_scop_heatmap(result$plot, stem, width, height),
      paste0("group_heatmap_render_", basename(stem))
    )
  }

  draw_group_heatmap(
    marker_spec,
    "Lineage_Cluster",
    cluster_colors,
    "Lineage cluster",
    file.path(figure_dir, paste0(lineage_id, "_positive_markers_by_cluster")),
    width = max(12, length(cluster_levels) * 0.22 + 6),
    height = max(9, nrow(marker_spec) * 0.14 + 4),
    suppress_group_legend = TRUE
  )
  draw_group_heatmap(
    marker_spec,
    "Lineage_Subtype",
    subtype_colors,
    "Supported broad lineage subclass/state",
    file.path(figure_dir, paste0(lineage_id, "_positive_markers_by_subtype")),
    width = max(9, length(subtype_levels) * 0.55 + 5),
    height = max(9, nrow(marker_spec) * 0.14 + 4)
  )
  draw_group_heatmap(
    exclusion_spec,
    "Lineage_Cluster",
    cluster_colors,
    "Lineage cluster",
    file.path(figure_dir, paste0(lineage_id, "_exclusion_markers_by_cluster")),
    width = max(12, length(cluster_levels) * 0.22 + 6),
    height = max(7, nrow(exclusion_spec) * 0.17 + 4),
    suppress_group_legend = TRUE
  )

  stability <- fread(paths$stability, sep = "\t")
  if (!all(c(
    "Assignment_A", "Assignment_B", "Adjusted_Rand_Index"
  ) %in% names(stability))) {
    stop("Lineage stability table schema changed")
  }
  stability[, Comparison := paste(Assignment_A, Assignment_B, sep = " vs ")]
  stability_plot <- ggplot(
    stability, aes(x = Comparison, y = Adjusted_Rand_Index)
  ) +
    geom_hline(yintercept = 0.8, linetype = 2, color = "grey55") +
    geom_col(width = 0.65, fill = "#3C5488") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = NULL, y = "Adjusted Rand index") +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  capture_warnings(save_scop_plot(
    stability_plot,
    file.path(figure_dir, paste0(lineage_id, "_cluster_stability")),
    width = 5.0, height = 3.8
  ), "stability_render")

  fwrite(
    unique(plot_warnings),
    file.path(figure_dir, "lineage_figure_warnings.tsv"),
    sep = "\t", quote = FALSE
  )

  figure_files <- sort(list.files(figure_dir, full.names = TRUE))
  figure_manifest <- data.table(
    File = basename(figure_files),
    Size_Bytes = as.character(file.info(figure_files)$size),
    Lineage_ID = lineage_id,
    Cells = nrow(assignments),
    Clusters = uniqueN(assignments$Lineage_Cluster),
    Supported_Subtypes = sum(
      subtype_clusters$Assignment_Status == "supported"
    ),
    Plotting_Function = ifelse(
      grepl("stability", basename(figure_files)),
      "ggplot2::ggplot",
      ifelse(
        grepl("markers", basename(figure_files)),
        "scop::GroupHeatmap",
        "scop::CellDimPlot"
      )
    )
  )
  fwrite(
    figure_manifest,
    file.path(figure_dir, "lineage_figure_manifest.tsv"),
    sep = "\t", quote = FALSE
  )
  previous_figure_dir <- paste0(final_figure_dir, ".previous.", Sys.getpid())
  unlink(previous_figure_dir, recursive = TRUE)
  if (dir.exists(final_figure_dir) &&
    !file.rename(final_figure_dir, previous_figure_dir)) {
    stop("Could not preserve the previous lineage figure directory")
  }
  if (!file.rename(figure_dir, final_figure_dir)) {
    if (dir.exists(previous_figure_dir)) {
      file.rename(previous_figure_dir, final_figure_dir)
    }
    stop("Could not atomically publish the complete lineage figure directory")
  }
  unlink(previous_figure_dir, recursive = TRUE)

  message(
    "[lineage-figures] completed ", lineage_id, " panels for ",
    nrow(assignments), " complete-lineage cells; no global annotation changed"
  )
}

lineage_root <- paste0(
  "../../data/BrainOmicsData/integration_25/",
  "lineage_analysis_20260825_v2"
)
lineage_directories <- c(
  excitatory = "excitatory_complete_primary",
  inhibitory = "inhibitory_complete_primary"
)
for (lineage_id in names(lineage_directories)) {
  plot_lineage(
    lineage_id,
    file.path(lineage_root, lineage_directories[[lineage_id]])
  )
  gc()
}
