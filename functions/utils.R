# Shared R helpers used across plotting and figure assembly.

read_tsv <- function(path) {
  utils::read.delim(
    path,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA")
  )
}

write_tsv <- function(x, path) {
  connection <- if (grepl("[.]gz$", path)) {
    gzfile(path, "wt")
  } else {
    file(path, "wt")
  }
  on.exit(close(connection))
  utils::write.table(
    x,
    connection,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
}

# Keep point order fixed across annotation plots.
reference_knn_projection <- function(
  query,
  reference,
  coordinates,
  labels = NULL,
  index = NULL,
  k = 30L
) {
  stopifnot(
    ncol(query) == ncol(reference),
    nrow(coordinates) == nrow(reference),
    identical(rownames(coordinates), rownames(reference))
  )
  neighbors <- Seurat::FindNeighbors(
    object = reference,
    query = query,
    k.param = k,
    nn.method = "annoy",
    annoy.metric = "cosine",
    n.trees = 50L,
    return.neighbor = TRUE,
    cache.index = TRUE,
    index = index,
    verbose = FALSE
  )
  idx <- neighbors@nn.idx
  projected <- vapply(
    seq_len(ncol(coordinates)),
    function(j) {
      rowMeans(matrix(coordinates[idx, j], nrow = nrow(idx)))
    },
    numeric(nrow(idx))
  )
  dimnames(projected) <- list(rownames(query), colnames(coordinates))
  result <- list(Projection = projected, Index = neighbors@alg.idx)
  if (!is.null(labels)) {
    stopifnot(length(labels) == nrow(reference), !anyNA(labels))
    levels <- unique(as.character(labels))
    votes <- vapply(
      levels,
      function(label) {
        rowMeans(matrix(labels[idx] == label, nrow = nrow(idx)))
      },
      numeric(nrow(idx))
    )
    best <- max.col(votes, ties.method = "first")
    support <- votes[cbind(seq_len(nrow(votes)), best)]
    votes[cbind(seq_len(nrow(votes)), best)] <- -Inf
    result$Prediction <- levels[best]
    result$Support <- support
    result$Margin <- support - apply(votes, 1L, max)
  }
  result
}

order_dim_plot_cells <- function(plot, cells) {
  stopifnot(
    is.data.frame(plot$data),
    !anyDuplicated(cells),
    !anyDuplicated(rownames(plot$data)),
    setequal(rownames(plot$data), cells)
  )
  plot$data <- plot$data[match(cells, rownames(plot$data)), , drop = FALSE]
  for (i in seq_along(plot$layers)) {
    layer_data <- plot$layers[[i]]$data
    if (
      is.data.frame(layer_data) &&
        nrow(layer_data) > 0L &&
        "group.by" %in% names(layer_data)
    ) {
      stopifnot(
        !anyDuplicated(rownames(layer_data)),
        all(rownames(layer_data) %in% cells)
      )
      plot$layers[[i]]$data <- layer_data[
        order(match(rownames(layer_data), cells)),
        ,
        drop = FALSE
      ]
    }
  }
  stopifnot(identical(rownames(plot$data), cells))
  plot
}

celltype_summary_group_heatmap <- function(summary, marker_spec, type_levels, flip = FALSE) {
  genes <- marker_spec$Gene
  marker_levels <- unique(as.character(marker_spec$Marker_Group))
  stopifnot(!anyDuplicated(genes), !anyDuplicated(summary[c("Group", "Gene")]))
  index <- match(
    paste(
      rep(type_levels, each = length(genes)),
      rep(genes, length(type_levels))
    ),
    paste(summary$Group, summary$Gene)
  )
  stopifnot(!anyNA(index))
  means <- matrix(
    summary$MeanLog[index],
    nrow = length(genes),
    dimnames = list(genes, type_levels)
  )
  fraction <- matrix(
    summary$PctPositive[index],
    nrow = length(genes),
    dimnames = dimnames(means)
  )
  stopifnot(
    all(is.finite(means)),
    all(is.finite(fraction)),
    all(fraction >= 0 & fraction <= 1)
  )
  carrier <- SeuratObject::CreateSeuratObject(
    counts = Matrix::Matrix(
      0,
      nrow = nrow(means),
      ncol = ncol(means),
      sparse = TRUE,
      dimnames = dimnames(means)
    ),
    assay = "Summary",
    meta.data = data.frame(
      CellType = factor(type_levels, levels = type_levels),
      Unit = "cell_type_summary",
      row.names = type_levels
    ),
    min.cells = 0,
    min.features = 0
  )
  SeuratObject::LayerData(carrier, assay = "Summary", layer = "data") <- means

  percent_legend <- ComplexHeatmap::Legend(
    title = "Percent",
    nrow = if (flip) 1 else NULL,
    title_position = if (flip) "leftcenter" else "topleft",
    labels = paste0(seq(20, 100, 20), "%"),
    type = "points",
    pch = 21,
    size = grid::unit(3 * seq(0.2, 1, 0.2), "mm"),
    legend_gp = grid::gpar(fill = "grey40", col = "black", lwd = 0.6),
    background = "transparent",
    border = FALSE,
    title_gp = grid::gpar(fontfamily = "Arial", fontsize = 9),
    labels_gp = grid::gpar(fontfamily = "Arial", fontsize = 8)
  )

  group_heatmap <- scop::GroupHeatmap
  parts <- as.list(body(group_heatmap))
  layout_at <- which(vapply(
    parts,
    function(x) {
      is.call(x) &&
        identical(x[[1]], as.name("<-")) &&
        identical(x[[2]], as.name("rendersize"))
    },
    logical(1)
  ))
  stopifnot(length(layout_at) == 1L)
  # Reserve capture space for the full cell-type row labels in the flipped panel.
  if (flip) parts[[layout_at]] <- substitute({
    original
    rendersize[["width_sum"]] <- rendersize[["width_sum"]] + 2
  }, list(original = parts[[layout_at]]))
  body(group_heatmap) <- as.call(append(
    parts,
    list(if (flip) quote(lgd <- list(ComplexHeatmap::Legend(
      title = "Z-score", col_fun = colors, direction = "horizontal",
      legend_width = grid::unit(22, "mm"), title_position = "leftcenter",
      at = c(-1, 0, 1), title_gp = grid::gpar(fontfamily = "Arial", fontsize = 8),
      labels_gp = grid::gpar(fontfamily = "Arial", fontsize = 7)), percent_legend)) else
      quote(lgd <- list(lgd[["ht"]], percent_legend, lgd[[group.by[1L]]]))),
    after = layout_at - 1L
  ))
  environment(group_heatmap) <- list2env(
    list(percent_legend = percent_legend),
    parent = environment(scop::GroupHeatmap)
  )
  p <- group_heatmap(
    carrier,
    features = genes,
    group.by = "CellType",
    assay = "Summary",
    layer = "data",
    lib_normalize = FALSE,
    exp_method = "zscore",
    exp_legend_title = if (flip) "Z-score" else "Z-score (mean log expression)",
    flip = flip,
    group_palcolor = brainomics_celltype_colors[type_levels],
    feature_split = factor(marker_spec$Marker_Group, levels = marker_levels),
    feature_split_palcolor = brainomics_celltype_colors[marker_levels],
    heatmap_palette = "Spectral",
    heatmap_border = TRUE,
    heatmap_border_palcolor = "grey",
    heatmap_border_size = 0.5,
    cell_annotation_border = FALSE,
    feature_annotation_border = FALSE,
    add_dot = FALSE,
    nlabel = 0,
    show_row_names = TRUE,
    show_column_names = flip,
    row_names_side = if (flip) "left" else "right",
    row_title = if (flip) "" else rep(" ", length(marker_levels)),
    row_title_side = "left",
    row_title_rot = 0,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    use_raster = FALSE,
    column_title = if (flip) rep(" ", length(marker_levels)) else "Celltype",
    width = if (flip) 5.4 else 2.5,
    height = if (flip) 1.4 else length(genes) * 0.12,
    legend.position = if (flip) "bottom" else "right",
    verbose = FALSE,
    ht_params = list(
      column_gap = grid::unit(if (flip) .3 else 2, "mm"),
      rect_gp = grid::gpar(type = "none"),
      row_names_gp = grid::gpar(
        fontfamily = "Arial",
        fontsize = 8,
        fontface = if (flip) "plain" else "italic"
      ),
      row_title_gp = grid::gpar(fontfamily = "Arial", fontsize = 7),
      column_names_gp = grid::gpar(fontfamily = "Arial", fontsize = if (flip) 8 else 7,
        fontface = if (flip) "italic" else "plain"),
      column_names_max_height = grid::unit(65, "mm"),
      column_title_gp = grid::gpar(fontfamily = "Arial", fontsize = 10),
      cell_fun = function(j, i, x, y, width, height, fill) {
        pct <- if (flip) fraction[j, i] else fraction[i, j]
        if (pct > 0) {
          grid::grid.points(
            x,
            y,
            pch = 21,
            size = grid::unit(3 * pct, "mm"),
            gp = grid::gpar(fill = fill, col = "black", lwd = 0.6)
          )
        }
      }
    )
  )
  expected <- t(scale(t(means)))
  stopifnot(isTRUE(all.equal(
    as.numeric(p$matrix_list[[1]][genes, type_levels]),
    as.numeric(expected),
    tolerance = 1e-7
  )))
  p
}

# Fixed manuscript CellType palette shared by reference and mapping UMAPs.
read_celltype_assignments <- function(cells = NULL) {
  path <- "../../data/BrainOmicsData/integration_25/annotation/celltype_assignments.rds"
  annotation <- as.data.frame(readRDS(path))
  stopifnot(
    identical(names(annotation), c("Cells", "Cluster", "CellType")),
    nrow(annotation) == 2602031L,
    !anyNA(annotation),
    !anyDuplicated(annotation$Cells),
    length(unique(annotation$Cluster)) == 75L,
    length(unique(annotation$CellType)) == 12L
  )
  mapping <- read_tsv("results/annotation/cluster_annotation.tsv")
  index <- match(annotation$Cluster, mapping$Cluster)
  stopifnot(
    nrow(mapping) == 75L,
    !anyDuplicated(mapping$Cluster),
    !anyNA(index),
    identical(annotation$CellType, mapping$CellType[index]),
    identical(
      as.integer(table(factor(annotation$Cluster, levels = mapping$Cluster))),
      as.integer(mapping$Cells)
    )
  )
  if (!is.null(cells)) {
    index <- match(cells, annotation$Cells)
    if (anyNA(index) || anyDuplicated(cells)) {
      stop("Cells do not match the final annotation")
    }
    annotation <- annotation[index, , drop = FALSE]
  }
  annotation
}

validate_celltype_metadata <- function(metadata, cells) {
  annotation <- read_celltype_assignments(cells)
  if (
    any(
      c("Detailed_CellType", "Main_CellType", "CellType_Broad") %in%
        names(metadata)
    )
  ) {
    stop(
      "Metadata contains obsolete cell-type fields; refresh the object first"
    )
  }
  for (field in "CellType") {
    if (
      !field %in% names(metadata) ||
        !identical(as.character(metadata[[field]]), annotation[[field]])
    ) {
      stop(
        "Stored ",
        field,
        " is not the final per-cell annotation; refresh the object first"
      )
    }
  }
  invisible(annotation)
}

# Check the result inventory; annotation freshness is reviewed before analysis/redraw.
validate_annotation_statistics <- function(directory) {
  path <- file.path(directory, "annotation_input.tsv")
  if (!file.exists(path)) {
    stop("Missing annotation result inventory: ", path)
  }
  input <- read_tsv(path)
  if (
    !"File" %in% names(input) ||
      !nrow(input) ||
      anyNA(input$File) ||
      any(!nzchar(input$File)) ||
      anyDuplicated(input$File)
  ) {
    stop("Invalid annotation result inventory: ", path)
  }
  files <- file.path(directory, input$File)
  if (any(!file.exists(files))) {
    stop("Annotation result files are missing: ", directory)
  }
  invisible(TRUE)
}

# Called by an analysis only after its output files have been written successfully.
write_annotation_input <- function(directory, files) {
  stopifnot(
    length(files) > 0L,
    !anyDuplicated(files),
    all(file.exists(file.path(directory, files)))
  )
  write_tsv(
    data.frame(File = files),
    file.path(directory, "annotation_input.tsv")
  )
  validate_annotation_statistics(directory)
}

brainomics_celltype_colors <- c(
  "Astrocytes" = "#D70440",
  "Endothelial cells" = "#5E7987",
  "Mural cells" = "#8C6D31",
  "Perivascular fibroblasts" = "#B07A4A",
  "Excitatory neurons" = "#02AD24",
  "MGE-derived inhibitory neurons" = "#174E7A",
  "CGE-derived inhibitory neurons" = "#4C9BCB",
  "Histaminergic neurons" = "#A24D70",
  "Microglia" = "#006D87",
  "Lymphocytes" = "#7E57C2",
  "Neuroblasts" = "#ED5736",
  "Oligodendrocyte progenitor cells" = "#F9BD10",
  "Oligodendrocyte lineage cells" = "#A6761D",
  "Differentiating oligodendrocytes" = "#E07A2D",
  "Oligodendrocytes" = "#B14B28",
  "Radial glia" = "#8076A3"
)

# Fixed dataset-name mapping shared by all dataset-coloured panels.
dataset_colors <- c(
  AllenM1 = "#377EB8",
  EGAS00001006537 = "#E64B35",
  GSE104276 = "#00A087",
  GSE144136 = "#8E63B0",
  GSE168408 = "#E69F00",
  GSE178175 = "#56B4E9",
  GSE186538 = "#D45087",
  GSE202210 = "#66A61E",
  GSE204683 = "#A65E2E",
  GSE207334 = "#008B8B",
  GSE212606 = "#E78AC3",
  GSE217511 = "#4055A8",
  GSE294786 = "#F28E2B",
  GSE296073 = "#4DAF4A",
  GSE67835 = "#AA4499",
  GSE81475 = "#17BECF",
  GSE97942 = "#B2182B",
  HYPOMAP = "#A6A832",
  Li_et_al_2018 = "#7B4173",
  Ma_et_al_2022 = "#1B7837",
  PRJCA015229 = "#C77CFF",
  ROSMAP = "#CC6677",
  SomaMut = "#0072B2",
  Velmeshev_2023 = "#B8860B",
  Wang_2025 = "#80B1D3"
)

method_levels <- c("Raw", "scVI", "Harmony", "RPCA")
method_colors <- c(
  Raw = "#E41A1C",
  scVI = "#984EA3",
  Harmony = "#4DAF4A",
  RPCA = "#377EB8"
)

# Shared UMAP theme

theme_blank_axis <- function(
  add_coord = TRUE,
  xlen_npc = 0.14,
  ylen_npc = 0.18,
  xlab = "UMAP_1",
  ylab = "UMAP_2",
  lab_size = 12,
  axis_lwd = 1.6,
  arrow_len = 0.02,
  ...
) {
  args1 <- list(
    panel.border = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    legend.background = ggplot2::element_blank(),
    legend.box.margin = ggplot2::margin(0, 0, 0, 0),
    legend.margin = ggplot2::margin(0, 0, 0, 0),
    legend.key = ggplot2::element_rect(
      fill = "transparent",
      color = "transparent"
    ),
    legend.key.size = grid::unit(10, "pt"),
    plot.margin = ggplot2::margin(
      lab_size + 10,
      lab_size + 10,
      lab_size + 10,
      lab_size + 10,
      unit = "points"
    )
  )
  args2 <- as.list(match.call())[-1]
  call_envir <- parent.frame(1)
  args2 <- lapply(args2, function(arg) {
    if (is.symbol(arg) || is.call(arg)) {
      eval(arg, envir = call_envir)
    } else {
      arg
    }
  })
  for (name in names(args2)) {
    args1[[name]] <- args2[[name]]
  }
  args <- args1[names(args1) %in% methods::formalArgs(ggplot2::theme)]
  out <- do.call(ggplot2::theme, args)
  if (isTRUE(add_coord)) {
    grob <- grid::grobTree(
      grid::gList(
        grid::linesGrob(
          x = grid::unit(c(0, xlen_npc), "npc"),
          y = grid::unit(c(0, 0), "npc"),
          arrow = grid::arrow(length = grid::unit(arrow_len, "npc")),
          gp = grid::gpar(lwd = axis_lwd, col = "black")
        ),
        grid::linesGrob(
          x = grid::unit(c(0, 0), "npc"),
          y = grid::unit(c(0, ylen_npc), "npc"),
          arrow = grid::arrow(length = grid::unit(arrow_len, "npc")),
          gp = grid::gpar(lwd = axis_lwd, col = "black")
        ),
        grid::textGrob(
          label = xlab,
          x = grid::unit(0, "npc"),
          y = grid::unit(-0.55, "lines"),
          vjust = 1,
          hjust = 0,
          gp = grid::gpar(fontsize = lab_size)
        ),
        grid::textGrob(
          label = ylab,
          x = grid::unit(-0.75, "lines"),
          y = grid::unit(0, "npc"),
          vjust = 0.5,
          hjust = 0,
          rot = 90,
          gp = grid::gpar(fontsize = lab_size)
        )
      )
    )
    return(list(
      list(ggplot2::annotation_custom(grob)),
      list(thisplot::theme_this() + out),
      list(ggplot2::coord_cartesian(clip = "off"))
    ))
  }
  list(list(thisplot::theme_this() + out))
}

# Vector PDF assembly

# Vector-PDF figure assembly adapted from the panel-first workflow used by
# multiCSN_manuscript.  Each panel is cropped and embedded as PDF, so the
# assembled figure does not rasterize text, legends or heatmaps.

pdf_panel <- function(
  id,
  x,
  y,
  width,
  height = NULL,
  label_x = x - 1.5,
  label_y = y,
  source = NULL,
  label = id
) {
  list(
    id = id,
    x = x,
    y = y,
    width = width,
    height = height,
    label_x = label_x,
    label_y = label_y,
    source = source,
    label = label
  )
}

assemble_pdf_figure <- function(
  stem,
  page,
  panels,
  repo,
  output_dir = file.path(repo, "figures"),
  work_root = tempdir()
) {
  dir.create(work_root, recursive = TRUE, showWarnings = FALSE)
  current_work <- tempfile(paste0(stem, "-"), tmpdir = work_root)
  on.exit(unlink(current_work, recursive = TRUE), add = TRUE)
  cropped_root <- file.path(current_work, "cropped")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cropped_root, recursive = TRUE, showWarnings = FALSE)

  programs <- Sys.which(c(pdfcrop = "pdfcrop", xelatex = "xelatex"))
  tinytex_bin <- Sys.glob(path.expand("~/.local/TinyTeX/bin/*"))
  if (length(tinytex_bin) == 1L) {
    local_programs <- file.path(tinytex_bin, names(programs))
    use_local <- !nzchar(programs) & file.exists(local_programs)
    programs[use_local] <- local_programs[use_local]
  }
  if (any(!nzchar(programs))) {
    stop(
      "PDF assembly requires: ",
      paste(names(programs)[!nzchar(programs)], collapse = ", "),
      call. = FALSE
    )
  }
  run_checked <- function(command, args) {
    output <- system2(command, args = args, stdout = TRUE, stderr = TRUE)
    status <- attr(output, "status")
    if (!is.null(status) && status != 0L) {
      stop(
        "Command failed: ",
        command,
        " ",
        paste(args, collapse = " "),
        "\n",
        paste(output, collapse = "\n"),
        call. = FALSE
      )
    }
    invisible(output)
  }

  panel_ids <- vapply(panels, `[[`, character(1), "id")
  if (anyDuplicated(panel_ids)) {
    stop("Duplicate panel ids in ", stem)
  }
  source_paths <- setNames(
    vapply(
      panels,
      function(current) {
        candidate <- current$source
        if (is.null(candidate)) {
          candidate <- paste0(stem, current$id, ".pdf")
        }
        if (!grepl("^/", candidate)) {
          candidate <- file.path(output_dir, candidate)
        }
        normalizePath(candidate, mustWork = FALSE)
      },
      character(1)
    ),
    panel_ids
  )
  missing <- source_paths[!file.exists(source_paths)]
  if (length(missing)) {
    stop("Missing panel PDFs: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  cropped_paths <- file.path(cropped_root, basename(source_paths))
  names(cropped_paths) <- names(source_paths)
  for (id in names(source_paths)) {
    run_checked(
      programs[["pdfcrop"]],
      c(
        "--margins",
        shQuote("3 3 3 3"),
        shQuote(source_paths[[id]]),
        shQuote(cropped_paths[[id]])
      )
    )
  }

  tex_path <- function(path) {
    normalizePath(path, winslash = "/", mustWork = TRUE)
  }
  nodes <- unlist(lapply(panels, function(current) {
    dimensions <- paste0("width=", current$width, "mm")
    if (!is.null(current$height)) {
      dimensions <- paste0(
        dimensions,
        ",height=",
        current$height,
        "mm,keepaspectratio"
      )
    }
    panel_node <- sprintf(
      paste0(
        "  \\node[anchor=north west,inner sep=0pt] at (%s,%s) ",
        "{\\includegraphics[%s]{\\detokenize{%s}}};"
      ),
      current$x,
      current$y,
      dimensions,
      tex_path(cropped_paths[[current$id]])
    )
    if (is.null(current$label) || !nzchar(current$label)) {
      panel_node
    } else {
      c(
        panel_node,
        sprintf(
          "  \\panellabel{%s}{%s}{%s}",
          current$label_x,
          current$label_y,
          current$label
        )
      )
    }
  }))

  tex_lines <- c(
    "\\documentclass{article}",
    sprintf(
      "\\usepackage[paperwidth=%smm,paperheight=%smm,margin=0mm]{geometry}",
      page[[1L]],
      page[[2L]]
    ),
    "\\usepackage{fontspec}",
    "\\setsansfont{Arial}",
    "\\usepackage{graphicx}",
    "\\usepackage{tikz}",
    "\\pagestyle{empty}",
    "\\setlength{\\parindent}{0pt}",
    "\\newcommand{\\panellabel}[3]{%",
    paste0(
      "  \\node[anchor=north west,inner sep=0pt,font=\\sffamily",
      "\\bfseries\\fontsize{10}{10}\\selectfont] at (#1,#2) {#3};%"
    ),
    "}",
    "\\begin{document}",
    "\\noindent",
    "\\begin{tikzpicture}[x=1mm,y=1mm]",
    sprintf(
      "  \\useasboundingbox (0,0) rectangle (%s,%s);",
      page[[1L]],
      page[[2L]]
    ),
    nodes,
    "\\end{tikzpicture}",
    "\\end{document}"
  )
  tex_file <- file.path(current_work, paste0(stem, ".tex"))
  writeLines(tex_lines, tex_file, useBytes = TRUE)
  run_checked(
    programs[["xelatex"]],
    c(
      "-halt-on-error",
      "-interaction=nonstopmode",
      paste0("-output-directory=", shQuote(current_work)),
      shQuote(tex_file)
    )
  )
  compiled <- file.path(current_work, paste0(stem, ".pdf"))
  output <- file.path(output_dir, paste0(stem, ".pdf"))
  run_checked(
    programs[["pdfcrop"]],
    c(
      "--margins",
      shQuote("6 6 6 6"),
      shQuote(compiled),
      shQuote(output)
    )
  )
  message("Saved assembled ", stem, " to ", output)
  invisible(output)
}

# Vector exports of an assembled figure plus the inverted night variant written by
# thisplot.  The PDF stays the vector master for the manuscript; the SVG is the
# shareable copy.
export_figure_assets <- function(pdf, suffix = "-night") {
  pdf <- normalizePath(pdf, mustWork = TRUE)
  stem <- sub("[.]pdf$", "", pdf)
  programs <- Sys.which(c(pdftocairo = "pdftocairo"))
  if (!nzchar(programs[["pdftocairo"]])) {
    stop("Figure export requires: pdftocairo", call. = FALSE)
  }
  args <- c("-svg", shQuote(pdf), shQuote(paste0(stem, ".svg")))
  output <- system2(programs[["pdftocairo"]], args = args, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop(
      "Command failed: pdftocairo ", paste(args, collapse = " "), "\n",
      paste(output, collapse = "\n"),
      call. = FALSE
    )
  }
  vector <- paste0(stem, ".svg")
  if (!file.exists(vector)) {
    stop("Figure export did not write: ", vector, call. = FALSE)
  }
  night <- unlist(thisplot::invert_figures(vector, suffix = suffix, quiet = TRUE))
  message("Saved figure exports to ", paste(c(vector, night), collapse = ", "))
  invisible(c(vector, night))
}

# Source-label palette shared by held-out projection and centroid heatmap.
brainomics_query_celltype_colors <- c(
  "Neurons" = "#8F6D9B",
  "Neuroblasts" = "#C47A6B",
  "Neuronal intermediate progenitor cells" = "#D6A24D",
  "Radial glia" = "#7B8F57",
  "Glioblasts" = "#4F8B7F",
  "Oligodendrocyte lineage" = "#7197B5",
  "Immune cells" = "#B15F7A",
  "Vascular cells" = "#355C7D",
  "Perivascular fibroblasts" = "#9B886F",
  "Erythrocytes" = "#A8473D",
  "Placodes" = "#B790B2",
  "Neural crest" = "#6F7686"
)
