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
library(thisutils)
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


color_stages <- c(
  colorRampPalette(
    c("#0AA344", "#006D87")
  )(7),
  colorRampPalette(
    c("#2B73AF", "#003D74")
  )(8)
)
names(color_stages) <- paste0("S", 1:15)

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
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.background = element_blank(),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0),
    legend.key = element_rect(fill = "transparent", color = "transparent"),
    legend.key.size = grid::unit(10, "pt"),
    plot.margin = margin(
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

  for (n in names(args2)) {
    args1[[n]] <- args2[[n]]
  }

  args <- args1[names(args1) %in% methods::formalArgs(theme)]
  out <- do.call(theme, args)

  if (isTRUE(add_coord)) {
    g <- grid::grobTree(
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
      list(ggplot2::annotation_custom(g)),
      list(thisplot::theme_this() + out),
      list(ggplot2::coord_cartesian(clip = "off"))
    ))
  }

  list(list(thisplot::theme_this() + out))
}
