#!/usr/bin/env Rscript
# Formal Table 2: user-approved Scheme 2 styling with verified age coordinates.
suppressPackageStartupMessages({
  library(grid)
  library(grDevices)
})

out <- "figures"
data_dir <- "data/table2_comparison"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
d <- read.delim(file.path(data_dir, "source_data.tsv"), sep = "\t", quote = "",
                stringsAsFactors = FALSE, check.names = FALSE)

# Keep a character column even when every count is exact (no ">" values).
d$Human_Cells_or_Nuclei <- as.character(d$Human_Cells_or_Nuclei)

# Format count helper
format_count <- function(s) {
  value <- format(as.numeric(sub("^>", "", s)), big.mark = ",", scientific = FALSE, trim = TRUE)
  paste0(if (startsWith(s, ">")) ">" else "", value)
}

labels <- read.delim(file.path(data_dir, "verified_plot_fields.tsv"), sep = "\t",
                     quote = "", stringsAsFactors = FALSE, check.names = FALSE)
source_key <- sub(" revised reference$", "", d$Resource)
stopifnot(nrow(labels) == 5L, !anyDuplicated(labels$Resource),
          !anyDuplicated(source_key), setequal(labels$Resource, source_key))
d <- d[match(labels$Resource, source_key), , drop = FALSE]
labels$Cell_count <- d$Human_Cells_or_Nuclei

# Shared, piecewise-linear scale. Prenatal 0–40 PCW uses 10% of the track;
# postnatal 0–18 years uses 45%; 18–104 years uses the remaining 45%.
# Prenatal and postnatal units are kept separate; no gestational duration is
# inferred for individual births. GW = PCW + 2 follows Kim et al.'s convention.
lifespan_position <- function(age, unit) {
  if (is.na(age)) return(NA_real_)
  stopifnot(unit %in% c("PCW", "GW", "years"), is.finite(age), age >= 0)
  if (unit %in% c("PCW", "GW")) {
    pcw <- age - if (unit == "GW") 2 else 0
    stopifnot(pcw >= 0, pcw <= 40)
    return(0.10 * pcw / 40)
  }
  stopifnot(age <= 104)
  if (age <= 18) return(0.10 + 0.45 * age / 18)
  0.55 + 0.45 * (age - 18) / (104 - 18)
}
age_start <- mapply(lifespan_position, labels$Age_min, labels$Age_min_unit)
age_end <- mapply(lifespan_position, labels$Age_max, labels$Age_max_unit)
stopifnot(identical(is.na(age_start), is.na(age_end)),
          all(age_start <= age_end, na.rm = TRUE),
          lifespan_position(5, "PCW") == lifespan_position(7, "GW"),
          isTRUE(all.equal(lifespan_position(18, "years"), 0.55)),
          isTRUE(all.equal(lifespan_position(104, "years"), 1)))
write.table(labels[c("Resource", "Cell_count", "Donors", "Age", "Regions", "Shared")],
            file.path(data_dir, "display_data.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)
write.table(data.frame(Resource = labels$Resource, Age_min = labels$Age_min,
  Age_min_unit = labels$Age_min_unit, Age_max = labels$Age_max,
  Age_max_unit = labels$Age_max_unit, Start = age_start, End = age_end,
  Coverage_scope = labels$Coverage_scope), file.path(data_dir, "age_geometry.tsv"),
  sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")

# Colors
ink        <- "#1E293B"  # Slate 800 - dark, crisp
ink_focal  <- "#0F4C81"  # Deep classic blue for BrainOmicsData title
muted      <- "#64748B"  # Slate 500 for secondary text
nr_color   <- "#94A3B8"  # Slate 400 for NR/missing
bar_track  <- "#EAEFF2"  # Subtle track background
bar_focal  <- "#1D6399"  # High-quality journal blue
bar_other  <- "#94A3B8"  # Balanced slate for other resources
accent_bg  <- "#F0F5FA"  # Soft elegant highlight for row 1
rule_heavy <- "#1E293B"  # Top rule
rule_mid   <- "#475569"  # Header divider
rule_light <- "#E2E8F0"  # Subtle row divider
grid_col   <- "#E2E8F0"  # Dotted vertical guideline

# ==============================================================================
# Variant 1: refined academic hybrid layout
# ==============================================================================
plot_refined <- function() {
  width_mm  <- 183
  height_mm <- 94
  family    <- "Arial"
  
  grid.newpage()
  
  # Canvas background
  grid.rect(x = unit(0, "mm"), y = unit(height_mm, "mm"),
            width = unit(width_mm, "mm"), height = unit(height_mm, "mm"),
            just = c("left", "top"), gp = gpar(fill = "white", col = NA))
  
  # Title
  grid.text("Human brain transcriptomic resources",
            x = unit(6, "mm"), y = unit(height_mm - 7.5, "mm"),
            just = c("left", "centre"),
            gp = gpar(fontfamily = family, fontsize = 10.5, fontface = "bold", col = ink))
  
  # Top heavy rule
  grid.lines(x = unit(c(6, 177), "mm"), y = unit(c(height_mm - 13, height_mm - 13), "mm"),
             gp = gpar(col = rule_heavy, lwd = 1.0))
  
  # Column layout:
  # Col 1: Resource (6 to 36 mm)
  # Col 2: Cells + Bar (36 to 88 mm), bar_x0 = 38, bar_w = 46
  # Col 3: Donors (88 to 102 mm)
  # Col 4: Age range (102 to 134 mm)
  # Col 5: Brain regions (134 to 160 mm)
  # Col 6: Shared sources (160 to 177 mm)
  
  # Headers
  grid.text("Resource", x = unit(8, "mm"), y = unit(height_mm - 18, "mm"),
            just = c("left", "centre"), gp = gpar(fontfamily = family, fontsize = 8.0, fontface = "bold", col = ink))
  
  grid.text("Human cells / nuclei", x = unit(38, "mm"), y = unit(height_mm - 16.5, "mm"),
            just = c("left", "centre"), gp = gpar(fontfamily = family, fontsize = 8.0, fontface = "bold", col = ink))
  
  bar_x0 <- 38
  bar_w  <- 46
  ticks  <- c(0, 1e6, 2e6, 3e6)
  tick_labels <- c("0", "1M", "2M", "3M")
  for (t in seq_along(ticks)) {
    tx <- bar_x0 + (ticks[t] / 3e6) * bar_w
    just_x <- if (t == 1) "left" else if (t == 4) "right" else "centre"
    grid.text(tick_labels[t], x = unit(tx, "mm"), y = unit(height_mm - 21.2, "mm"),
              just = c(just_x, "centre"), gp = gpar(fontfamily = family, fontsize = 6.8, col = muted))
  }
  
  grid.text("Donors", x = unit(100, "mm"), y = unit(height_mm - 18, "mm"),
            just = c("right", "centre"), gp = gpar(fontfamily = family, fontsize = 8.0, fontface = "bold", col = ink))
  
  grid.text("Age range", x = unit(105, "mm"), y = unit(height_mm - 18, "mm"),
            just = c("left", "centre"), gp = gpar(fontfamily = family, fontsize = 8.0, fontface = "bold", col = ink))
  
  grid.text("Brain regions", x = unit(135, "mm"), y = unit(height_mm - 18, "mm"),
            just = c("left", "centre"), gp = gpar(fontfamily = family, fontsize = 8.0, fontface = "bold", col = ink))
  
  grid.text("Shared\nsources", x = unit(175, "mm"), y = unit(height_mm - 18, "mm"),
            just = c("right", "centre"), gp = gpar(fontfamily = family, fontsize = 7.5, fontface = "bold", col = ink, lineheight = 0.95))
  
  # Header bottom rule
  grid.lines(x = unit(c(6, 177), "mm"), y = unit(c(height_mm - 24, height_mm - 24), "mm"),
             gp = gpar(col = rule_mid, lwd = 0.6))
  
  # Vertical subtle grid lines across rows at 1M, 2M, 3M
  for (t in 2:4) {
    tx <- bar_x0 + (ticks[t] / 3e6) * bar_w
    grid.lines(x = unit(rep(tx, 2), "mm"),
               y = unit(c(height_mm - 24, height_mm - 81.5), "mm"),
               gp = gpar(col = grid_col, lwd = 0.5, lty = "dotted"))
  }
  
  # Rows
  row_h <- 11.5
  y0    <- 24
  for (i in 1:nrow(labels)) {
    top <- y0 + (i - 1) * row_h
    mid <- top + row_h / 2
    focal <- (i == 1L)
    
    if (focal) {
      grid.rect(x = unit(6, "mm"), y = unit(height_mm - top, "mm"),
                width = unit(171, "mm"), height = unit(row_h, "mm"),
                just = c("left", "top"), gp = gpar(fill = accent_bg, col = NA))
      grid.rect(x = unit(6, "mm"), y = unit(height_mm - top, "mm"),
                width = unit(1.2, "mm"), height = unit(row_h, "mm"),
                just = c("left", "top"), gp = gpar(fill = bar_focal, col = NA))
    }
    
    # Col 1: Resource name
    grid.text(labels$Resource[i], x = unit(8, "mm"), y = unit(height_mm - mid, "mm"),
              just = c("left", "centre"),
              gp = gpar(fontfamily = family, fontsize = 8.2,
                        fontface = if (focal) "bold" else "plain",
                        col = if (focal) ink_focal else ink))
    
    # Col 2: Cells - text + mini bar
    count_raw <- d$Human_Cells_or_Nuclei[i]
    count_val <- as.numeric(sub("^>", "", count_raw))
    is_lower  <- startsWith(count_raw, ">")
    
    grid.text(format_count(count_raw), x = unit(bar_x0, "mm"), y = unit(height_mm - (mid - 2.1), "mm"),
              just = c("left", "centre"),
              gp = gpar(fontfamily = family, fontsize = 7.8,
                        fontface = if (focal) "bold" else "plain",
                        col = if (focal) ink_focal else ink))
    
    bar_y <- mid + 2.5
    bar_h <- 2.2
    grid.rect(x = unit(bar_x0, "mm"), y = unit(height_mm - (bar_y - bar_h/2), "mm"),
              width = unit(bar_w, "mm"), height = unit(bar_h, "mm"),
              just = c("left", "top"), gp = gpar(fill = bar_track, col = NA))
    
    val_w <- min((count_val / 3e6) * bar_w, bar_w)
    grid.rect(x = unit(bar_x0, "mm"), y = unit(height_mm - (bar_y - bar_h/2), "mm"),
              width = unit(val_w, "mm"), height = unit(bar_h, "mm"),
              just = c("left", "top"),
              gp = gpar(fill = if (focal) bar_focal else bar_other, col = NA))
    
    if (is_lower) {
      grid.lines(x = unit(c(bar_x0 + val_w + 0.3, bar_x0 + val_w + 4.2), "mm"),
                 y = unit(c(height_mm - bar_y, height_mm - bar_y), "mm"),
                 arrow = arrow(length = unit(1.1, "mm"), type = "open"),
                 gp = gpar(col = "#475569", lwd = 0.85))
    }
    
    # Col 3: Donors (right aligned)
    is_nr_donor <- labels$Donors[i] == "NR"
    grid.text(labels$Donors[i], x = unit(100, "mm"), y = unit(height_mm - mid, "mm"),
              just = c("right", "centre"),
              gp = gpar(fontfamily = family, fontsize = 8.0,
                        fontface = if (focal) "bold" else "plain",
                        col = if (is_nr_donor) nr_color else if (focal) ink_focal else ink))
    
    # Col 4: Age range
    is_nr_age <- labels$Age[i] == "NR"
    grid.text(labels$Age[i], x = unit(105, "mm"), y = unit(height_mm - mid, "mm"),
              just = c("left", "centre"),
              gp = gpar(fontfamily = family, fontsize = 8.0,
                        col = if (is_nr_age) nr_color else ink))
    
    # Col 5: Brain regions
    grid.text(labels$Regions[i], x = unit(135, "mm"), y = unit(height_mm - mid, "mm"),
              just = c("left", "centre"),
              gp = gpar(fontfamily = family, fontsize = 7.8, col = ink))
    
    # Col 6: Shared sources (right aligned)
    grid.text(labels$Shared[i], x = unit(175, "mm"), y = unit(height_mm - mid, "mm"),
              just = c("right", "centre"),
              gp = gpar(fontfamily = family, fontsize = 8.0,
                        col = if (focal) ink_focal else ink))
    
    # Row dividing line
    grid.lines(x = unit(c(6, 177), "mm"),
               y = unit(c(height_mm - (top + row_h), height_mm - (top + row_h)), "mm"),
               gp = gpar(col = if (i == nrow(labels)) rule_heavy else rule_light,
                         lwd = if (i == nrow(labels)) 0.8 else 0.4))
  }
  
  # Footnotes
  notes <- c(
    "PCW, post-conception weeks; GW, gestational weeks; yr, years; NR, overall value not reported or reliably resolved.",
    "* Human + mouse. † Original 10x release. ‡ AllenM1 cohort link. § Public metadata subset."
  )
  grid.text(notes[1], x = unit(6, "mm"), y = unit(height_mm - 86, "mm"),
            just = c("left", "centre"), gp = gpar(fontfamily = family, fontsize = 6.8, col = muted))
  grid.text(notes[2], x = unit(6, "mm"), y = unit(height_mm - 89.5, "mm"),
            just = c("left", "centre"), gp = gpar(fontfamily = family, fontsize = 6.8, col = muted))
}

# ==============================================================================
# Variant 2: infographic table with lifespan timeline
# ==============================================================================
plot_timeline <- function() {
  width_mm  <- 183
  height_mm <- 96
  family    <- "Arial"
  
  grid.newpage()
  
  # Canvas background
  grid.rect(x = unit(0, "mm"), y = unit(height_mm, "mm"),
            width = unit(width_mm, "mm"), height = unit(height_mm, "mm"),
            just = c("left", "top"), gp = gpar(fill = "white", col = NA))
  
  # Title
  grid.text("Human brain transcriptomic resources: scale and cross-lifespan coverage",
            x = unit(6, "mm"), y = unit(height_mm - 7.5, "mm"),
            just = c("left", "centre"),
            gp = gpar(fontfamily = family, fontsize = 10.2, fontface = "bold", col = ink))
  
  # Top heavy rule
  grid.lines(x = unit(c(6, 177), "mm"), y = unit(c(height_mm - 13, height_mm - 13), "mm"),
             gp = gpar(col = rule_heavy, lwd = 1.0))
  
  # Layout columns:
  # Col 1: Resource (6 to 35 mm) -> 29 mm
  # Col 2: Human cells (35 to 76 mm) -> 41 mm (bar_x0 = 36, bar_w = 37 mm)
  # Col 3: Donors (76 to 89 mm) -> 13 mm (right aligned at 87 mm)
  # Col 4: Lifespan timeline (89 to 133 mm) -> 44 mm (age_x0 = 91, age_w = 39 mm)
  # Col 5: Brain regions (133 to 159 mm) -> 26 mm (left aligned at 134 mm)
  # Col 6: Shared (159 to 177 mm) -> 18 mm (right aligned at 175 mm)
  
  # Headers
  grid.text("Resource", x = unit(8, "mm"), y = unit(height_mm - 18, "mm"),
            just = c("left", "centre"), gp = gpar(fontfamily = family, fontsize = 8.0, fontface = "bold", col = ink))
  
  # Col 2: Cells header + axis ticks
  grid.text("Human cells / nuclei", x = unit(36, "mm"), y = unit(height_mm - 16.5, "mm"),
            just = c("left", "centre"), gp = gpar(fontfamily = family, fontsize = 8.0, fontface = "bold", col = ink))
  cell_x0 <- 36
  cell_w  <- 37
  cell_ticks <- c(0, 1e6, 2e6, 3e6)
  cell_tick_lbls <- c("0", "1M", "2M", "3M")
  for (t in seq_along(cell_ticks)) {
    tx <- cell_x0 + (cell_ticks[t] / 3e6) * cell_w
    just_x <- if (t == 1) "left" else if (t == 4) "right" else "centre"
    grid.text(cell_tick_lbls[t], x = unit(tx, "mm"), y = unit(height_mm - 21.2, "mm"),
              just = c(just_x, "centre"), gp = gpar(fontfamily = family, fontsize = 6.6, col = muted))
  }
  
  # Col 3: Donors
  grid.text("Donors", x = unit(87, "mm"), y = unit(height_mm - 18, "mm"),
            just = c("right", "centre"), gp = gpar(fontfamily = family, fontsize = 8.0, fontface = "bold", col = ink))
  
  # Col 4: Age range & timeline header + milestones
  grid.text("Lifespan coverage", x = unit(91, "mm"), y = unit(height_mm - 16.5, "mm"),
            just = c("left", "centre"), gp = gpar(fontfamily = family, fontsize = 8.0, fontface = "bold", col = ink))
  age_x0 <- 91
  age_w  <- 39
  age_ticks <- c(0, lifespan_position(18, "years"), lifespan_position(100, "years"))
  age_tick_lbls <- c("Fetal", "Adult", "100 yr")
  for (t in seq_along(age_ticks)) {
    ax <- age_x0 + age_ticks[t] * age_w
    just_ax <- if (t == 1) "left" else if (t == 3) "right" else "centre"
    grid.text(age_tick_lbls[t], x = unit(ax, "mm"), y = unit(height_mm - 21.2, "mm"),
              just = c(just_ax, "centre"),
              gp = gpar(fontfamily = family, fontsize = 6.5, col = muted))
  }
  
  # Col 5: Brain regions
  grid.text("Brain regions", x = unit(134, "mm"), y = unit(height_mm - 18, "mm"),
            just = c("left", "centre"), gp = gpar(fontfamily = family, fontsize = 8.0, fontface = "bold", col = ink))
  
  # Col 6: Shared sources
  grid.text("Shared\nsources", x = unit(175, "mm"), y = unit(height_mm - 18, "mm"),
            just = c("right", "centre"), gp = gpar(fontfamily = family, fontsize = 7.5, fontface = "bold", col = ink, lineheight = 0.95))
  
  # Header bottom rule
  grid.lines(x = unit(c(6, 177), "mm"), y = unit(c(height_mm - 24, height_mm - 24), "mm"),
             gp = gpar(col = rule_mid, lwd = 0.6))
  
  # Vertical subtle grid lines across rows at 1M, 2M, 3M
  for (t in 2:4) {
    tx <- cell_x0 + (cell_ticks[t] / 3e6) * cell_w
    grid.lines(x = unit(rep(tx, 2), "mm"),
               y = unit(c(height_mm - 24, height_mm - 82.5), "mm"),
               gp = gpar(col = grid_col, lwd = 0.5, lty = "dotted"))
  }
  
  # Rows
  row_h <- 11.7
  y0    <- 24
  for (i in 1:nrow(labels)) {
    top <- y0 + (i - 1) * row_h
    mid <- top + row_h / 2
    focal <- (i == 1L)
    
    if (focal) {
      grid.rect(x = unit(6, "mm"), y = unit(height_mm - top, "mm"),
                width = unit(171, "mm"), height = unit(row_h, "mm"),
                just = c("left", "top"), gp = gpar(fill = accent_bg, col = NA))
      grid.rect(x = unit(6, "mm"), y = unit(height_mm - top, "mm"),
                width = unit(1.2, "mm"), height = unit(row_h, "mm"),
                just = c("left", "top"), gp = gpar(fill = bar_focal, col = NA))
    }
    
    # Col 1: Resource
    grid.text(labels$Resource[i], x = unit(8, "mm"), y = unit(height_mm - mid, "mm"),
              just = c("left", "centre"),
              gp = gpar(fontfamily = family, fontsize = 8.2,
                        fontface = if (focal) "bold" else "plain",
                        col = if (focal) ink_focal else ink))
    
    # Col 2: Cells
    count_raw <- d$Human_Cells_or_Nuclei[i]
    count_val <- as.numeric(sub("^>", "", count_raw))
    is_lower  <- startsWith(count_raw, ">")
    
    grid.text(format_count(count_raw), x = unit(cell_x0, "mm"), y = unit(height_mm - (mid - 2.1), "mm"),
              just = c("left", "centre"),
              gp = gpar(fontfamily = family, fontsize = 7.8,
                        fontface = if (focal) "bold" else "plain",
                        col = if (focal) ink_focal else ink))
    
    bar_y <- mid + 2.5
    bar_h <- 2.2
    grid.rect(x = unit(cell_x0, "mm"), y = unit(height_mm - (bar_y - bar_h/2), "mm"),
              width = unit(cell_w, "mm"), height = unit(bar_h, "mm"),
              just = c("left", "top"), gp = gpar(fill = bar_track, col = NA))
    
    val_w <- min((count_val / 3e6) * cell_w, cell_w)
    grid.rect(x = unit(cell_x0, "mm"), y = unit(height_mm - (bar_y - bar_h/2), "mm"),
              width = unit(val_w, "mm"), height = unit(bar_h, "mm"),
              just = c("left", "top"),
              gp = gpar(fill = if (focal) bar_focal else bar_other, col = NA))
    
    if (is_lower) {
      grid.lines(x = unit(c(cell_x0 + val_w + 0.3, cell_x0 + val_w + 3.8), "mm"),
                 y = unit(c(height_mm - bar_y, height_mm - bar_y), "mm"),
                 arrow = arrow(length = unit(1.1, "mm"), type = "open"),
                 gp = gpar(col = "#475569", lwd = 0.85))
    }
    
    # Col 3: Donors
    is_nr_donor <- labels$Donors[i] == "NR"
    grid.text(labels$Donors[i], x = unit(87, "mm"), y = unit(height_mm - mid, "mm"),
              just = c("right", "centre"),
              gp = gpar(fontfamily = family, fontsize = 8.0,
                        fontface = if (focal) "bold" else "plain",
                        col = if (is_nr_donor) nr_color else if (focal) ink_focal else ink))
    
    # Col 4: Age coverage & Lifespan timeline
    grid.text(labels$Age[i], x = unit(age_x0, "mm"), y = unit(height_mm - (mid - 2.1), "mm"),
              just = c("left", "centre"),
              gp = gpar(fontfamily = family, fontsize = 7.5,
                        col = if (labels$Age[i] == "NR") nr_color else ink))
    
    grid.rect(x = unit(age_x0, "mm"), y = unit(height_mm - (bar_y - bar_h/2), "mm"),
              width = unit(age_w, "mm"), height = unit(bar_h, "mm"),
              just = c("left", "top"), gp = gpar(fill = bar_track, col = NA))
    
    span_coords <- if (is.na(age_start[i])) NULL else c(age_start[i], age_end[i])
    
    if (!is.null(span_coords)) {
      sx <- age_x0 + span_coords[1] * age_w
      sw <- (span_coords[2] - span_coords[1]) * age_w
      grid.rect(x = unit(sx, "mm"), y = unit(height_mm - (bar_y - bar_h/2), "mm"),
                width = unit(sw, "mm"), height = unit(bar_h, "mm"),
                just = c("left", "top"),
                gp = gpar(fill = if (focal) bar_focal else "#78909C", col = NA))
    } else {
      grid.text("—", x = unit(age_x0 + age_w/2, "mm"), y = unit(height_mm - bar_y, "mm"),
                just = c("centre", "centre"), gp = gpar(fontfamily = family, fontsize = 7.0, col = nr_color))
    }
    
    # Col 5: Brain regions
    grid.text(labels$Regions[i], x = unit(134, "mm"), y = unit(height_mm - mid, "mm"),
              just = c("left", "centre"),
              gp = gpar(fontfamily = family, fontsize = 7.8, col = ink))
    
    # Col 6: Shared sources
    grid.text(labels$Shared[i], x = unit(175, "mm"), y = unit(height_mm - mid, "mm"),
              just = c("right", "centre"),
              gp = gpar(fontfamily = family, fontsize = 8.0,
                        col = if (focal) ink_focal else ink))
    
    # Row dividing line
    grid.lines(x = unit(c(6, 177), "mm"),
               y = unit(c(height_mm - (top + row_h), height_mm - (top + row_h)), "mm"),
               gp = gpar(col = if (i == nrow(labels)) rule_heavy else rule_light,
                         lwd = if (i == nrow(labels)) 0.8 else 0.4))
  }
  
  # Footnotes
  notes <- c(
    "PCW, post-conception weeks; GW, gestational weeks; yr, years; NR, overall value not reported or reliably resolved.",
    "* Human + mouse. † Original 10x release. ‡ AllenM1 cohort link. § Public metadata subset."
  )
  grid.text(notes[1], x = unit(6, "mm"), y = unit(height_mm - 87.5, "mm"),
            just = c("left", "centre"), gp = gpar(fontfamily = family, fontsize = 6.8, col = muted))
  grid.text(notes[2], x = unit(6, "mm"), y = unit(height_mm - 91, "mm"),
            just = c("left", "centre"), gp = gpar(fontfamily = family, fontsize = 6.8, col = muted))
}

# Formal export: same approved drawing, named for its permanent location.
grDevices::cairo_pdf(file.path(out, "table2_comparison.pdf"), width = 183/25.4, height = 96/25.4, family = "Arial")
plot_timeline()
grDevices::dev.off()

svglite::svglite(file.path(out, "table2_comparison.svg"), width = 183/25.4, height = 96/25.4, bg = "white")
plot_timeline()
grDevices::dev.off()

cat("Formal Table 2 PDF and SVG generated; display_data.tsv synchronized.\n")
