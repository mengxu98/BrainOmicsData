#!/usr/bin/env Rscript
# Figure 3: marker evidence, heatmap and source-label concordance.
for (script in c("functions/cluster_marker_balance.R", "functions/group_heatmap.R",
                 "plotting/fig3_annotation_panels.R")) {
  status <- system2("Rscript", c("--vanilla", script))
  if (status != 0L) stop("Figure 3 step failed: ", script)
}
source("functions/export_png.R")
export_pdf_png("figures/fig3.pdf", "figures/fig3.png", 9.45)
