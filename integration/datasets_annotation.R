library(Seurat)
library(ggplot2)
library(patchwork)
library(thisutils)

res_dir <- "../../data/BrainOmicsData/integration/"

objects <- readRDS(file.path(res_dir, "objects_integrated.rds"))

marker_genes <- c(
  # Radial glia
  "PAX6", "VIM", "GLI3",
  # Endothelial cells
  "CLDN5", "PECAM1", "VWF", "FLT1",
  # Inhibitory neurons
  "GAD1", "GAD2", "SLC6A1",
  # Oligodendrocyte progenitor cells (OPCs)
  "PDGFRA", "CSPG4", "OLIG1", "OLIG2", "SOX10",
  # Microglia
  "CX3CR1", "P2RY12", "CSF1R",
  # Neuroblasts
  "STMN2",
  # Excitatory neurons
  "SLC17A7", "CAMK2A", "SATB2",
  # Astrocytes
  "GFAP", "AQP4", "ALDH1L1", "FGFR3", "GJA1",
  # Oligodendrocytes
  "MOG", "MAG", "CLDN11"
)

cluster2celltype <- c(
  "0" = "Oligodendrocytes",
  "1" = "Oligodendrocytes",
  "2" = "Astrocytes",
  "3" = "Oligodendrocytes",
  "4" = "Excitatory neurons",
  "5" = "Neuroblasts",
  "6" = "Microglia",
  "7" = "Oligodendrocyte progenitor cells",
  "8" = "Inhibitory neurons",
  "9" = "Excitatory neurons",
  "10" = "Inhibitory neurons",
  "11" = "Inhibitory neurons",
  "12" = "Excitatory neurons",
  "13" = "Oligodendrocytes",
  "14" = "Astrocytes",
  "15" = "Inhibitory neurons",
  "16" = "Excitatory neurons",
  "17" = "Excitatory neurons",
  "18" = "Inhibitory neurons",
  "19" = "Excitatory neurons",
  "20" = "Oligodendrocytes",
  "21" = "Endothelial cells",
  "22" = "Excitatory neurons",
  "23" = "Endothelial cells",
  "24" = "Neuroblasts",
  "25" = "Inhibitory neurons",
  "26" = "Inhibitory neurons",
  "27" = "Excitatory neurons",
  "28" = "Excitatory neurons",
  "29" = "Radial glia",
  "30" = "Excitatory neurons",
  "31" = "Astrocytes",
  "32" = "Inhibitory neurons",
  "33" = "Microglia",
  "34" = "Inhibitory neurons",
  "35" = "Astrocytes",
  "36" = "Inhibitory neurons",
  "37" = "Oligodendrocyte progenitor cells",
  "38" = "Excitatory neurons",
  "39" = "Excitatory neurons",
  "40" = "Astrocytes",
  "41" = "Radial glia",
  "42" = "Astrocytes",
  "43" = "Excitatory neurons",
  "44" = "Excitatory neurons",
  "45" = "Excitatory neurons",
  "46" = "Oligodendrocyte progenitor cells",
  "47" = "Excitatory neurons",
  "48" = "Excitatory neurons",
  "49" = "Astrocytes",
  "50" = "Excitatory neurons",
  "51" = "Inhibitory neurons",
  "52" = "Excitatory neurons",
  "53" = "Excitatory neurons",
  "54" = "Excitatory neurons",
  "55" = "Radial glia",
  "56" = "Oligodendrocyte progenitor cells",
  "57" = "Excitatory neurons",
  "58" = "Excitatory neurons",
  "59" = "Inhibitory neurons",
  "60" = "Excitatory neurons",
  "61" = "Excitatory neurons",
  "62" = "Excitatory neurons",
  "63" = "Inhibitory neurons",
  "64" = "Inhibitory neurons",
  "65" = "Excitatory neurons",
  "66" = "Astrocytes",
  "67" = "Radial glia",
  "68" = "Excitatory neurons",
  "69" = "Excitatory neurons",
  "70" = "Excitatory neurons",
  "71" = "Excitatory neurons",
  "72" = "Excitatory neurons",
  "73" = "Inhibitory neurons",
  "74" = "Excitatory neurons",
  "75" = "Neuroblasts",
  "76" = "Excitatory neurons",
  "77" = "Radial glia",
  "78" = "Excitatory neurons",
  "79" = "Excitatory neurons",
  "80" = "Oligodendrocytes",
  "81" = "Inhibitory neurons",
  "82" = "Excitatory neurons",
  "83" = "Excitatory neurons",
  "84" = "Excitatory neurons",
  "85" = "Inhibitory neurons",
  "86" = "Inhibitory neurons",
  "87" = "Astrocytes",
  "88" = "Excitatory neurons",
  "89" = "Excitatory neurons",
  "90" = "Excitatory neurons",
  "91" = "Excitatory neurons",
  "92" = "Excitatory neurons",
  "93" = "Inhibitory neurons",
  "94" = "Excitatory neurons",
  "95" = "Oligodendrocytes",
  "96" = "Microglia",
  "97" = "Excitatory neurons",
  "98" = "Excitatory neurons",
  "99" = "Oligodendrocyte progenitor cells",
  "100" = "Excitatory neurons",
  "101" = "Oligodendrocytes",
  "102" = "Excitatory neurons",
  "103" = "Inhibitory neurons",
  "104" = "Oligodendrocytes",
  "105" = "Excitatory neurons",
  "106" = "Excitatory neurons",
  "107" = "Inhibitory neurons",
  "108" = "Excitatory neurons",
  "109" = "Excitatory neurons",
  "110" = "Excitatory neurons",
  "111" = "Neuroblasts",
  "112" = "Inhibitory neurons",
  "113" = "Excitatory neurons",
  "114" = "Excitatory neurons",
  "115" = "Radial glia",
  "116" = "Excitatory neurons",
  "117" = "Neuroblasts",
  "118" = "Oligodendrocytes",
  "119" = "Oligodendrocytes"
)
objects$CellType <- plyr::mapvalues(
  x    = objects$seurat_clusters,
  from = names(cluster2celltype),
  to   = cluster2celltype
)

Idents(objects) <- "CellType"

log_message("Save data...")
saveRDS(objects, file.path(res_dir, "objects_celltypes.rds"))

objects_plot <- objects[intersect(rownames(objects), marker_genes), ]
saveRDS(objects_plot, file.path(res_dir, "objects_celltype_plot.rds"))
