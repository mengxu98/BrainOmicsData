#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1]]
}

object_file <- normalizePath(
  value_after("--object-file", "/work/home/mengxu1310/data/BrainOmicsData/integration/objects_celltypes.rds"),
  mustWork = TRUE
)
out_file <- value_after(
  "--out-file",
  "/work/home/mengxu1310/data/ScienceDB/HumanBrain_sc_snRNAseq_age_resource/04_processed_objects/human_brain_integrated_full.h5ad"
)

if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat is required.")
if (!requireNamespace("SeuratObject", quietly = TRUE)) stop("SeuratObject is required.")
if (!requireNamespace("scop", quietly = TRUE)) stop("scop is required.")
if (!requireNamespace("reticulate", quietly = TRUE)) stop("reticulate is required.")
if (!reticulate::py_module_available("anndata")) stop("Python module anndata is required.")
if (!reticulate::py_module_available("numpy")) stop("Python module numpy is required.")

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

message("Loading Seurat object: ", object_file)
object <- readRDS(object_file)

message("Writing h5ad: ", out_file)
scop::srt_to_h5ad(
  object,
  path = out_file,
  assay_x = "RNA",
  layer_x = "counts",
  reductions = names(object@reductions),
  graphs = character(0),
  neighbors = character(0),
  convert_tools = FALSE,
  convert_misc = FALSE,
  overwrite = TRUE
)

message("h5ad export complete: ", out_file)
