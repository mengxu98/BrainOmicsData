#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1]]
}

repo_dir <- normalizePath(value_after("--repo-dir", "."), mustWork = TRUE)
object_file <- normalizePath(
  value_after("--object-file", "../../data/BrainOmicsData/integration/objects_celltype_plot.rds"),
  mustWork = TRUE
)
out_dir <- value_after("--out-dir", "../../data/BrainOmicsData/ScienceDB")
reduction_pca <- value_after("--pca-reduction", "integrated.rpca")
reduction_umap <- value_after("--umap-reduction", "umap.rpca")
reduction_unintegrated_umap <- value_after("--unintegrated-umap-reduction", "umap.unintegrated")
lisi_file <- value_after("--lisi-file", file.path(repo_dir, "results", "lisi_results.rds"))

source(file.path(repo_dir, "functions", "metadata_schema.R"))
source(file.path(repo_dir, "functions", "sample_schema.R"))

dirs <- list(
  metadata = file.path(out_dir, "metadata"),
  embeddings = file.path(out_dir, "embeddings"),
  validation = file.path(out_dir, "validation"),
  objects = file.path(out_dir, "objects"),
  scripts = file.path(out_dir, "scripts"),
  provenance = file.path(out_dir, "provenance")
)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

write_tsv <- function(x, file, gzip = grepl("\\.gz$", file)) {
  con <- if (gzip) gzfile(file, "wt") else file(file, "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "NA")
}

message("Loading metadata object: ", object_file)
object <- readRDS(object_file)
if (!inherits(object, "Seurat")) stop("object is not a Seurat object")
meta <- object@meta.data
meta[] <- lapply(meta, as.character)
meta <- add_metadata_schema(meta)
meta <- add_sample_schema(meta)
cells <- rownames(meta)

metadata <- data.frame(
  Cells = meta$Cells,
  Dataset = meta$Dataset,
  Technology = meta$sequencing_technology,
  Sequence = meta$sequencing_modality_standardized,
  Sample = meta$Sample_ID,
  Sample_ID = meta$Sample_ID,
  Donor_ID = meta$Donor_ID,
  Library_ID = meta$Library_ID,
  BrainRegion = meta$BrainRegion,
  Age = meta$age_raw,
  Sex = meta$sex_standardized,
  RNA_snn_res.1 = if ("RNA_snn_res.1" %in% names(meta)) meta$RNA_snn_res.1 else NA_character_,
  seurat_clusters = if ("seurat_clusters" %in% names(meta)) meta$seurat_clusters else NA_character_,
  CellType = if ("CellType" %in% names(meta)) meta$CellType else NA_character_,
  AgeIntervalID = meta$AgeIntervalID,
  AgeInterval = meta$AgeInterval,
  AgeRange = meta$AgeRange,
  percent_mito = if ("percent.mt" %in% names(meta)) meta$percent.mt else NA_character_,
  n_counts = if ("nCount_RNA" %in% names(meta)) meta$nCount_RNA else NA_character_,
  n_genes = if ("nFeature_RNA" %in% names(meta)) meta$nFeature_RNA else NA_character_,
  stringsAsFactors = FALSE
)
write_tsv(metadata, file.path(dirs$metadata, "metadata.tsv.gz"))
write_tsv(sample_schema_audit(meta), file.path(dirs$metadata, "sample_schema_audit.tsv"))

object_meta <- meta[, setdiff(names(meta), c("Original_Sample", "Original_Sample_ID", "sample_schema_rule")), drop = FALSE]
object@meta.data <- object_meta
saveRDS(object, file.path(dirs$objects, "objects_celltype_plot.rds"), compress = TRUE)

dataset_summary <- aggregate(Cells ~ Dataset + Technology + Sequence, metadata, length)
names(dataset_summary)[4] <- "cell_count"
donor_counts <- aggregate(Donor_ID ~ Dataset, metadata, function(x) length(unique(x)))
names(donor_counts)[2] <- "donor_count"
sample_counts <- aggregate(Sample_ID ~ Dataset, metadata, function(x) length(unique(x)))
names(sample_counts)[2] <- "biological_sample_count"
library_counts <- aggregate(Library_ID ~ Dataset, metadata, function(x) length(unique(x)))
names(library_counts)[2] <- "library_count"
dataset_summary <- merge(dataset_summary, donor_counts, by = "Dataset", all.x = TRUE)
dataset_summary <- merge(dataset_summary, sample_counts, by = "Dataset", all.x = TRUE)
dataset_summary <- merge(dataset_summary, library_counts, by = "Dataset", all.x = TRUE)
write_tsv(dataset_summary[order(dataset_summary$Dataset), ], file.path(dirs$metadata, "dataset_summary.tsv"))

write_reduction <- function(reduction, file) {
  if (!reduction %in% names(object@reductions)) {
    warning("missing reduction: ", reduction)
    return(invisible(FALSE))
  }
  emb <- as.data.frame(object@reductions[[reduction]]@cell.embeddings)
  emb <- emb[cells, , drop = FALSE]
  emb <- cbind(cell_id = rownames(emb), emb)
  write_tsv(emb, file)
  TRUE
}
unlink(file.path(dirs$embeddings, c("pca.tsv.gz", "umap.tsv.gz")))
write_reduction(reduction_pca, file.path(dirs$embeddings, "integrated_pca.tsv.gz"))
write_reduction(reduction_umap, file.path(dirs$embeddings, "integrated_umap.tsv.gz"))
write_reduction(reduction_unintegrated_umap, file.path(dirs$embeddings, "unintegrated_umap.tsv.gz"))

if (file.exists(lisi_file)) {
  lisi <- readRDS(lisi_file)
  lisi <- as.data.frame(lisi)
  lisi <- lisi[cells, , drop = FALSE]
  lisi <- cbind(cell_id = rownames(lisi), lisi)
  write_tsv(lisi, file.path(dirs$validation, "lisi.tsv.gz"))
}

read_seurat <- c(
  "#!/usr/bin/env Rscript",
  "",
  "suppressPackageStartupMessages({",
  "  library(Matrix)",
  "  library(Seurat)",
  "})",
  "",
  "args <- commandArgs(trailingOnly = TRUE)",
  "cmd <- commandArgs(FALSE)",
  "file_arg <- grep('^--file=', cmd, value = TRUE)",
  "script_file <- if (length(file_arg) > 0) normalizePath(sub('^--file=', '', file_arg[[1]])) else NA_character_",
  "package_dir <- if (length(args) >= 1) args[[1]] else if (!is.na(script_file)) dirname(dirname(script_file)) else getwd()",
  "package_dir <- normalizePath(package_dir, mustWork = TRUE)",
  "",
  "expr_dir <- file.path(package_dir, 'expression')",
  "metadata_file <- file.path(package_dir, 'metadata', 'metadata.tsv.gz')",
  "integrated_umap_file <- file.path(package_dir, 'embeddings', 'integrated_umap.tsv.gz')",
  "integrated_pca_file <- file.path(package_dir, 'embeddings', 'integrated_pca.tsv.gz')",
  "unintegrated_umap_file <- file.path(package_dir, 'embeddings', 'unintegrated_umap.tsv.gz')",
  "",
  "counts <- ReadMtx(",
  "  mtx = file.path(expr_dir, 'matrix.mtx.gz'),",
  "  features = file.path(expr_dir, 'features.tsv.gz'),",
  "  cells = file.path(expr_dir, 'barcodes.tsv.gz'),",
  "  feature.column = 2",
  ")",
  "metadata <- read.delim(metadata_file, sep = '\\t', stringsAsFactors = FALSE, check.names = FALSE)",
  "rownames(metadata) <- metadata$Cells",
  "metadata <- metadata[colnames(counts), , drop = FALSE]",
  "obj <- CreateSeuratObject(counts = counts, meta.data = metadata, assay = 'RNA')",
  "",
  "add_reduction <- function(obj, file, key, name) {",
  "  if (!file.exists(file)) return(obj)",
  "  emb <- read.delim(file, sep = '\\t', stringsAsFactors = FALSE, check.names = FALSE)",
  "  rownames(emb) <- emb$cell_id",
  "  emb$cell_id <- NULL",
  "  emb <- as.matrix(emb[colnames(obj), , drop = FALSE])",
  "  colnames(emb) <- paste0(key, seq_len(ncol(emb)))",
  "  obj[[name]] <- CreateDimReducObject(embeddings = emb, key = key, assay = DefaultAssay(obj))",
  "  obj",
  "}",
  "obj <- add_reduction(obj, integrated_pca_file, 'IRPCA_', 'integrated_pca')",
  "obj <- add_reduction(obj, integrated_umap_file, 'UMAP_', 'integrated_umap')",
  "obj <- add_reduction(obj, unintegrated_umap_file, 'RAWUMAP_', 'unintegrated_umap')",
  "obj"
)
writeLines(read_seurat, file.path(dirs$scripts, "read_seurat.R"))

read_h5ad <- c(
  "#!/usr/bin/env python3",
  "from pathlib import Path",
  "import sys",
  "import pandas as pd",
  "import scanpy as sc",
  "",
  "package_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parents[1]",
  "expr_dir = package_dir / 'expression'",
  "adata = sc.read_10x_mtx(expr_dir, var_names='gene_symbols', make_unique=True)",
  "metadata = pd.read_csv(package_dir / 'metadata' / 'metadata.tsv.gz', sep='\\t')",
  "metadata = metadata.set_index('Cells').loc[adata.obs_names]",
  "adata.obs = metadata",
  "",
  "def add_embedding(name, key):",
  "    path = package_dir / 'embeddings' / f'{name}.tsv.gz'",
  "    if not path.exists():",
  "        return",
  "    emb = pd.read_csv(path, sep='\\t').set_index('cell_id').loc[adata.obs_names]",
  "    adata.obsm[key] = emb.to_numpy()",
  "",
  "add_embedding('integrated_pca', 'X_integrated_pca')",
  "add_embedding('integrated_umap', 'X_integrated_umap')",
  "add_embedding('unintegrated_umap', 'X_unintegrated_umap')",
  "if len(sys.argv) > 2:",
  "    adata.write_h5ad(sys.argv[2])",
  "else:",
  "    print(adata)",
  ""
)
writeLines(read_h5ad, file.path(dirs$scripts, "read_h5ad.py"))

reference_audit <- file.path(repo_dir, "manuscript", "revised_assets", "reference_audit.tsv")
if (file.exists(reference_audit)) {
  file.copy(reference_audit, file.path(dirs$provenance, "references.tsv"), overwrite = TRUE)
}

readme <- c(
  "# An integrated single-cell and single-nucleus transcriptomic dataset of the human brain across age intervals",
  "",
  "This package provides a compact, directly reusable release of an integrated human brain scRNA-seq and snRNA-seq dataset across age intervals.",
  "",
  "## Contents",
  "",
  "- `expression/matrix.mtx.gz`, `expression/features.tsv.gz`, `expression/barcodes.tsv.gz`: 10X-compatible count matrix.",
  "- `metadata/metadata.tsv.gz`: minimal cell-level metadata using corrected harmonized labels.",
  "- `metadata/sample_schema_audit.tsv`: donor, biological-sample and library-count audit by source dataset.",
  "- `embeddings/integrated_pca.tsv.gz` and `embeddings/integrated_umap.tsv.gz`: RPCA-integrated PCA and UMAP coordinates.",
  "- `embeddings/unintegrated_umap.tsv.gz`: pre-integration UMAP coordinates used for integration comparison.",
  "- `validation/lisi.tsv.gz`: per-cell dataset-label LISI values for raw and RPCA embeddings.",
  "- `objects/objects_celltype_plot.rds`: lightweight Seurat object for plotting, metadata inspection and cell-type annotation reuse.",
  "- `scripts/read_seurat.R` and `scripts/read_h5ad.py`: example readers for R and Python.",
  "- `provenance/references.tsv`: source dataset references, access links and citation provenance.",
  "",
  "The `Sequence` field is standardized to `scRNA-seq` or `snRNA-seq`; `Technology` records sequencing technology; `Sex` is standardized to `Female`, `Male` or `Not reported`. `Sample` is retained as the biological-sample field and matches `Sample_ID`; `Donor_ID`, `Sample_ID` and `Library_ID` separate donor/sample/library levels. Source-label provenance and curation rules are documented in `metadata/sample_schema_audit.tsv`.",
  "",
  "The expression matrix contains raw count values from Seurat counts layers. `percent_ribo` is not exported because it cannot be reliably reconstructed for the current integrated object.",
  "",
  "The lightweight Seurat object in `objects/` is provided for convenient metadata, clustering, embedding and cell-type visualization workflows. It is not the full raw-count expression object; users requiring complete expression values should start from the 10X-compatible files in `expression/`."
)
writeLines(readme, file.path(out_dir, "README.md"))

manifest_files <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
manifest_files <- manifest_files[file.info(manifest_files)$isdir == FALSE]
manifest_files <- manifest_files[!basename(manifest_files) %in% c("file_manifest.tsv", "md5sum.txt")]
rel <- sub(paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", normalizePath(out_dir, mustWork = TRUE)), "/?"), "", normalizePath(manifest_files))
manifest <- data.frame(
  path = rel,
  file_name = basename(rel),
  size_bytes = as.numeric(file.info(manifest_files)$size),
  md5 = unname(tools::md5sum(manifest_files)),
  stringsAsFactors = FALSE
)
write_tsv(manifest[order(manifest$path), ], file.path(out_dir, "file_manifest.tsv"))
writeLines(paste(manifest$md5, manifest$path), file.path(out_dir, "md5sum.txt"))
message("ScienceDB sidecar files written to: ", out_dir)
