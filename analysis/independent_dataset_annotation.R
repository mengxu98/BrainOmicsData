#!/usr/bin/env Rscript
# Dataset-specific, expression-only annotation. Never import atlas/source labels.
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
})
source("functions/data_paths.R")
args <- commandArgs(trailingOnly = TRUE)
dataset <- args[1]
stopifnot(dataset %in% c("GSE104276", "GSE178175", "GSE202210"))
out <- brainomics_data_path("independent_annotation", dataset)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(42)
data.table::setDTthreads(2)
future::plan("sequential")

if (length(args) > 1L && args[2] == "finalize") {
  cells <- fread(file.path(out, "cells.tsv.gz"), colClasses = c(Cluster = "character"))
  decisions <- fread("annotation/independent_dataset_cluster_labels.tsv", colClasses = c(Cluster = "character"))
  decisions <- decisions[Dataset == dataset]
  stopifnot(!anyDuplicated(decisions$Cluster),
            setequal(cells$Cluster, decisions$Cluster))
  result <- merge(cells, decisions, by = "Cluster", sort = FALSE)
  result <- result[match(cells$Cells, result$Cells)]
  stopifnot(identical(cells$Cells, result$Cells), !anyNA(result$Independent_CellType))
  result[, Cluster_CellType := Independent_CellType]
  result[, QC_Flag := fifelse(nFeature < 200, "Fewer than 200 detected genes",
                             fifelse(percent_mt > 20, "Mitochondrial fraction above 20%", "None of these flags"))]
  result[nFeature < 200, `:=`(Independent_CellType = "Unresolved (low coverage)",
                             Evidence_Status = "Unresolved",
                             Evidence = "Insufficient detected genes for reliable individual assignment; cell retained")]
  result[, Annotation_Origin := "BrainOmicsData independent within-dataset reannotation"]
  result[, Annotation_Method := "Expression-only PCA/SNN clusters; reviewed differential and canonical markers; seed 42"]
  fwrite(result, file.path(out, "independent_labels.tsv.gz"), sep = "\t")
  saveRDS(result, file.path(out, "independent_labels.rds"))
  fwrite(result[, .N, by = .(Independent_CellType, Evidence_Status)],
         file.path(out, "label_counts.tsv"), sep = "\t")
  quit(save = "no")
}

refine <- length(args) > 1L && args[2] == "refine"
evidence_only <- length(args) > 1L && args[2] == "evidence"
if (!refine && !evidence_only && file.exists(file.path(out, "expression_clusters.rds"))) {
  stop("Existing independent analysis: archive it before rerunning")
}
input <- brainomics_data_path("processed", dataset, paste0(dataset, "_processed.rds"))
if (refine || evidence_only) {
  parent <- args[3]
  full <- readRDS(file.path(out, "expression_clusters.rds"))
  input_cells <- colnames(full)
  sample_id <- full$Sample
  object <- if (refine) subset(full, idents = parent) else full
} else {
  message(dataset, ": reading expression")
  original <- readRDS(input)
  input_cells <- colnames(original)
  # Read only sample identity for diagnostics; labels/embeddings never enter analysis.
  sample_id <- as.character(original[[]]$Sample)
  if (length(Layers(original, assay = "RNA", search = "^counts")) > 1L) {
    original <- JoinLayers(original, assay = "RNA")
  }
  counts <- GetAssayData(original, assay = "RNA", layer = "counts")
  stopifnot(identical(colnames(counts), input_cells), !anyDuplicated(input_cells),
            all(Matrix::colSums(counts) > 0), all(counts@x >= 0))
  object <- CreateSeuratObject(counts, min.cells = 0, min.features = 0)
  rm(original, counts)
  gc()
  object$Sample <- sample_id
  object$percent_mt <- PercentageFeatureSet(object, pattern = "^MT-")
}
if (!evidence_only) {
  message(dataset, ": normalize, select variable genes, PCA and clustering")
  object <- NormalizeData(object, verbose = FALSE)
  object <- FindVariableFeatures(object, nfeatures = 3000, verbose = FALSE)
  # Prevent sex/mitochondrial/ribosomal abundance from defining the PCA basis.
  features <- setdiff(VariableFeatures(object),
                     grep("^MT-|^RP[SL][0-9]|^XIST$|^RPS4Y1$|^DDX3Y$|^KDM5D$|^UTY$",
                          rownames(object), value = TRUE))
  object <- ScaleData(object, features = features, verbose = FALSE)
  object <- RunPCA(object, features = features, npcs = 30, seed.use = 42, verbose = FALSE)
  object <- FindNeighbors(object, dims = 1:30, verbose = FALSE)
  object <- FindClusters(object, resolution = 0.6, random.seed = 42, verbose = FALSE)
}
if (refine) {
  identities <- setNames(as.character(Idents(full)), colnames(full))
  identities[colnames(object)] <- paste0(parent, ".", as.character(Idents(object)))
  Idents(full) <- identities
  object <- full
  rm(full)
}
stopifnot(identical(colnames(object), input_cells))
cells <- data.table(Cells = paste0(dataset, "::", input_cells),
                    Cluster = as.character(Idents(object)), Sample = sample_id,
                    nCount = object$nCount_RNA, nFeature = object$nFeature_RNA,
                    percent_mt = object$percent_mt)
if (any(grepl(paste0("^", dataset, "::"), input_cells))) {
  cells[, Cells := input_cells]
}
fwrite(cells, file.path(out, "cells.tsv.gz"), sep = "\t")
fwrite(cells[, .(N = .N, median_genes = median(nFeature), median_mt = median(percent_mt)),
             by = .(Cluster, Sample)], file.path(out, "cluster_sample_qc.tsv"), sep = "\t")
if (!evidence_only) saveRDS(object, file.path(out, "expression_clusters.rds"), compress = FALSE)

# Canonical markers span adult cortex and fetal neural/glial lineages.
# Sources: Zhong et al. DOI 10.1038/nature25980; Hardwick et al.
# DOI 10.1038/s41587-022-01231-3; Zhu et al. DOI 10.1126/scitranslmed.abo1997.
panel <- c("SLC17A7", "SLC17A6", "SATB2", "CUX1", "CUX2", "TBR1", "BCL11B",
           "GAD1", "GAD2", "SLC32A1", "DLX1", "DLX2", "DLX5", "DLX6", "LHX6", "SST", "PVALB", "VIP", "RELN",
           "RBFOX3", "SNAP25", "SYT1", "STMN2", "DCX", "NEUROD1", "NEUROD2",
           "SOX2", "PAX6", "HES1", "HES5", "VIM", "NES", "HOPX", "EOMES", "NEUROG2", "MKI67", "TOP2A",
           "AQP4", "ALDH1L1", "SLC1A2", "SLC1A3", "GFAP", "FGFR3", "GJA1",
           "PDGFRA", "CSPG4", "VCAN", "OLIG1", "OLIG2", "SOX10", "GPR17", "BCAS1", "ENPP6",
           "MBP", "PLP1", "MOG", "MOBP", "MAG", "CNP",
           "P2RY12", "CX3CR1", "CSF1R", "AIF1", "TYROBP", "C1QA", "C1QB", "CD74", "CD163", "MRC1",
           "CLDN5", "PECAM1", "VWF", "FLT1", "KDR", "RGS5", "PDGFRB", "ACTA2", "TAGLN", "COL1A1", "DCN", "LUM",
           "FOXJ1", "PIFO", "TTR", "KRT8", "KRT18", "CD3D", "CD3E", "TRAC", "PTPRC", "NKG7", "MS4A1", "HBB")
data <- GetAssayData(object, layer = "data")
panel <- intersect(panel, rownames(data))
evidence <- rbindlist(lapply(levels(Idents(object)), function(cl) {
  ii <- which(as.character(Idents(object)) == cl)
  data.table(Cluster = cl, Gene = panel,
             Mean_LogExpression = Matrix::rowMeans(data[panel, ii, drop = FALSE]),
             Pct_Expressing = 100 * Matrix::rowMeans(data[panel, ii, drop = FALSE] > 0))
}))
fwrite(evidence, file.path(out, "canonical_marker_evidence.tsv"), sep = "\t")
message(dataset, ": finding differential markers")
markers <- as.data.table(presto::wilcoxauc(data, as.character(Idents(object))))
# Compute fold changes on linear normalized expression, not differences of logs.
membership <- sparseMatrix(i = seq_len(ncol(data)), j = as.integer(Idents(object)),
                          x = 1, dims = c(ncol(data), nlevels(Idents(object))))
sizes <- as.numeric(table(Idents(object)))
linear <- data
linear@x <- expm1(linear@x)
sums <- as.matrix(linear %*% membership)
inside <- sweep(sums, 2, sizes, "/")
outside <- sweep(Matrix::rowSums(linear) - sums, 2, ncol(data) - sizes, "/")
fc <- log2((inside + 1) / (outside + 1))
markers[, avg_log2FC := fc[cbind(match(feature, rownames(data)),
                                match(group, levels(Idents(object))))]]
markers <- markers[pct_in >= 10 & avg_log2FC >= 0.25]
setnames(markers, c("feature", "group", "pval", "padj", "pct_in", "pct_out"),
         c("gene", "cluster", "p_val", "p_val_adj", "pct.1", "pct.2"))
markers[, `:=`(pct.1 = pct.1 / 100, pct.2 = pct.2 / 100)]
fwrite(markers, file.path(out, "cluster_markers.tsv.gz"), sep = "\t")
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
refined_parents <- unique(sub("\\..*$", "", grep("\\.", levels(Idents(object)), value = TRUE)))
writeLines(c(paste("Input:", input), "Input fields used: RNA counts, cell IDs, Sample only",
             "No current/source/external labels, atlas embeddings, or label transfer used.",
             "All input cells retained; independent PCA/SNN per dataset; seed=42; resolution=0.6.",
             if (length(refined_parents)) paste("Expression-based subclustering of ambiguous parent cluster:", paste(refined_parents, collapse = ",")) else "Initial dataset clustering",
             "Marker statistics are descriptive annotation evidence, not donor-level inference."),
           file.path(out, "method.txt"))
cat("Presto batched Wilcoxon; BH-adjusted p values; log2FC from mean linear normalized expression with pseudocount 1.\n",
    file = file.path(out, "method.txt"), append = TRUE)
message(dataset, ": evidence complete; cluster labels await marker review")
