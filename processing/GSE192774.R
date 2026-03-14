source("functions/prepare_env.R")

log_message("Starting GSE192774 integration (Human and Chimpanzee only)...")

data_dir <- "../../data/BrainOmicsData/raw"
res_dir <- check_dir("../../data/BrainOmicsData/processed/GSE192774/")

gse192772_dir <- file.path(data_dir, "GSE192772")
gse192772_files <- list(
  "Astrocyte_ATAC_Chimp" = "GSE192772_Seurat_Astrocyte_ATAC_Chimp.RDS",
  "Astrocyte_ATAC_Human" = "GSE192772_Seurat_Astrocyte_ATAC_Human.RDS",
  "Excitatory_ATAC_Chimp" = "GSE192772_Seurat_Excitatory_ATAC_Chimp.RDS",
  "Excitatory_ATAC_Human" = "GSE192772_Seurat_Excitatory_ATAC_Human.RDS",
  "Inhibitory_ATAC_Chimp" = "GSE192772_Seurat_Inhibitory_ATAC_Chimp.RDS",
  "Inhibitory_ATAC_Human" = "GSE192772_Seurat_Inhibitory_ATAC_Human.RDS",
  "Microglia_ATAC_Chimp" = "GSE192772_Seurat_Microglia_ATAC_Chimp.RDS",
  "Microglia_ATAC_Human" = "GSE192772_Seurat_Microglia_ATAC_Human.RDS",
  "OPCOli_ATAC_Chimp" = "GSE192772_Seurat_OPCOli_ATAC_Chimp.RDS",
  "OPCOli_ATAC_Human" = "GSE192772_Seurat_OPCOli_ATAC_Human.RDS"
)

objects_list_atac <- list()
log_message("Loading GSE192772 (snATACseq) objects...")
for (name in names(gse192772_files)) {
  file_path <- file.path(gse192772_dir, gse192772_files[[name]])

  log_message("  Loading {.val {name}}...")
  obj <- readRDS(file_path)

  if (grepl("Chimp", name)) {
    species <- "Chimpanzee"
  } else if (grepl("Human", name)) {
    species <- "Human"
  }
  obj$CellType_raw <- obj$CellType
  obj$CellType <- obj$broadAnnot

  obj$Species <- species
  obj <- subset(obj, CellType != "MOL")
  obj$Sequence <- "snATACseq"
  obj$Dataset <- "GSE192772"
  obj$Technology <- "10X Genomics"
  obj$Cells <- rownames(obj@meta.data)
  obj$Sample_ID <- obj$Sample
  obj$Brain_Region <- "Posterior cingulate cortex"
  obj$Region <- "Posterior cingulate cortex"

  objects_list_atac[[name]] <- obj
}

log_message("Merging GSE192772 (snATACseq) objects...")
objects_atac <- merge(
  objects_list_atac[[1]],
  y = objects_list_atac[-1],
  project = "GSE192774"
)
rm(objects_list_atac)
gc()

celltype_mapping <- c(
  "Astrocyte" = "Astrocytes",
  "Excitatory" = "Excitatory neurons",
  "Inhibitory" = "Inhibitory neurons",
  "Microglia" = "Microglia",
  "OPC" = "Oligodendrocyte progenitor cells",
  "OPCOli" = "Oligodendrocyte progenitor cells"
)
objects_atac$CellType_raw <- objects_atac$CellType
objects_atac$CellType <- ifelse(
  objects_atac$CellType %in% names(celltype_mapping),
  celltype_mapping[objects_atac$CellType],
  objects_atac$CellType
)
objects_atac$Sex <- ifelse(objects_atac$sex == "M", "Male", "Female")

sample_info <- data.frame(
  Sample = c("c1", "c2", "c3", "c4", "h1", "h2", "h3", "h4"),
  Species_sample = c(rep("Chimpanzee", 4), rep("Human", 4)),
  Sex = c("Male", "Female", "Male", "Female", "Female", "Male", "Male", "Female"),
  age = c(27, 34, 40, 21, 53, 32, 66, 45),
  humanized_age = c(53, 68, 80, 41, 53, 32, 66, 45),
  stringsAsFactors = FALSE
)

metadata_atac <- objects_atac@meta.data
metadata_atac$Age <- NA_character_
metadata_atac$Age_raw <- NA_character_
metadata_atac$Sex <- NA_character_
for (i in seq_len(nrow(sample_info))) {
  sample_id <- sample_info$Sample[i]
  sample_rows <- metadata_atac$Sample == sample_id
  if (sum(sample_rows) > 0) {
    metadata_atac$Sex[sample_rows] <- sample_info$Sex[i]
    if (sample_info$Species_sample[i] == "Chimpanzee") {
      metadata_atac$Age[sample_rows] <- as.character(sample_info$humanized_age[i])
      metadata_atac$Age_raw[sample_rows] <- as.character(sample_info$age[i])
    } else {
      metadata_atac$Age[sample_rows] <- as.character(sample_info$age[i])
      metadata_atac$Age_raw[sample_rows] <- as.character(sample_info$age[i])
    }
  }
}

column_order <- c(
  "Cells", "Dataset", "Technology", "Sequence", "Sample",
  "Sample_ID", "CellType", "CellType_raw", "Brain_Region",
  "Region", "Age", "Age_raw", "Sex", "Species"
)
existing_atac <- intersect(column_order, colnames(metadata_atac))
metadata_atac <- metadata_atac[, existing_atac]
objects_atac@meta.data <- metadata_atac

output_atac <- file.path(res_dir, "GSE192774_ATAC_processed.rds")
log_message("Saving ATAC object to {.file {output_atac}}...")
saveRDS(objects_atac, output_atac)

gse192773_dir <- file.path(data_dir, "GSE192773")
gse192773_files <- list(
  "Astrocyte_RNA" = "GSE192773_Seurat_Astrocyte_RNA.RDS",
  "Excitatory_RNA" = "GSE192773_Seurat_Excitatory_RNA.RDS",
  "Inhibitory_RNA" = "GSE192773_Seurat_Inhibitory_RNA.RDS",
  "Microglia_RNA" = "GSE192773_Seurat_Microglia_RNA.RDS",
  "OPCOli_RNA" = "GSE192773_Seurat_OPCOli_RNA.RDS"
)

log_message("Loading GSE192773 (snRNAseq) objects...")
objects_list_rna <- list()
for (name in names(gse192773_files)) {
  file_path <- file.path(gse192773_dir, gse192773_files[[name]])

  log_message("  Loading {.val {name}}...")
  obj <- readRDS(file_path)

  keep_cells <- grepl(
    "human|chimp|chimpanzee", obj$Species,
    ignore.case = TRUE
  )
  if (sum(keep_cells) == 0) {
    log_message("    Warning: No Human/Chimp cells found, skipping")
    rm(obj)
    gc()
    next
  }

  obj <- obj[, keep_cells]
  log_message(
    "    Filtered to {.val {ncol(obj)}} / {.val {length(keep_cells)}} cells (Human and Chimpanzee only)"
  )

  obj$Species <- ifelse(
    grepl("human", obj$Species, ignore.case = TRUE),
    "Human",
    ifelse(
      grepl("chimp|chimpanzee", obj$Species, ignore.case = TRUE),
      "Chimpanzee",
      obj$Species
    )
  )

  obj$CellType_raw <- obj$CellType
  obj$CellType <- obj$broadAnnot

  obj <- subset(
    obj,
    CellType != "MOL"
  )
  obj$Sequence <- "snRNAseq"
  obj$Dataset <- "GSE192773"
  obj$Technology <- "10X Genomics"
  obj$Cells <- rownames(obj@meta.data)
  obj$Sample_ID <- obj$Sample
  obj$Brain_Region <- "Posterior cingulate cortex"
  obj$Region <- "Posterior cingulate cortex"

  objects_list_rna[[name]] <- obj

  rm(obj)
  gc()
}

log_message("Merging all Seurat objects...")
objects_rna <- merge(
  objects_list_rna[[1]],
  y = objects_list_rna[-1],
  project = "GSE192774"
)
rm(objects_list_rna)

objects_rna$CellType <- ifelse(
  objects_rna$CellType %in% names(celltype_mapping),
  celltype_mapping[objects_rna$CellType],
  objects_rna$CellType
)

metadata <- objects_rna@meta.data

sample_info <- data.frame(
  Sample = c("c1", "c2", "c3", "c4", "h1", "h2", "h3", "h4"),
  Species_sample = c(rep("Chimpanzee", 4), rep("Human", 4)),
  Sex = c("Male", "Female", "Male", "Female", "Female", "Male", "Male", "Female"),
  age = c(27, 34, 40, 21, 53, 32, 66, 45),
  humanized_age = c(53, 68, 80, 41, 53, 32, 66, 45),
  stringsAsFactors = FALSE
)

metadata$Age <- NA_character_
metadata$Age_raw <- NA_character_
metadata$Sex <- NA_character_

for (i in seq_len(nrow(sample_info))) {
  sample_id <- sample_info$Sample[i]
  sample_rows <- metadata$Sample == sample_id

  if (sum(sample_rows) > 0) {
    metadata$Sex[sample_rows] <- sample_info$Sex[i]

    if (sample_info$Species_sample[i] == "Chimpanzee") {
      metadata$Age[sample_rows] <- as.character(sample_info$humanized_age[i])
      metadata$Age_raw[sample_rows] <- as.character(sample_info$age[i])
    } else {
      metadata$Age[sample_rows] <- as.character(sample_info$age[i])
      metadata$Age_raw[sample_rows] <- as.character(sample_info$age[i])
    }
  }
}

column_order <- c(
  "Cells", "Dataset", "Technology", "Sequence", "Sample",
  "Sample_ID", "CellType", "CellType_raw", "Brain_Region",
  "Region", "Age", "Age_raw", "Sex", "Species"
)
existing_columns <- intersect(column_order, colnames(metadata))
metadata <- metadata[, existing_columns]

objects_rna@meta.data <- metadata

output_file <- file.path(res_dir, "GSE192774_processed.rds")
log_message("Saving integrated object to {.file {output_file}}...")
saveRDS(objects_rna, output_file)
