source("functions/prepare_env.R")

data_dir <- "../../data/BrainOmicsData/raw/ROSMAP"
res_dir <- check_dir("../../data/BrainOmicsData/processed/ROSMAP/")

if (!file.exists(file.path(data_dir, "individual_metadata_deidentified.tsv"))) {
  download.file(
    "https://personal.broadinstitute.org/cboix/ad427_data/Data/Metadata/individual_metadata_deidentified.tsv",
    file.path(data_dir, "individual_metadata_deidentified.tsv")
  )
}

log_message("Start loading data...")

metadata1 <- read.csv(
  file.path(
    data_dir, "individual_metadata_deidentified.tsv"
  ),
  sep = "\t"
)

metadata1$Age <- sapply(
  metadata1$age_death, function(x) {
    if (x == "90+") {
      return(90)
    } else if (grepl("\\(", x) && grepl("\\]", x)) {
      numbers <- as.numeric(unlist(regmatches(x, gregexpr("\\d+", x))))
      if (length(numbers) == 2) {
        return(mean(numbers))
      }
    }
    return(NA)
  }
)
metadata1 <- metadata1[metadata1$Pathologic_diagnosis_of_AD == "no", ]
metadata1$Sex <- ifelse(metadata1$msex == "0", "Female", "Male")

metadata2 <- read.csv(
  file.path(
    data_dir, "RNA/meta.tsv"
  ),
  sep = "\t"
)
metadata2 <- metadata2[metadata2$ADdiag3types == "nonAD", ]

metadata <- merge(
  metadata2,
  metadata1,
  by.x = "Individual",
  by.y = "subject",
  all.x = TRUE
)

rownames(metadata) <- metadata$cellId
metadata$Cells <- metadata$cellId
metadata$Dataset <- "ROSMAP"
metadata$Technology <- "10X Genomics"
metadata$Sequence <- "snRNA-seq"
metadata$Sample <- metadata$Individual
metadata$Sample_ID <- metadata$Individual
metadata$CellType_raw <- metadata$Celltype
metadata$Brain_Region <- "Prefrontal cortex"
metadata$Region <- "Prefrontal cortex"

column_order <- c(
  "Cells", "Dataset", "Technology", "Sequence", "Sample",
  "Sample_ID", "CellType_raw", "Brain_Region", "Region", "Age", "Sex"
)
metadata <- metadata[, column_order]
metadata <- na.omit(metadata)

counts <- Read10X(
  file.path(
    data_dir, "RNA"
  ),
  gene.column = 1
)

intersect_cells <- intersect(colnames(counts), rownames(metadata))
counts <- counts[, intersect_cells]

metadata <- metadata[intersect_cells, ]

object <- CreateSeuratObject(
  counts = counts,
  meta.data = metadata
)

log_message("Save data...")
saveRDS(
  object,
  file.path(res_dir, "ROSMAP_processed.rds")
)
