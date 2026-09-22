brainomics_source_first_authors <- function() {
  c(
    AllenM1 = "Bakken", EGAD00001006049 = "Braun",
    EGAS00001006537 = "Cameron", GSE103723 = "Fan",
    GSE104276 = "Zhong", GSE144136 = "Nagy",
    GSE168408 = "Herring", GSE178175 = "Hardwick",
    GSE186538 = "Franjic", GSE199762 = "Nascimento",
    GSE202210 = "Zhu", GSE204683 = "Zhu", GSE207334 = "Ma",
    GSE212606 = "Sziraki", GSE217511 = "Ramos",
    GSE261983 = "Emani", GSE296073 = "Yu", GSE67835 = "Darmanis",
    GSE81475 = "Onorati", GSE97942 = "Lake", HYPOMAP = "Tadross",
    Li_et_al_2018 = "Li", Ma_et_al_2022 = "Ma",
    Nowakowski_et_al_2017 = "Nowakowski", PRJCA015229 = "Yuan",
    ROSMAP = "Mathys; Xiong", SomaMut = "Jeffries", GSE294786 = "Klavert",
    Velmeshev_2023 = "Velmeshev", Wang_2025 = "Wang"
  )
}

brainomics_source_input_contract <- function() {
  datasets <- c(
    "AllenM1", "EGAS00001006537", "GSE104276", "GSE168408", "GSE186538", "GSE204683",
    "GSE207334", "GSE212606", "GSE217511", "GSE296073", "GSE67835", "GSE81475",
    "GSE97942", "HYPOMAP", "Li_et_al_2018", "Ma_et_al_2022", "PRJCA015229", "ROSMAP",
    "SomaMut", "GSE294786", "Velmeshev_2023", "Wang_2025"
  )
  technology <- c(
    "10x Genomics", "10x Genomics", "modified Smart-seq2 with UMI barcodes",
    "10x Genomics Chromium", "10x Genomics", "10x Genomics Multiome",
    "10x Genomics Multiome", "EasySci-RNA", "10x Genomics", "10x Genomics",
    "Fluidigm C1", "Fluidigm C1", "snDrop-seq", "10x Genomics Chromium 3-prime v3.1",
    "Fluidigm C1; 10x Genomics", "10x Genomics",
    "10x Genomics 3-prime Gene Expression; 10x Genomics Multiome", "10x Genomics",
    "10x Genomics", "10x Genomics", "10x Genomics", "10x Genomics Multiome"
  )
  modality <- c(
    "snRNA-seq", "snRNA-seq", "scRNA-seq", "snRNA-seq", "snRNA-seq",
    "snRNA-seq (RNA arm of multiome)", "snRNA-seq (RNA arm of multiome)", "snRNA-seq",
    "snRNA-seq", "snRNA-seq", "scRNA-seq", "scRNA-seq", "snRNA-seq", "snRNA-seq",
    "scRNA-seq; snRNA-seq", "snRNA-seq",
    "snRNA-seq (standalone RNA plus RNA arm of multiome)", "snRNA-seq", "snRNA-seq",
    "snRNA-seq", "snRNA-seq", "snRNA-seq (RNA arm of multiome)"
  )
  sequencing_platform <- c(
    "Illumina NovaSeq 6000", "Illumina NovaSeq 6000", "Illumina HiSeq 4000",
    "Illumina NextSeq 550; Illumina NovaSeq 6000", "Illumina HiSeq 4000",
    "Illumina NovaSeq 6000", "Illumina NovaSeq 6000", "Illumina NovaSeq 6000",
    "Illumina NovaSeq 6000", "Illumina NovaSeq 6000",
    "Illumina MiSeq; Illumina NextSeq 500", "Illumina HiSeq 2000",
    "Illumina HiSeq 2500", "Illumina NovaSeq 6000",
    "fetal: Illumina HiSeq 2000; adult: Not reported in retained source files",
    "Illumina NovaSeq 6000", "Illumina NovaSeq 6000",
    "NextSeq 500/550 or NovaSeq 6000; library-level platform mapping unavailable",
    "Illumina NovaSeq 6000", "Illumina HiSeq 2500", "Illumina NovaSeq 6000",
    "Illumina NovaSeq 6000"
  )
  library_chemistry <- c(
    "10x Chromium 3-prime v3", "10x Chromium 3-prime v3",
    "modified Smart-seq2 with UMI barcodes",
    "10x Chromium 3-prime v2; 10x Chromium 3-prime v3", "10x Chromium 3-prime v3",
    "10x Chromium Next GEM Single Cell Multiome ATAC plus Gene Expression",
    "10x Chromium Single Cell Multiome ATAC plus Gene Expression",
    "EasySci-RNA combinatorial indexing workflow",
    "10x Chromium 3-prime Gene Expression v3",
    "10x Chromium Next GEM Single Cell 3-prime v3.1",
    "Fluidigm C1 + Clontech SMARTer/SMART-seq + Nextera XT",
    "Fluidigm C1 capture; SMARTer Ultra Low RNA cDNA; Illumina Nextera XT library",
    "snDrop-seq with Nextera XT", "10x Genomics Chromium Single Cell 3-prime v3.1",
    "fetal: Fluidigm C1 + SMARTer Ultra Low RNA + Nextera XT; adult: 10x Genomics 3-prime gene expression, version not reported",
    "10x Chromium 3-prime v3",
    "10x Chromium 3-prime Gene Expression; 10x Chromium Single Cell Multiome ATAC plus Gene Expression",
    "10x Chromium Single Cell 3-prime v3",
    "10x Next GEM Single Cell 3-prime v3 or v3.1",
    "10x Chromium v3 Single Cell RNA Gene Expression", "10x Chromium 3-prime v2",
    "10x Chromium Next GEM Single Cell Multiome ATAC plus Gene Expression"
  )
  original_file_type <- c(
    "author CSV gene-by-nucleus count matrix",
    "five gzipped tab-delimited gene-by-nucleus raw-count matrices",
    "gzipped author XLS gene-by-cell UMI-count matrix, cross-audited against 39 GEO per-GSM tables",
    "author AnnData H5AD with the full-count matrix in X",
    "gzipped Matrix Market gene-by-nucleus matrix plus gene table",
    "author RDS sparse gene-by-nucleus count matrix",
    "gzipped Matrix Market RNA-count matrix plus gene table (Multiome RNA arm)",
    "gzipped tab-delimited EasySci-RNA sparse gene-count matrix plus gene annotation",
    "33 GEO 10x feature-barcode Matrix Market triplets",
    "author Seurat RDS counts layer", "466 gzipped per-cell CSV gene-count tables",
    "tab-delimited gene-by-cell count matrix",
    "three tab-delimited gene-by-nucleus snDrop-seq UMI-count matrices",
    "CELLxGENE AnnData H5AD counts layer",
    "two PsychENCODE RData objects containing fetal count2 and adult umi.raw matrices",
    "author Seurat RDS counts assay (Ma_Sestan_mat.rds)",
    "six 10x filtered feature-barcode matrices from OMIX archives",
    "10x Matrix Market gene-count matrix",
    "public processed Seurat RDS with raw counts",
    "cross-species tab-delimited count matrix converted losslessly to H5AD",
    "CELLxGENE AnnData H5AD raw/count matrix plus UCSC metadata",
    "CELLxGENE AnnData H5AD RNA count matrix (Multiome RNA arm)"
  )
  expression_units <- c(
    "raw integer UMI counts including intronic reads", "raw integer UMI counts",
    "raw integer UMI transcript counts",
    "raw integer UMI counts (full counts; not ds_norm_cts)",
    "raw integer UMI counts (exonic plus intronic nuclear reads)",
    "raw integer RNA UMI counts", "raw integer RNA UMI counts",
    "raw integer molecule counts", "raw integer UMI counts from pre-mRNA",
    "raw integer UMI counts including introns",
    "raw integer read counts (HTSeq; not UMI)", "raw integer read counts",
    "raw integer UMI counts", "raw integer UMI counts including introns",
    "raw integer counts (fetal read counts; adult UMI counts)",
    "raw integer UMI counts", "raw integer RNA UMI counts",
    "raw integer UMI counts including pre-mRNA",
    "raw integer UMI counts including introns", "raw integer UMI counts",
    "raw integer UMI counts", "raw integer RNA UMI counts"
  )
  genome_build <- c(
    "GRCh38 / NCBI RefSeq GRCh38.p2", "hg38", "hg19", "hg19", "hg38",
    "hg38 / GRCh38 (Cell Ranger ARC 1.0.0)",
    "study-specific reciprocal-exon-liftOver consensus genome", "hg38",
    "GRCh38 pre-mRNA reference", "10x GRCh38-2020-A", "hg19", "hg38", "GENCODE GRCh38",
    "GRCh38 / Ensembl 98",
    "fetal arm: hg38; adult arm: Not reported in retained source files",
    "study-specific reciprocal-exon-liftOver consensus genome", "GRCh38", "GRCh38",
    "GRCh38 / GENCODE v32", "humanized consensus genome in hg38 coordinates", "GRCh38",
    "GRCh38; GENCODE v32 / Ensembl 98"
  )
  result <- data.frame(
    Dataset = datasets,
    Technology = technology,
    Modality = modality,
    Sequencing_Platform = sequencing_platform,
    Library_Chemistry = library_chemistry,
    Original_File_Type = original_file_type,
    Expression_Units = expression_units,
    Genome_Transcriptome_Build = genome_build,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (nrow(result) != 22L || anyDuplicated(result$Dataset) ||
    anyNA(result) || any(!nzchar(as.matrix(result)))) {
    stop("ScienceDB source-input contract is incomplete")
  }
  result
}

load_sciencedb_publication_summary <- function(repo_dir, datasets) {
  if (all(grepl("^Fixture_", datasets))) {
    return(NULL)
  }
  file <- file.path(
    repo_dir, "results", "integration_25", "publication_tables",
    "source_dataset_summary.tsv"
  )
  if (!file.exists(file)) {
    stop("Missing formal source dataset summary: ", file)
  }
  result <- utils::read.delim(
    file,
    sep = "\t", quote = "", stringsAsFactors = FALSE,
    check.names = FALSE, na.strings = c("", "NA")
  )
  required <- c(
    "Dataset", "Technologies", "Modalities", "Sequencing_Platforms",
    "Library_Chemistries"
  )
  if (any(!required %in% names(result)) || anyDuplicated(result$Dataset) ||
    !setequal(as.character(result$Dataset), datasets)) {
    stop("Formal source dataset summary does not exactly cover the release cohort")
  }
  result[match(datasets, result$Dataset), required, drop = FALSE]
}

metadata_value <- function(meta, candidates, default = NA_character_) {
  out <- rep(default, nrow(meta))
  for (candidate in candidates) {
    if (!candidate %in% names(meta)) next
    values <- as.character(meta[[candidate]])
    out_missing <- is.na(out) | !nzchar(out) |
      out %in% c("NA", "Not reported")
    value_available <- !is.na(values) & nzchar(values) &
      !values %in% c("NA", "Not reported")
    out[out_missing & value_available] <- values[out_missing & value_available]
  }
  out
}

format_age_number <- function(value) {
  digits <- ifelse(abs(value) < 1, 3L, 2L)
  formatted <- mapply(
    function(current, current_digits) {
      formatC(
        round(current, current_digits),
        format = "f", digits = current_digits
      )
    },
    value, digits,
    USE.NAMES = FALSE
  )
  sub("[.]?0+$", "", formatted)
}

standardized_age_label <- function(value, unit) {
  missing <- is.na(value) | is.na(unit) | !nzchar(as.character(unit))
  number <- format_age_number(value)
  result <- ifelse(
    tolower(as.character(unit)) == "pcw",
    paste(number, "PCW"),
    number
  )
  result[missing] <- NA_character_
  result
}

build_sciencedb_cell_metadata <- function(meta, cells) {
  cell_meta <- meta[cells, , drop = FALSE]
  if (anyNA(rownames(cell_meta)) || !identical(rownames(cell_meta), cells)) {
    stop("ScienceDB metadata could not be aligned to the expression columns")
  }
  age_value_field <- intersect(c("Age_num", "age_value"), names(cell_meta))
  age_unit_field <- intersect(c("Unit", "age_unit"), names(cell_meta))
  age <- if (length(age_value_field) && length(age_unit_field)) {
    standardized_age_label(
      suppressWarnings(as.numeric(cell_meta[[age_value_field[[1L]]]])),
      cell_meta[[age_unit_field[[1L]]]]
    )
  } else {
    metadata_value(cell_meta, c("Age", "age_raw"))
  }
  age_fallback <- metadata_value(cell_meta, c("Age", "age_raw"))
  age_missing <- is.na(age) | !nzchar(age)
  age[age_missing] <- age_fallback[age_missing]

  result <- data.frame(
    Cells = cells,
    Dataset = metadata_value(cell_meta, c("Dataset")),
    Original_Sample_ID = metadata_value(
      cell_meta, c("Original_Sample_ID", "Original_Source_Sample_ID")
    ),
    Donor_ID = metadata_value(cell_meta, c("Donor_ID", "Global_Donor_ID")),
    Sample_ID = metadata_value(cell_meta, c("Sample_ID")),
    Library_ID = metadata_value(cell_meta, c("Library_ID")),
    Technology = metadata_value(
      cell_meta, c("sequencing_technology", "Technology")
    ),
    Modality = metadata_value(
      cell_meta, c(
        "sequencing_modality_standardized", "Sequence", "Modality"
      )
    ),
    Age = age,
    AgeIntervalID = metadata_value(cell_meta, c("AgeIntervalID")),
    Sex = metadata_value(cell_meta, c("sex_standardized", "Sex")),
    BrainRegion = metadata_value(cell_meta, c("BrainRegion")),
    Source_CellType = metadata_value(
      cell_meta, c("Source_CellType", "source_cell_type_label", "CellType_raw")
    ),
    Cluster = metadata_value(
      cell_meta, c("Cluster", "Resolution_Cluster", "seurat_clusters")
    ),
    CellType = metadata_value(cell_meta, c("CellType")),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!identical(result$Cells, cells) || ncol(result) != 15L ||
    anyDuplicated(result$Cells)) {
    stop("ScienceDB 15-column metadata contract failed")
  }
  result
}

collapse_reported <- function(values) {
  values <- trimws(as.character(values))
  values <- sort(unique(values[!is.na(values) & nzchar(values)]))
  if (length(values) == 0L) "Not reported" else paste(values, collapse = ";")
}

known_unique <- function(values) {
  values <- trimws(as.character(values))
  length(unique(values[!is.na(values) & nzchar(values)]))
}

build_sciencedb_dataset_summary <- function(meta, metadata, repo_dir) {
  source_file <- file.path(repo_dir, "data", "source_access_summary.tsv")
  if (!file.exists(source_file)) stop("Missing source access summary: ", source_file)
  source <- read.delim(
    source_file,
    sep = "\t", quote = "", stringsAsFactors = FALSE,
    check.names = FALSE, na.strings = c("", "NA")
  )
  authors <- brainomics_source_first_authors()
  datasets <- sort(unique(metadata$Dataset))
  publication_summary <- load_sciencedb_publication_summary(repo_dir, datasets)
  input_contract <- brainomics_source_input_contract()
  real_release <- !all(grepl("^Fixture_", datasets))
  if (real_release && !setequal(input_contract$Dataset, datasets)) {
    stop("ScienceDB source-input contract does not match the release cohort")
  }
  rows <- lapply(datasets, function(dataset) {
    selected <- metadata$Dataset == dataset
    source_row <- source[source$dataset == dataset, , drop = FALSE]
    if (nrow(source_row) == 0L && grepl("^Fixture_", dataset)) {
      source_row <- data.frame(
        source_accession = NA_character_, publication_doi = NA_character_,
        verified_title = "Synthetic clean-room fixture",
        access_status = "synthetic test data", stringsAsFactors = FALSE
      )
    }
    if (nrow(source_row) != 1L) {
      stop("Source access summary does not uniquely cover ", dataset)
    }
    source_cells <- rownames(meta)[as.character(meta$Dataset) == dataset]
    source_meta <- meta[source_cells, , drop = FALSE]
    first_author <- unname(authors[dataset])
    if (length(first_author) != 1L || is.na(first_author)) {
      first_author <- "Not reported"
    }
    publication_row <- if (is.null(publication_summary)) {
      NULL
    } else {
      publication_summary[publication_summary$Dataset == dataset, , drop = FALSE]
    }
    contract_row <- if (!real_release) {
      NULL
    } else {
      input_contract[input_contract$Dataset == dataset, , drop = FALSE]
    }
    reported_or_metadata <- function(publication_field, metadata_fields) {
      if (!is.null(publication_row)) {
        value <- collapse_reported(publication_row[[publication_field]])
        if (!identical(value, "Not reported")) {
          return(value)
        }
      }
      collapse_reported(metadata_value(source_meta, metadata_fields))
    }
    platform <- if (is.null(contract_row)) {
      reported_or_metadata(
        "Sequencing_Platforms", c("Sequencing_Platform")
      )
    } else {
      contract_row$Sequencing_Platform
    }
    chemistry <- if (is.null(contract_row)) {
      reported_or_metadata(
        "Library_Chemistries", c("Library_Chemistry")
      )
    } else {
      contract_row$Library_Chemistry
    }
    data.frame(
      Dataset = dataset,
      First_Author = first_author,
      Source_Accession = collapse_reported(source_row$source_accession),
      Publication_DOI = collapse_reported(source_row$publication_doi),
      Publication_Title = collapse_reported(source_row$verified_title),
      Technology = if (!is.null(contract_row)) {
        contract_row$Technology
      } else if (is.null(publication_row)) {
        collapse_reported(metadata$Technology[selected])
      } else {
        collapse_reported(publication_row$Technologies)
      },
      Modality = if (!is.null(contract_row)) {
        contract_row$Modality
      } else if (is.null(publication_row)) {
        collapse_reported(metadata$Modality[selected])
      } else {
        collapse_reported(publication_row$Modalities)
      },
      Sequencing_Platform = platform,
      Library_Chemistry = chemistry,
      Original_File_Type = if (is.null(contract_row)) {
        collapse_reported(metadata_value(source_meta, c("Original_File_Type")))
      } else {
        contract_row$Original_File_Type
      },
      Expression_Units = if (is.null(contract_row)) {
        collapse_reported(metadata_value(source_meta, c("Expression_Units")))
      } else {
        contract_row$Expression_Units
      },
      Genome_Transcriptome_Build = if (is.null(contract_row)) {
        collapse_reported(metadata_value(
          source_meta, c("Genome_Transcriptome_Build")
        ))
      } else {
        contract_row$Genome_Transcriptome_Build
      },
      Cells_or_Nuclei = sum(selected),
      Known_Donors = known_unique(metadata$Donor_ID[selected]),
      Biological_Samples = known_unique(metadata$Sample_ID[selected]),
      Verified_Sequence_Libraries = known_unique(metadata$Library_ID[selected]),
      Access_Status = collapse_reported(source_row$access_status),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  result <- do.call(rbind, rows)
  if (!identical(result$Dataset, datasets) ||
    sum(result$Cells_or_Nuclei) != nrow(metadata)) {
    stop("ScienceDB dataset summary failed dataset or cell coverage")
  }
  result
}

build_sciencedb_references <- function(dataset_summary) {
  dataset_summary[, c(
    "Dataset", "First_Author", "Source_Accession", "Publication_DOI",
    "Publication_Title", "Access_Status"
  ), drop = FALSE]
}

brainomics_verification_status_dictionary <- function() {
  data.frame(
    Verification_Status = c(
      "verified", "source-reported", "derived", "curated",
      "unresolved", "excluded", "not-applicable"
    ),
    Match_Rule = c(
      "status begins with 'verified:' or an exact-source audit states verified",
      "value is explicitly attributed to the source publication or repository",
      "value is produced by a documented deterministic conversion from source fields",
      "value is assigned by a documented manual crosswalk with retained evidence",
      "status states unavailable, unknown, ambiguous, unresolved or not independently identifiable",
      "record has Analysis_Include=FALSE with an explicit exclusion reason",
      "field is not biologically or technically applicable to the record"
    ),
    Definition = c(
      "The identifier or value was matched to an original author/repository field or independently checked against a retained source artifact.",
      "The value is copied from the source study but was not independently reconstructed or cross-validated.",
      "The value follows reproducibly from retained source evidence and the recorded conversion rule; it is not a direct source field.",
      "The standardized value reflects expert curation supported by a retained source-to-standard crosswalk and rationale.",
      "Available evidence is insufficient for a known assignment; the record remains explicit rather than being imputed.",
      "The record remains auditable but is not part of the stated analysis cohort.",
      "No value is expected for this record under the released schema."
    ),
    Known_Unit_Rule = c(
      "May contribute to a known biological or technical-unit total when the field semantics also establish the unit.",
      "Contributes only when the source semantics establish the biological or technical unit; reporting alone is not a license or independent verification.",
      "May contribute only when the deterministic rule preserves the relevant unit identity.",
      "May contribute only when the curated mapping preserves the relevant unit identity and is not an unknown placeholder.",
      "Never contributes to known donor, specimen or verified-library totals.",
      "Never contributes to the analysis-cohort total.",
      "Never contributes to a known-unit total for that field."
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

build_sciencedb_attrition_summary <- function(
  repo_dir, dataset_summary, integration_features = 3000L
) {
  datasets <- as.character(dataset_summary$Dataset)
  if (all(grepl("^Fixture_", datasets))) {
    return(data.frame(
      Dataset = datasets,
      Source_Cells = dataset_summary$Cells_or_Nuclei,
      Diagnosis_Excluded_Cells = 0,
      Duplicate_Input_or_Donor_Excluded_Cells = 0,
      Species_or_Tissue_Excluded_Cells = 0,
      Source_Quality_Excluded_Cells = 0,
      Doublet_Excluded_Cells = 0,
      Metadata_Excluded_Cells = 0,
      Region_Excluded_Cells = 0,
      Other_Excluded_Cells = 0,
      Analysis_Cells = dataset_summary$Cells_or_Nuclei,
      Exported_Cells = dataset_summary$Cells_or_Nuclei,
      Source_Features = NA_integer_,
      Integration_Features = NA_integer_,
      Count_Semantics = "synthetic clean-room fixture",
      Feature_Filtering_Description = "not applicable to the fixture",
      Attrition_Boundary = "fixture cells generated for workflow testing",
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }

  ledger_file <- file.path(
    repo_dir, "..", "..", "data", "BrainOmicsData", "integration_25", "evaluation",
    "inclusion_ledger", "datasets_inclusion_ledger.tsv"
  )
  publication_file <- file.path(
    repo_dir, "..", "..", "data", "BrainOmicsData", "integration_25", "publication_tables",
    "source_dataset_summary.tsv"
  )
  if (!file.exists(ledger_file) || !file.exists(publication_file)) {
    stop("Missing formal inclusion ledger or source dataset summary")
  }
  ledger <- utils::read.delim(
    ledger_file,
    sep = "\t", quote = "", stringsAsFactors = FALSE,
    check.names = FALSE
  )
  publication <- utils::read.delim(
    publication_file,
    sep = "\t", quote = "", stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required_ledger <- c(
    "Dataset", "Analysis_Include", "Exclusion_Reason", "Cells_or_Nuclei"
  )
  required_publication <- c(
    "Dataset", "Source_Cells", "Analysis_Cells", "Source_Features"
  )
  if (any(!required_ledger %in% names(ledger)) ||
    any(!required_publication %in% names(publication)) ||
    anyDuplicated(publication$Dataset) ||
    !setequal(datasets, publication$Dataset)) {
    stop("Formal attrition inputs do not match the release cohort")
  }
  publication <- publication[match(datasets, publication$Dataset), , drop = FALSE]
  ledger <- ledger[ledger$Dataset %in% datasets, , drop = FALSE]
  ledger$Analysis_Include <- as.logical(ledger$Analysis_Include)
  ledger$Cells_or_Nuclei <- as.numeric(ledger$Cells_or_Nuclei)

  classify_exclusion <- function(reason) {
    reason <- tolower(as.character(reason))
    result <- rep("other", length(reason))
    result[grepl("major depressive|parkinson|disease source donor", reason)] <-
      "diagnosis"
    result[grepl(
      "overlapping donor|duplicate donor|previously integrated dataset",
      reason
    )] <- "duplicate"
    result[grepl("non-human|head compartment|non-brain|species", reason)] <-
      "species_tissue"
    result[grepl("poor-quality|low quality|quality-excluded", reason)] <-
      "source_quality"
    result[grepl("doublet", reason)] <- "doublet"
    result[grepl("no matching released cell metadata|metadata", reason)] <-
      "metadata"
    result[grepl("region exclusion|excluded brain region", reason)] <- "region"
    result
  }

  rows <- lapply(datasets, function(dataset) {
    selected <- ledger[ledger$Dataset == dataset, , drop = FALSE]
    if (nrow(selected) == 0L || anyNA(selected$Analysis_Include) ||
      anyNA(selected$Cells_or_Nuclei)) {
      stop("Incomplete inclusion ledger for ", dataset)
    }
    excluded <- selected[!selected$Analysis_Include, , drop = FALSE]
    category <- classify_exclusion(excluded$Exclusion_Reason)
    category_count <- function(value) {
      sum(excluded$Cells_or_Nuclei[category == value], na.rm = TRUE)
    }
    source_cells <- sum(selected$Cells_or_Nuclei)
    analysis_cells <- sum(
      selected$Cells_or_Nuclei[selected$Analysis_Include],
      na.rm = TRUE
    )
    published <- publication[publication$Dataset == dataset, , drop = FALSE]
    if (source_cells != published$Source_Cells ||
      analysis_cells != published$Analysis_Cells) {
      stop("Attrition ledger disagrees with the publication summary for ", dataset)
    }
    data.frame(
      Dataset = dataset,
      Source_Cells = source_cells,
      Diagnosis_Excluded_Cells = category_count("diagnosis"),
      Duplicate_Input_or_Donor_Excluded_Cells = category_count("duplicate"),
      Species_or_Tissue_Excluded_Cells = category_count("species_tissue"),
      Source_Quality_Excluded_Cells = category_count("source_quality"),
      Doublet_Excluded_Cells = category_count("doublet"),
      Metadata_Excluded_Cells = category_count("metadata"),
      Region_Excluded_Cells = category_count("region"),
      Other_Excluded_Cells = category_count("other"),
      Analysis_Cells = analysis_cells,
      Exported_Cells = analysis_cells,
      Source_Features = published$Source_Features,
      Integration_Features = as.integer(integration_features),
      Count_Semantics = paste0(
        "Source_Cells is the complete processed source object entering the ",
        "audited project workflow; upstream author filtering is not silently ",
        "reclassified as project QC."
      ),
      Feature_Filtering_Description = paste0(
        "All released count shards retain complete harmonized source-count ",
        "features; a fixed ", integration_features,
        "-gene consensus HVG set was used only for integration."
      ),
      Attrition_Boundary = paste0(
        "Diagnosis, duplicate-input/donor, species/tissue, source-quality, ",
        "doublet, metadata, region and other exclusions are mutually ",
        "exclusive inclusion-ledger rows."
      ),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  result <- do.call(rbind, rows)
  excluded_columns <- grep("_Excluded_Cells$", names(result), value = TRUE)
  if (!identical(result$Dataset, datasets) ||
    any(result$Source_Cells - rowSums(result[, excluded_columns]) !=
      result$Analysis_Cells) ||
    sum(result$Analysis_Cells) != sum(dataset_summary$Cells_or_Nuclei)) {
    stop("ScienceDB attrition summary failed its conservation checks")
  }
  result
}
