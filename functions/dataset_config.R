dataset_config <- function(dataset) {
  configs <- list(
    Catching_2026 = list(
      input_root = "raw",
      input_file = "final_rna_data.h5ad",
      expected_input_bytes = 18439154935,
      matrix_group = "X",
      expected_cells = 1501089L,
      expected_features = 38606L,
      expected_analysis_cells = 1493496L,
      expected_analysis_donors = 355L,
      matrix_storage = "dgCMatrix_shards",
      layer_by = "Original_Technical_Batch_ID",
      expected_layers = 10L,
      refresh_shard_metadata_on_metadata_only = TRUE,
      metadata_only_protected_changes = list(
        Library_ID = list(
          expected_previous_unique = 10L,
          expected_new_unique = 357L,
          previous_source_column = "Original_Library_ID",
          expected_previous_source_unique = 10L,
          new_source_column = "Original_Library_ID",
          expected_new_source_unique = 357L,
          evidence = paste(
            "replace the legacy pooled sequencing-batch Library_ID with",
            "the source SampleID-specific gene-expression library; retain",
            "the ten source batch values only as Technical_Batch_ID"
          )
        )
      ),
      default_analysis_role = "holdout",
      external_metadata = list(
        neuron_subtype = list(
          root = "raw",
          file = "neuron_subtype.csv",
          sep = ",",
          key_column = "Catching_Cell_ID",
          match_column = "Catching_Neuron_Metadata_Matched",
          require_all_matrix_cells = TRUE,
          require_all_metadata_cells = TRUE,
          transform = function(meta) {
            required <- c(
              "cell_barcode", "SampleID", "cell_type",
              "InN_subtype", "ExN_subtype"
            )
            missing <- setdiff(required, names(meta))
            if (length(missing) > 0L) {
              stop(
                "Catching neuron-subtype metadata is missing columns: ",
                paste(missing, collapse = ", ")
              )
            }
            data.frame(
              Catching_Cell_ID = paste(
                meta$cell_barcode,
                meta$SampleID,
                sep = "_"
              ),
              Catching_Source_Sample_ID = as.character(meta$SampleID),
              Catching_Source_Cell_Type = as.character(meta$cell_type),
              Catching_InN_Subtype = as.character(meta$InN_subtype),
              Catching_ExN_Subtype = as.character(meta$ExN_subtype),
              stringsAsFactors = FALSE,
              check.names = FALSE
            )
          }
        ),
        cohort_demographics = list(
          root = "raw",
          file = "TableS1_CohortDemographics.xlsx",
          key_column = "Catching_Demographics_Sample_ID",
          matrix_key_column = "SampleID",
          match_column = "Catching_Demographics_Matched",
          require_all_matrix_cells = TRUE,
          require_all_metadata_cells = FALSE,
          transform = function(meta) {
            required <- c(
              "SampleID", "Age", "Sex", "PMI", "Ethnicity", "Race",
              "BrainBank"
            )
            missing <- setdiff(required, names(meta))
            if (length(missing) > 0L) {
              stop(
                "Catching Table S1 is missing columns: ",
                paste(missing, collapse = ", ")
              )
            }
            sample_id <- as.character(meta$SampleID)
            numeric_hbcc <- grepl("^[0-9]+(?:[.]0+)?$", sample_id)
            sample_id[numeric_hbcc] <- sub(
              "[.]0+$",
              "",
              sample_id[numeric_hbcc]
            )
            sample_id[numeric_hbcc] <- paste0(
              "HBCC-",
              sample_id[numeric_hbcc]
            )
            data.frame(
              Catching_Demographics_Sample_ID = sample_id,
              Catching_TableS1_Age = as.numeric(meta$Age),
              Catching_TableS1_Sex = as.character(meta$Sex),
              Catching_TableS1_PMI = as.numeric(meta$PMI),
              Catching_TableS1_Ethnicity = as.character(meta$Ethnicity),
              Catching_TableS1_Race = as.character(meta$Race),
              Catching_TableS1_Brain_Bank = as.character(meta$BrainBank),
              stringsAsFactors = FALSE,
              check.names = FALSE
            )
          }
        )
      ),
      raw_metadata_transform = function(meta) {
        required <- c(
          "Original_Cell_ID", "SampleID", "cell_type",
          "Catching_Neuron_Metadata_Matched",
          "Catching_Demographics_Matched",
          "Catching_Source_Sample_ID",
          "Catching_Source_Cell_Type", "Catching_InN_Subtype",
          "Catching_ExN_Subtype", "Catching_TableS1_Age",
          "Catching_TableS1_Sex", "Catching_TableS1_PMI",
          "Catching_TableS1_Ethnicity",
          "Catching_TableS1_Brain_Bank"
        )
        missing <- setdiff(required, names(meta))
        if (length(missing) > 0L) {
          stop(
            "Catching joined metadata is missing columns: ",
            paste(missing, collapse = ", ")
          )
        }
        if (any(!meta$Catching_Neuron_Metadata_Matched)) {
          stop("Catching neuron-subtype metadata did not match every cell")
        }
        if (any(!meta$Catching_Demographics_Matched)) {
          stop("Catching Table S1 metadata did not match every cell")
        }
        meta$External_Metadata_Matched <-
          meta$Catching_Neuron_Metadata_Matched
        if (!identical(
          as.character(meta$SampleID),
          as.character(meta$Catching_Source_Sample_ID)
        )) {
          stop("Catching sample IDs differ after neuron-subtype join")
        }
        if (!identical(
          as.character(meta$cell_type),
          as.character(meta$Catching_Source_Cell_Type)
        )) {
          stop("Catching broad cell types differ after subtype join")
        }
        if (any(as.character(meta$Sex) !=
          as.character(meta$Catching_TableS1_Sex)) ||
          any(abs(as.numeric(meta$PMI) -
            as.numeric(meta$Catching_TableS1_PMI)) > 1e-8) ||
          any(as.character(meta$Ethnicity) !=
            as.character(meta$Catching_TableS1_Ethnicity)) ||
          any(as.character(meta$Brain_bank) !=
            as.character(meta$Catching_TableS1_Brain_Bank))) {
          stop("Catching H5AD donor metadata differs from Table S1")
        }
        h5ad_age <- as.numeric(meta$Age)
        table_age <- as.numeric(meta$Catching_TableS1_Age)
        age_discrepancy <- abs(h5ad_age - table_age) > 1e-8
        if (sum(age_discrepancy) == 0L ||
          length(unique(meta$SampleID[age_discrepancy])) != 16L ||
          any(h5ad_age[age_discrepancy] != 89.9) ||
          any(table_age[age_discrepancy] <= 90 | table_age[age_discrepancy] > 100)) {
          stop("Catching age protection pattern differs from the source audit")
        }
        meta$Catching_H5AD_Age <- h5ad_age
        meta$Age_Source_Raw <- as.character(table_age)
        meta$Age_Source_Unit <- "years"
        meta$Age_Source_Basis <- "postnatal age at death"
        meta$Age_Harmonization_Input <- paste(
          format(table_age, scientific = FALSE, trim = TRUE),
          "years"
        )
        meta$Age_Conversion_Formula <-
          "Table S1 age in years retained without conversion"
        meta$Age_Conversion_Confidence <-
          "source supplementary donor table"
        meta$Age_Conversion_Applied <- FALSE
        meta$Age_Source_Discrepancy <- ifelse(
          age_discrepancy,
          paste(
            "processed H5AD age is 89.9; public Table S1 reports",
            table_age
          ),
          NA_character_
        )
        inn <- normalize_missing_metadata(meta$Catching_InN_Subtype)
        exn <- normalize_missing_metadata(meta$Catching_ExN_Subtype)
        if (any(!is.na(inn) & !is.na(exn))) {
          stop("Catching cell has both inhibitory and excitatory subtypes")
        }
        neuron_subtype <- ifelse(!is.na(inn), inn, exn)
        broad <- normalize_missing_metadata(
          meta$Catching_Source_Cell_Type
        )
        meta$Catching_Source_Neuron_Subtype <- neuron_subtype
        meta$Catching_Source_Fine_Cell_Type <- ifelse(
          is.na(neuron_subtype),
          broad,
          neuron_subtype
        )
        meta
      },
      field_map = list(
        cell_id = "Original_Cell_ID",
        donor_id = "SampleID",
        sample_id = "SampleID",
        specimen_id = "SampleID",
        library_id = "SampleID",
        source_record_id = "SampleID",
        technical_batch_id = "batch",
        species = NULL,
        species_raw = NULL,
        age = "Age_Harmonization_Input",
        sex = "Sex",
        brain_region = NULL,
        brain_region_source = NULL,
        brain_region_ontology_label = NULL,
        brain_region_ontology_id = NULL,
        diagnosis = NULL,
        cell_type_original = "Catching_Source_Fine_Cell_Type",
        cell_type = "Catching_Source_Fine_Cell_Type",
        cell_type_level_1 = "Catching_Source_Cell_Type",
        cell_type_level_2 = "Catching_Source_Neuron_Subtype",
        cell_type_level_3 = NULL,
        cell_type_cluster_id = "leiden_2",
        cell_type_ontology_label = NULL,
        cell_type_ontology_id = NULL,
        technology = NULL,
        modality = NULL,
        assay = NULL,
        sequencing_platform = NULL,
        chemistry = NULL
      ),
      constants = list(
        species = "Homo sapiens",
        brain_region = "Dorsolateral prefrontal cortex",
        brain_region_source = "Dorsolateral prefrontal cortex",
        brain_region_ontology_label = "dorsolateral prefrontal cortex",
        brain_region_ontology_id = "UBERON:0009834",
        brain_region_ontology_source =
          "UBERON 2026-04-01 exact-label lookup",
        brain_region_anatomical_level = "cortical region",
        brain_region_mapping_method =
          "exact source anatomical description mapped to UBERON",
        brain_region_mapping_confidence = "publication methods verified",
        diagnosis = "Neurologically unaffected",
        technology = "10x Genomics Multiome",
        modality = "snRNA-seq",
        donor_id_semantics = "brain-bank donor",
        specimen_id_semantics = "donor-level prefrontal cortex specimen",
        specimen_id_verification_status = "source reported",
        library_id_semantics = paste(
          "donor-level prefrontal-cortex gene-expression library from the",
          "paired 10x Multiome assay"
        ),
        library_id_verification_status =
          "verified: Cell Reports 2026 methods and Zenodo metadata",
        library_id_evidence = paste(
          "DOI 10.1016/j.celrep.2026.117110 describes construction of the",
          "gene-expression library for each SampleID before pooling the",
          "libraries into sequencing batches; Zenodo record 20834804 retains",
          "the exact SampleID field for every RNA nucleus"
        ),
        source_record_id_semantics = "source SampleID record",
        technical_batch_semantics = paste(
          "gene-expression sequencing batch spanning approximately 40",
          "paired-multiome donor samples"
        ),
        technical_batch_verification_status =
          "verified: Cell Reports 2026 methods",
        technical_batch_evidence = paste(
          "DOI 10.1016/j.celrep.2026.117110 reports gene-expression",
          "libraries sequenced in batches of 40 on three NovaSeq 6000 S4",
          "lanes; ten batch labels are retained in final_rna_data.h5ad",
          "(Zenodo record 20834804)"
        ),
        technical_batch_use_status =
          "eligible for technical-batch sensitivity model",
        technical_batch_use_evidence = paste(
          "all ten source batches are complete in the RNA object and each",
          "spans 23-45 donors; use only as the common sensitivity batch key"
        ),
        assay = "10x Multiome gene expression",
        sequencing_platform = "Illumina NovaSeq 6000",
        chemistry = "10x Genomics Chromium Single Cell Multiome",
        cell_type_source_column =
          "cell_type;InN_subtype;ExN_subtype;leiden_2",
        cell_type_source_file =
          "final_rna_data.h5ad obs + neuron_subtype.csv",
        cell_type_original_label_semantics = paste(
          "original Catching et al. broad or neuron-subtype annotation;",
          "the most specific reported label is retained"
        ),
        cell_type_label_semantics =
          "preferred original-study label without atlas harmonization"
      ),
      required_filters = list()
    ),
    GSE168408 = list(
      input_root = "raw",
      input_file =
        "RNA-all_full-counts-and-downsampled-CPM.h5ad",
      expected_input_bytes = 4963815084,
      expected_input_sha256 =
        "5ca0b1ebe638f13e4302b721b80f41755a759c87f7fe43da6aeb7e8f9d0fd9ce",
      matrix_group = "X",
      expected_cells = 154748L,
      expected_features = 26747L,
      expected_nonzero_values = 399668744,
      expected_total_counts = 1020935295,
      expected_analysis_cells = 153473L,
      expected_analysis_donors = 27L,
      matrix_storage = "dgCMatrix",
      default_analysis_role = "reference",
      external_metadata = list(
        author_barcode_metadata = list(
          root = "raw",
          file = "RNA-all_BCs-meta-data.csv",
          sep = ",",
          key_column = "Herring_Author_Cell_ID",
          matrix_key_column = "Original_Cell_ID",
          match_column = "Herring_Author_Metadata_Matched",
          require_all_matrix_cells = TRUE,
          require_all_metadata_cells = TRUE,
          transform = function(meta) {
            required <- c(
              "batch", "RL#", "age", "chem", "numerical_age",
              "stage_id", "Sex", "Brain Regions*", "cell_type",
              "major_clust", "sub_clust", "leiden"
            )
            missing <- setdiff(required, names(meta))
            if (length(missing) > 0L) {
              stop(
                "Herring author barcode metadata is missing columns: ",
                paste(missing, collapse = ", ")
              )
            }
            first_column <- names(meta)[[1L]]
            cell_id <- as.character(meta[[first_column]])
            result <- data.frame(
              Herring_Author_Cell_ID = cell_id,
              stringsAsFactors = FALSE,
              check.names = FALSE
            )
            for (column in required) {
              result[[paste0("Herring_CSV_", column)]] <- meta[[column]]
            }
            result
          }
        ),
        geo_library_metadata = list(
          root = "raw",
          file = "source_library_metadata.tsv",
          sep = "\t",
          key_column = "Source_RL_ID",
          matrix_key_column = "RL#",
          match_column = "Herring_GEO_Library_Matched",
          require_all_matrix_cells = TRUE,
          require_all_metadata_cells = FALSE,
          transform = function(meta) {
            required <- c(
              "Source_RL_ID", "Source_GEO_RNA_Record_ID",
              "Source_GEO_Sample_Title", "Source_GEO_Platform_ID",
              "Source_Sequencing_Platform",
              "Source_Filtered_Matrix_File"
            )
            missing <- setdiff(required, names(meta))
            if (length(missing) > 0L) {
              stop(
                "Herring GEO library metadata is missing columns: ",
                paste(missing, collapse = ", ")
              )
            }
            meta[, required, drop = FALSE]
          }
        )
      ),
      raw_metadata_transform = function(meta) {
        required <- c(
          "Original_Cell_ID", "batch", "RL#", "age", "chem",
          "numerical_age", "stage_id", "Sex", "Brain Regions*",
          "cell_type", "major_clust", "sub_clust", "leiden",
          "Herring_Author_Metadata_Matched",
          "Herring_GEO_Library_Matched",
          "Source_GEO_RNA_Record_ID", "Source_GEO_Platform_ID",
          "Source_Sequencing_Platform"
        )
        missing <- setdiff(required, names(meta))
        if (length(missing) > 0L) {
          stop(
            "Herring joined source metadata is missing columns: ",
            paste(missing, collapse = ", ")
          )
        }
        if (any(!meta$Herring_Author_Metadata_Matched) ||
          any(!meta$Herring_GEO_Library_Matched)) {
          stop("Herring source metadata did not match every nucleus")
        }

        compared_fields <- c(
          "batch", "RL#", "age", "chem", "stage_id", "Sex",
          "Brain Regions*", "cell_type", "major_clust", "sub_clust",
          "leiden"
        )
        for (column in compared_fields) {
          csv_column <- paste0("Herring_CSV_", column)
          if (!identical(
            as.character(meta[[column]]),
            as.character(meta[[csv_column]])
          )) {
            stop(
              "Herring H5AD obs differs from author barcode CSV: ",
              column
            )
          }
        }
        numerical_age <- suppressWarnings(as.numeric(meta$numerical_age))
        csv_numerical_age <- suppressWarnings(as.numeric(
          meta$Herring_CSV_numerical_age
        ))
        if (anyNA(numerical_age) || anyNA(csv_numerical_age) ||
          any(abs(numerical_age - csv_numerical_age) > 1e-12)) {
          stop("Herring numerical age differs from author barcode CSV")
        }
        if (nrow(meta) != 154748L ||
          length(unique(meta$batch)) != 27L ||
          length(unique(meta[["RL#"]])) != 27L ||
          sum(meta$cell_type == "Poor-Quality") != 1275L ||
          !setequal(unique(meta$chem), c("v2", "v3"))) {
          stop("Herring author object dimensions or labels changed")
        }

        source_age <- as.character(meta$age)
        prenatal <- grepl(
          "^ga[0-9]+(?:[.][0-9]+)?$",
          source_age,
          perl = TRUE
        )
        postnatal_days <- grepl(
          "^[0-9]+(?:[.][0-9]+)?d$",
          source_age,
          perl = TRUE
        )
        postnatal_years <- grepl(
          "^[0-9]+(?:[.][0-9]+)?yr$",
          source_age,
          perl = TRUE
        )
        if (any(!(prenatal | postnatal_days | postnatal_years))) {
          stop("Herring contains an unreviewed source age label")
        }
        source_number <- suppressWarnings(as.numeric(sub(
          "^(?:ga)?([0-9]+(?:[.][0-9]+)?)(?:d|yr)?$",
          "\\1",
          source_age,
          perl = TRUE
        )))
        if (anyNA(source_number)) {
          stop("Herring source age cannot be converted")
        }
        canonical_value <- ifelse(
          prenatal,
          source_number - 2,
          ifelse(postnatal_days, source_number / 365, source_number)
        )
        canonical_unit <- ifelse(prenatal, "PCW", "years")
        meta$Age_Source_Raw <- source_age
        meta$Age_Source_Unit <- ifelse(
          prenatal,
          "gestational weeks",
          ifelse(postnatal_days, "postnatal days", "postnatal years")
        )
        meta$Age_Source_Basis <- ifelse(
          prenatal,
          "gestational age",
          "postnatal age at death"
        )
        meta$Age_Source_Reference <- paste(
          "Herring et al. Cell 2022 author barcode metadata;",
          "corrected public metadata dated 2023-01-18"
        )
        meta$Age_Harmonization_Input <- paste(
          format(
            canonical_value,
            scientific = FALSE,
            trim = TRUE,
            digits = 10
          ),
          canonical_unit
        )
        meta$Age_Conversion_Formula <- ifelse(
          prenatal,
          "reported gestational weeks - 2 weeks = PCW",
          ifelse(
            postnatal_days,
            "reported postnatal days / 365 = years",
            "reported postnatal years retained"
          )
        )
        meta$Age_Conversion_Confidence <-
          "high: author age label and explicit unit"
        meta$Age_Conversion_Applied <- prenatal | postnatal_days

        region_map <- data.frame(
          source = c("BA8", "BA9", "BA10", "BA46"),
          label = c(
            "Brodmann (1909) area 8", "Brodmann (1909) area 9",
            "Brodmann (1909) area 10", "Brodmann (1909) area 46"
          ),
          ontology_id = c(
            "UBERON:0013539", "UBERON:0013540",
            "UBERON:0013541", "UBERON:0006483"
          ),
          stringsAsFactors = FALSE
        )
        region_index <- match(meta[["Brain Regions*"]], region_map$source)
        if (anyNA(region_index)) {
          stop("Herring contains an unreviewed Brodmann area")
        }
        meta$Herring_Brain_Region <- region_map$label[region_index]
        meta$Herring_Brain_Region_Ontology_Label <-
          region_map$label[region_index]
        meta$Herring_Brain_Region_Ontology_ID <-
          region_map$ontology_id[region_index]
        meta$brain_region_anatomical_level <- "Brodmann cortical area"
        meta$brain_region_mapping_method <- paste(
          "exact source Brodmann-area label mapped to UBERON",
          "2026-04-01"
        )
        meta$brain_region_mapping_confidence <- "high"

        meta$Herring_Technical_Batch_Source <- paste(
          meta$Source_GEO_Platform_ID,
          meta$chem,
          sep = "_"
        )
        batch_donors <- tapply(
          meta[["RL#"]],
          meta$Herring_Technical_Batch_Source,
          function(value) length(unique(value))
        )
        if (length(batch_donors) != 3L || any(batch_donors < 2L)) {
          stop(
            "Herring platform-chemistry levels do not span the expected",
            " donors"
          )
        }
        meta
      },
      field_map = list(
        cell_id = "Original_Cell_ID",
        donor_id = "RL#",
        sample_id = "batch",
        specimen_id = "batch",
        library_id = "batch",
        source_record_id = "Source_GEO_RNA_Record_ID",
        technical_batch_id = "Herring_Technical_Batch_Source",
        species = NULL,
        species_raw = NULL,
        age = "Age_Harmonization_Input",
        sex = "Sex",
        brain_region = "Herring_Brain_Region",
        brain_region_source = "Brain Regions*",
        brain_region_ontology_label =
          "Herring_Brain_Region_Ontology_Label",
        brain_region_ontology_id = "Herring_Brain_Region_Ontology_ID",
        diagnosis = NULL,
        cell_type_original = "sub_clust",
        cell_type = "sub_clust",
        cell_type_level_1 = "cell_type",
        cell_type_level_2 = "major_clust",
        cell_type_level_3 = "sub_clust",
        cell_type_cluster_id = "leiden",
        cell_type_ontology_label = NULL,
        cell_type_ontology_id = NULL,
        technology = NULL,
        modality = NULL,
        assay = NULL,
        sequencing_platform = "Source_Sequencing_Platform",
        chemistry = "chem"
      ),
      constants = list(
        species = "Homo sapiens",
        diagnosis = "Neurotypical",
        technology = "10x Genomics Chromium",
        modality = "snRNA-seq",
        assay = "10x 3-prime single-nucleus RNA-seq",
        donor_id_semantics = "author RL-number brain donor",
        specimen_id_semantics =
          "author prefrontal-cortex tissue sample",
        specimen_id_verification_status = "source reported",
        library_id_semantics = paste(
          "author 10x gene-expression library identified by donor,",
          "age label and chemistry"
        ),
        library_id_verification_status = paste(
          "verified: Herring et al. author H5AD, publication code and GEO"
        ),
        library_id_evidence = paste(
          "the author batch field follows RL-number_age_chemistry and maps",
          "one-to-one to 27 RNA libraries and GEO RNA records"
        ),
        source_record_id_semantics = "GEO snRNA-seq sample accession",
        technical_batch_semantics = paste(
          "cross-donor sequencing-platform and 10x chemistry combination"
        ),
        technical_batch_verification_status = paste(
          "verified: author chemistry field and GEO platform record"
        ),
        technical_batch_evidence = paste(
          "author metadata retains v2/v3 chemistry; GEO assigns each RNA",
          "library to GPL21697 NextSeq 550 or GPL24676 NovaSeq 6000"
        ),
        technical_batch_derivation_rule = paste(
          "Source_GEO_Platform_ID + '_' + author chem"
        ),
        technical_batch_use_status = paste(
          "eligible for technical-batch sensitivity model"
        ),
        technical_batch_use_evidence = paste(
          "all three platform-chemistry combinations span at least two",
          "donors; retain Dataset separately as the primary study batch"
        ),
        brain_region_ontology_source =
          "UBERON 2026-04-01",
        cell_type_source_column =
          "cell_type;major_clust;sub_clust;leiden",
        cell_type_source_file = paste(
          "Lister Lab RNA-all_full-counts-and-downsampled-CPM.h5ad obs",
          "+ RNA-all_BCs-meta-data.csv"
        ),
        cell_type_original_label_semantics = paste(
          "original Herring et al. final subcluster annotation"
        ),
        cell_type_label_semantics = paste(
          "preferred finest author annotation; Poor-Quality remains in the",
          "complete object but is excluded from reference analysis"
        )
      ),
      required_filters = list(
        list(
          column = "cell_type",
          values = c("Non-Neu", "PN", "IN"),
          reason = "author final label is Poor-Quality",
          excluded_role = "quality_excluded"
        )
      )
    ),
    Wang_2025 = list(
      input_root = "raw",
      input_file = "Wang_2025_RNA.h5ad",
      expected_input_bytes = 2782116565,
      matrix_group = "raw/X",
      expected_cells = 232328L,
      expected_features = 35477L,
      expected_analysis_cells = 196466L,
      expected_analysis_donors = 23L,
      default_analysis_role = "reference",
      raw_metadata_transform = function(meta) {
        required_qc <- c(
          "nCount_RNA", "nFeature_RNA", "ATAC_fragments_in_peaks",
          "TSS.enrichment", "Nucleosome_signal",
          "Scrublet_doublet_score"
        )
        missing_qc <- setdiff(required_qc, names(meta))
        if (length(missing_qc) > 0L) {
          stop(
            "Wang source QC fields are missing: ",
            paste(missing_qc, collapse = ", ")
          )
        }
        qc_numeric <- lapply(meta[required_qc], function(value) {
          suppressWarnings(as.numeric(value))
        })
        if (any(vapply(qc_numeric, anyNA, logical(1)))) {
          stop("Wang source QC fields contain non-numeric values")
        }
        if (any(qc_numeric$nCount_RNA < 1000 |
          qc_numeric$nCount_RNA > 25000) ||
          any(qc_numeric$nFeature_RNA <= 400) ||
          any(qc_numeric$ATAC_fragments_in_peaks < 100 |
            qc_numeric$ATAC_fragments_in_peaks > 100000) ||
          any(qc_numeric$TSS.enrichment <= 1) ||
          any(qc_numeric$Nucleosome_signal >= 2) ||
          any(qc_numeric$Scrublet_doublet_score > 0.3)) {
          stop("Wang retained source object violates a published QC bound")
        }
        source_age <- as.character(meta$development_stage)
        postconception_days <- suppressWarnings(as.numeric(
          meta$Estimated_postconceptional_age_in_days
        ))
        prenatal <- !is.na(source_age) & grepl(
          "post-fertilization|Carnegie stage",
          source_age,
          ignore.case = TRUE
        )
        if (any(prenatal & is.na(postconception_days))) {
          stop(
            "Wang prenatal rows require estimated postconceptional age in days"
          )
        }
        harmonization_input <- source_age
        harmonization_input[prenatal] <- paste(
          format(
            postconception_days[prenatal] / 7,
            scientific = FALSE,
            trim = TRUE,
            digits = 8
          ),
          "PCW"
        )
        meta$Age_Source_Raw <- source_age
        meta$Age_Source_Unit <- ifelse(
          prenatal,
          "postconceptional days",
          "controlled-vocabulary postnatal stage"
        )
        meta$Age_Source_Basis <- ifelse(
          prenatal,
          "postconceptional age",
          "postnatal age"
        )
        meta$Age_Harmonization_Input <- harmonization_input
        meta$Age_Conversion_Formula <- ifelse(
          prenatal,
          "Estimated_postconceptional_age_in_days / 7 = PCW",
          "explicit month/year value parsed from development_stage"
        )
        meta$Age_Conversion_Confidence <- ifelse(
          prenatal,
          "source numeric postconceptional days",
          "source controlled-vocabulary developmental stage"
        )
        meta$Age_Conversion_Applied <- prenatal
        explicit_source_label <- function(value) {
          value <- as.character(value)
          value[value == "Unknown"] <- "Unknown (source-reported)"
          value
        }
        meta$Source_Cell_Type_Class <- explicit_source_label(meta$Class)
        meta$Source_Cell_Type_Subclass <-
          explicit_source_label(meta$Subclass)
        meta$Source_Cell_Type_Type <- explicit_source_label(meta$Type)
        meta$Source_Cell_Type_Type_Updated <-
          explicit_source_label(meta$Type_updated)
        meta$Source_Cell_Type_Cluster <- as.character(meta$Cluster)
        meta$Source_Cell_Type_Unknown <-
          as.character(meta$Type_updated) == "Unknown"
        meta$Source_Scrublet_Threshold <- 0.3
        meta$Source_Scrublet_Class <- "singlet"
        meta$source_cell_type_doublet_flag <- FALSE
        meta$Source_Post_QC_Before_Neocortex_Filter_Nuclei <- 243535L
        meta$Source_Upstream_Non_Neocortical_Excluded_Nuclei <- 11207L
        meta$Source_Final_Neocortical_Nuclei <- 232328L
        meta$Genome_Transcriptome_Build <-
          "GRCh38; GENCODE v32 / Ensembl 98"
        meta$Source_Alignment_Pipeline <- "Cell Ranger ARC v2.0.0"
        source_tissue <- as.character(meta$tissue)
        brodmann_area <- grepl("^Brodmann", source_tissue)
        developmental_division <- source_tissue %in%
          c("telencephalon", "forebrain")
        meta$brain_region_anatomical_level <- ifelse(
          brodmann_area,
          "Brodmann area",
          ifelse(
            developmental_division,
            "developmental brain division",
            "cortical region"
          )
        )
        meta$brain_region_mapping_method <-
          "exact released CELLxGENE tissue ontology mapping retained"
        meta$brain_region_mapping_confidence <-
          "source-reported ontology mapping"
        meta$brain_region_source_column <- paste(
          "CELLxGENE tissue; source Region retained as the broader",
          "study region group"
        )
        meta$Source_Region_Group <- as.character(meta$Region)
        meta
      },
      field_map = list(
        cell_id = "Original_Cell_ID",
        donor_id = "donor_id",
        sample_id = "Sample_ID",
        specimen_id = "Sample_ID",
        library_id = "Sample_ID",
        source_record_id = "Sample_ID",
        technical_batch_id = "Sample_ID",
        species = NULL,
        species_raw = NULL,
        age = "Age_Harmonization_Input",
        sex = "sex",
        brain_region = "tissue",
        brain_region_source = "tissue",
        brain_region_ontology_label = "tissue",
        brain_region_ontology_id = "tissue_ontology_term_id",
        diagnosis = "disease",
        cell_type_original = "Source_Cell_Type_Type_Updated",
        cell_type = "Source_Cell_Type_Type_Updated",
        cell_type_level_1 = "Source_Cell_Type_Class",
        cell_type_level_2 = "Source_Cell_Type_Subclass",
        cell_type_level_3 = "Source_Cell_Type_Type_Updated",
        cell_type_cluster_id = "Source_Cell_Type_Cluster",
        cell_type_ontology_label = "cell_type",
        cell_type_ontology_id = "cell_type_ontology_term_id",
        technology = NULL,
        modality = NULL,
        assay = "assay",
        sequencing_platform = NULL,
        chemistry = NULL
      ),
      constants = list(
        species = "Homo sapiens",
        technology = "10x Genomics Multiome",
        modality = "snRNA-seq",
        donor_id_semantics = "source-reported donor",
        specimen_id_semantics = "source-reported tissue sample",
        specimen_id_verification_status = "source reported",
        library_id_semantics = paste(
          "sample-specific 10x Multiome ATAC plus Gene Expression reaction",
          "and pooled gene-expression library"
        ),
        library_id_verification_status =
          "verified: Nature 2025 methods and source sample metadata",
        library_id_evidence = paste(
          "DOI 10.1038/s41586-024-08351-7 reports one 10x Multiome",
          "reaction per sample, followed by pooling of libraries from",
          "individual samples; CELLxGENE retains the exact Sample_ID"
        ),
        source_record_id_semantics = "paired-multiome source sample record",
        technical_batch_semantics =
          "sample-specific 10x Multiome reaction and RNA library",
        technical_batch_verification_status =
          "verified: Nature 2025 methods and source sample metadata",
        technical_batch_evidence = paste(
          "DOI 10.1038/s41586-024-08351-7 reports 10,000 nuclei targeted",
          "per sample per reaction and libraries from individual samples",
          "pooled for NovaSeq 6000 sequencing"
        ),
        technical_batch_use_status =
          "record only: ineligible for technical-batch correction",
        technical_batch_use_evidence = paste(
          "all 38 source Sample_ID levels, and all 32 retained levels after",
          "donor de-duplication, contain one donor; the cross-donor gate fails"
        ),
        brain_region_ontology_source = "CELLxGENE/UBERON",
        sequencing_platform = "Illumina NovaSeq 6000",
        chemistry = paste(
          "10x Genomics Chromium Next GEM Single Cell Multiome ATAC plus",
          "Gene Expression"
        ),
        cell_type_source_column = paste(
          "Class;Subclass;Type;Type_updated;Cluster;cell_type;",
          "cell_type_ontology_term_id",
          sep = ""
        ),
        cell_type_source_file = "Wang_2025_RNA.h5ad obs",
        cell_type_original_label_semantics =
          "original Wang et al. Type_updated annotation",
        cell_type_label_semantics =
          "preferred original-study label without atlas harmonization"
      ),
      required_filters = list()
    ),
    GSE294786 = list(
      input_root = "processed",
      input_file = "GSE294786_all_counts_full.h5ad",
      matrix_group = "X",
      expected_cells = 37994L,
      expected_features = 23841L,
      expected_analysis_cells = 11842L,
      expected_analysis_donors = 5L,
      default_analysis_role = "reference",
      external_metadata = list(
        root = "raw",
        file = "GSE294786_all_meta_progenitor.tsv.gz",
        key_column = 1L,
        require_all_matrix_cells = FALSE,
        require_all_metadata_cells = TRUE
      ),
      raw_metadata_transform = function(meta) {
        species_map <- c(
          Hg = "Homo sapiens",
          Pt = "Pan troglodytes",
          rheMac8 = "Macaca mulatta"
        )
        meta$Species_standardized <- unname(
          species_map[as.character(meta$Species)]
        )
        paper_age_value <- c(
          `1452` = 2, `1490` = 4, `1623` = 4,
          `1350` = 24, `1385` = 36,
          Chimp2 = 0, Chimp3 = 5, Frits = 50, Ruben = 18,
          Bart = 35, R12121 = 0, R17043 = 0
        )
        paper_age_unit <- c(
          `1452` = "months", `1490` = "months", `1623` = "months",
          `1350` = "years", `1385` = "years",
          Chimp2 = "months", Chimp3 = "months", Frits = "years",
          Ruben = "years", Bart = "years", R12121 = "months",
          R17043 = "months"
        )
        paper_pmd_hours <- c(
          `1452` = 36, `1490` = 29, `1623` = 20,
          `1350` = 36, `1385` = 22,
          Chimp2 = 5, Chimp3 = 5, Frits = 7, Ruben = 7,
          Bart = 7, R12121 = 1, R17043 = 1
        )
        paper_tissue_source <- c(
          `1452` = "HBCC", `1490` = "HBCC", `1623` = "HBCC",
          `1350` = "HBCC", `1385` = "HBCC",
          Chimp2 = "BPRC", Chimp3 = "BPRC", Frits = "BPRC",
          Ruben = "BPRC", Bart = "BPRC", R12121 = "BPRC",
          R17043 = "BPRC"
        )
        paper_cause_of_death <- c(
          `1452` = "asphyxia", `1490` = "asphyxia",
          `1623` = "asphyxia", `1350` = "bronchial asthma",
          `1385` = "blunt force injury (chest)",
          Chimp2 = "euthanasia due to parental neglect",
          Chimp3 = "euthanasia due to parental neglect",
          Frits = "cardiomyopathy",
          Ruben = "cardiomyopathy, pneumonia",
          Bart = "spontaneous death",
          R12121 = "euthanasia due to parental neglect",
          R17043 = "euthanasia due to parental neglect"
        )
        sample_id <- as.character(meta$SampleID)
        age_value <- unname(paper_age_value[sample_id])
        age_unit <- unname(paper_age_unit[sample_id])
        paper_age <- ifelse(
          is.na(age_value) | is.na(age_unit),
          NA_character_,
          paste(age_value, age_unit)
        )
        metadata_age_days <- suppressWarnings(as.numeric(meta$Age2))
        paper_age_days <- ifelse(
          age_unit == "months",
          age_value * (365 / 12),
          age_value * 365
        )
        age_agrees <- !is.na(metadata_age_days) &
          !is.na(paper_age_days) &
          abs(metadata_age_days - paper_age_days) <= 2
        age_conflict <- !is.na(metadata_age_days) &
          !is.na(paper_age_days) &
          !age_agrees
        meta$Source_Metadata_Age2_Days <- as.character(meta$Age2)
        meta$Source_Paper_Table_S1_Age <- paper_age
        meta$Source_Paper_Table_S1_Sex <- ifelse(
          is.na(sample_id),
          NA_character_,
          "Male"
        )
        meta$Source_Paper_Table_S1_PMD_Hours <- unname(
          paper_pmd_hours[sample_id]
        )
        meta$Source_Paper_Table_S1_Tissue_Source <- unname(
          paper_tissue_source[sample_id]
        )
        meta$Source_Paper_Table_S1_Cause_of_Death <- unname(
          paper_cause_of_death[sample_id]
        )
        meta$Source_Age_Metadata_Paper_Agreement <- ifelse(
          is.na(paper_age),
          NA_character_,
          ifelse(age_agrees, "agree", "conflict")
        )
        meta$Age_Source_Raw <- paper_age
        meta$Age_Source_Unit <- age_unit
        meta$Age_Source_Basis <- "postnatal age"
        meta$Age_Harmonization_Input <- paper_age
        meta$Age_Conversion_Formula <- ifelse(
          is.na(age_unit),
          NA_character_,
          ifelse(
            age_unit == "months",
            "publication Table S1 months / 12 = years",
            "publication Table S1 reported years retained"
          )
        )
        meta$Age_Conversion_Confidence <- ifelse(
          age_conflict,
          paste(
            "verified publication Table S1 selected over conflicting",
            "GEO processed metadata Age2"
          ),
          "verified publication Table S1; agrees with GEO Age2"
        )
        meta$Age_Conversion_Applied <- !is.na(paper_age)
        meta$Source_Paper_Diagnosis <- ifelse(
          is.na(sample_id),
          NA_character_,
          "No diagnosed developmental defect or disease"
        )
        meta$Source_Original_Cell_Type <- meta$CellType
        meta$Source_Original_Cell_Type_2 <- meta$CellType2
        meta$Source_Integrated_Cell_Type <- meta$cell_type
        meta$Source_Final_Cell_Type <- meta$cell_type_final
        meta$Source_Final_Cell_Type_Broad <- meta$cell_type_final_broad
        meta$Source_Integrated_Seurat_Cluster <-
          meta$`integrated_snn_res.0.8`
        meta$Source_Seurat_Cluster <- meta$seurat_clusters
        meta$Source_Leiden_Cluster <- meta$leiden
        meta$Source_Leiden_0_5_Cluster <- meta$leiden_0_5
        meta$Genome_Transcriptome_Build <-
          "humanized consensus genome in hg38 coordinates"
        meta
      },
      field_map = list(
        cell_id = "Original_Cell_ID",
        donor_id = "SampleID",
        sample_id = "SampleID",
        specimen_id = "SampleID",
        library_id = "orig.ident",
        source_record_id = "orig.ident",
        technical_batch_id = "orig.ident",
        species = "Species_standardized",
        species_raw = "Species",
        age = "Age_Harmonization_Input",
        sex = "Source_Paper_Table_S1_Sex",
        brain_region = NULL,
        brain_region_source = NULL,
        brain_region_ontology_label = NULL,
        brain_region_ontology_id = NULL,
        diagnosis = "Source_Paper_Diagnosis",
        cell_type_original = "Source_Final_Cell_Type",
        cell_type = "Source_Final_Cell_Type",
        cell_type_level_1 = "Source_Final_Cell_Type_Broad",
        cell_type_level_2 = "Source_Final_Cell_Type",
        cell_type_level_3 = NULL,
        cell_type_cluster_id = "Source_Leiden_Cluster",
        cell_type_ontology_label = NULL,
        cell_type_ontology_id = NULL,
        technology = NULL,
        modality = NULL,
        assay = NULL,
        sequencing_platform = NULL,
        chemistry = NULL
      ),
      constants = list(
        brain_region = "Dorsolateral prefrontal cortex",
        brain_region_source = "DLPFC, approximately Brodmann area 46/9",
        brain_region_ontology_label = "dorsolateral prefrontal cortex",
        brain_region_ontology_id = "UBERON:0009834",
        brain_region_ontology_source =
          "UBERON 2026-04-01 exact-label lookup",
        technology = "10x Genomics",
        modality = "snRNA-seq",
        assay = "Gene Expression",
        sequencing_platform = "Illumina HiSeq 2500",
        chemistry = "10x Chromium v3 Single Cell RNA Gene Expression",
        donor_id_semantics = "source-reported donor",
        specimen_id_semantics = "source-reported tissue sample",
        specimen_id_verification_status = "source reported",
        library_id_semantics = paste(
          "10x Chromium gene-expression library for one cross-species",
          "nuclei pool and sequencing sample"
        ),
        library_id_verification_status =
          "verified: Science Advances 2026 methods and GEO metadata",
        library_id_evidence = paste(
          "DOI 10.1126/sciadv.aea3316 states that nuclei from separate",
          "species were combined 1:1 per sequencing sample and loaded on",
          "10x Chromium; GEO GSE294786 retains the six exact sequencing-pool",
          "identifiers in orig.ident"
        ),
        source_record_id_semantics =
          "source cross-species sequencing-pool record",
        technical_batch_semantics = paste(
          "cross-species 10x processing and sequencing pool identified by",
          "orig.ident"
        ),
        technical_batch_verification_status =
          "verified: Science Advances 2026 methods and GEO metadata",
        technical_batch_evidence = paste(
          "DOI 10.1126/sciadv.aea3316 reports simultaneous 1:1 processing",
          "of different species per sequencing sample; GEO metadata maps",
          "all 37,994 nuclei to six orig.ident pool values"
        ),
        technical_batch_use_status =
          "record only: ineligible for technical-batch correction",
        technical_batch_use_evidence = paste(
          "the human-only analysis retains five pools but each contains",
          "exactly one human donor; cross-donor correction is therefore",
          "not identifiable after the required species filter"
        ),
        cell_type_source_column = paste(
          "CellType;CellType2;cell_type;cell_type_final_broad;",
          "cell_type_final;integrated_snn_res.0.8;seurat_clusters;",
          "leiden;leiden_0_5",
          sep = ""
        ),
        cell_type_source_file = "GSE294786 source metadata",
        cell_type_original_label_semantics =
          "original Klavert et al. cell_type_final annotation",
        cell_type_label_semantics =
          "preferred original-study label without atlas harmonization"
      ),
      required_filters = list(
        list(
          column = "Species",
          values = "Homo sapiens",
          reason = paste(
            "non-human species excluded from the human-only reference;",
            "all species remain auditable in the complete full object"
          ),
          missing_reason = paste(
            "source count row has no matching released cell metadata;",
            "retained only in the complete full object"
          ),
          excluded_role = "reference_excluded"
        )
      )
    ),
    Velmeshev_2023 = list(
      input_root = "raw",
      input_file = "Velmeshev_2023_full_cellxgene.h5ad",
      expected_input_bytes = 5976097241,
      matrix_group = "raw/X",
      expected_cells = 709372L,
      expected_features = 17631L,
      expected_analysis_cells = 356566L,
      expected_analysis_donors = 59L,
      default_analysis_role = "reference",
      external_metadata = list(
        root = "raw",
        file = "Velmeshev_2023_ucsc_metadata.tsv",
        key_column = "UCSC_Cell_ID",
        require_all_matrix_cells = TRUE,
        require_all_metadata_cells = TRUE,
        transform = function(meta) {
          required <- c(
            "Cell_ID", "Age", "Age_(days)", "Age_Range",
            "Chemistry", "Dataset", "Individual", "PMI", "Region",
            "Region_Broad", "Sample", "Seurat_clusters", "Sex",
            "Sex_Original", "Lineage"
          )
          missing <- setdiff(required, names(meta))
          if (length(missing) > 0L) {
            stop(
              "Velmeshev UCSC metadata is missing columns: ",
              paste(missing, collapse = ", ")
            )
          }
          selected <- meta[, required, drop = FALSE]
          names(selected) <- paste0(
            "UCSC_",
            c(
              "Cell_ID", "Age", "Age_Days", "Age_Range",
              "Chemistry", "Dataset", "Individual", "PMI", "Region",
              "Region_Broad", "Sample", "Seurat_Clusters", "Sex",
              "Sex_Original", "Lineage"
            )
          )
          selected
        }
      ),
      raw_metadata_transform = function(meta) {
        required <- c(
          "External_Metadata_Matched", "dataset", "donor_id", "sample",
          "UCSC_Age", "UCSC_Age_Days", "UCSC_Age_Range",
          "UCSC_Dataset", "UCSC_Individual", "UCSC_Sample",
          "UCSC_PMI", "UCSC_Region", "UCSC_Sex",
          "UCSC_Seurat_Clusters", "UCSC_Lineage"
        )
        missing <- setdiff(required, names(meta))
        if (length(missing) > 0L) {
          stop(
            "Velmeshev joined metadata is missing columns: ",
            paste(missing, collapse = ", ")
          )
        }
        if (any(!meta$External_Metadata_Matched)) {
          stop("Velmeshev UCSC metadata did not match every matrix cell")
        }
        source_dataset <- as.character(meta$dataset)
        source_donor <- as.character(meta$donor_id)
        source_sample <- as.character(meta$sample)
        external_dataset <- as.character(meta$UCSC_Dataset)
        external_donor <- as.character(meta$UCSC_Individual)
        external_sample <- as.character(meta$UCSC_Sample)
        if (!identical(source_dataset, external_dataset)) {
          stop("Velmeshev UCSC Dataset differs from the H5AD dataset field")
        }
        if (!identical(source_donor, external_donor)) {
          stop("Velmeshev UCSC Individual differs from the H5AD donor_id")
        }
        if (!identical(source_sample, external_sample)) {
          stop("Velmeshev UCSC Sample differs from the H5AD sample field")
        }
        if (!identical(
          as.character(meta$region),
          as.character(meta$UCSC_Region)
        )) {
          stop("Velmeshev UCSC Region differs from the H5AD region field")
        }
        if (!identical(
          tolower(as.character(meta$sex)),
          tolower(as.character(meta$UCSC_Sex))
        )) {
          stop("Velmeshev UCSC Sex differs from the H5AD sex field")
        }
        if (!identical(
          as.character(meta$PMI),
          as.character(meta$UCSC_PMI)
        )) {
          stop("Velmeshev UCSC PMI differs from the H5AD PMI field")
        }

        expected_dataset_counts <- c(
          Velmeshev = 358663L,
          Ramos = 186980L,
          Herring = 154748L,
          Trevino = 8981L
        )
        observed_dataset_counts <- table(source_dataset)
        observed_dataset_counts <- observed_dataset_counts[
          names(expected_dataset_counts)
        ]
        if (anyNA(observed_dataset_counts) ||
          !identical(
            as.integer(observed_dataset_counts),
            unname(expected_dataset_counts)
          )) {
          stop(
            "Velmeshev source-dataset counts differ from the reviewed ",
            "complete object"
          )
        }
        velmeshev_rows <- source_dataset == "Velmeshev"
        velmeshev_donors <- unique(stats::na.omit(
          source_donor[velmeshev_rows]
        ))
        if (length(velmeshev_donors) != 60L) {
          stop("Velmeshev source subset must contain 60 donors")
        }
        duplicate_donor_cells <- sum(
          velmeshev_rows & source_donor == "5936",
          na.rm = TRUE
        )
        if (duplicate_donor_cells != 2097L) {
          stop("Velmeshev donor 5936 must contain 2,097 source cells")
        }
        if (any(is.na(normalize_missing_metadata(
          meta$UCSC_Lineage
        )))) {
          stop("Velmeshev UCSC Lineage contains missing values")
        }

        source_age <- as.character(meta$UCSC_Age)
        source_age_days <- suppressWarnings(as.numeric(
          meta$UCSC_Age_Days
        ))
        if (anyNA(source_age_days)) {
          stop("Velmeshev UCSC Age_(days) contains non-numeric values")
        }
        prenatal <- grepl(
          "trimester",
          as.character(meta$UCSC_Age_Range),
          ignore.case = TRUE
        )
        if (any(!prenatal & source_age_days < 266)) {
          stop(
            "Velmeshev postnatal Age_(days) precedes the source birth offset"
          )
        }
        harmonized_value <- ifelse(
          prenatal,
          source_age_days / 7,
          (source_age_days - 266) / 365
        )
        harmonized_unit <- ifelse(prenatal, "PCW", "years")
        meta$Age_Source_Raw <- source_age
        meta$Age_Source_Unit <- ifelse(
          prenatal,
          paste(
            "reported obstetric gestational-week label; companion",
            "Age_(days) is days post-fertilization"
          ),
          paste(
            "reported postnatal age label; companion Age_(days)",
            "includes the 266-day source birth offset"
          )
        )
        meta$Age_Source_Basis <-
          "continuous post-fertilization age axis"
        meta$Age_Source_Days_Post_Fertilization <- source_age_days
        meta$Age_Source_Birth_Offset_Days <- 266
        meta$Age_Harmonization_Input <- paste(
          format(
            harmonized_value,
            scientific = FALSE,
            trim = TRUE,
            digits = 10
          ),
          harmonized_unit
        )
        meta$Age_Conversion_Formula <- ifelse(
          prenatal,
          "UCSC Age_(days) / 7 = PCW",
          "(UCSC Age_(days) - 266 source birth offset) / 365 = years"
        )
        meta$Age_Conversion_Confidence <-
          "paper-linked UCSC continuous post-fertilization day axis"
        meta$Age_Conversion_Applied <- TRUE
        detailed_region <- as.character(meta$UCSC_Region)
        brodmann_area <- grepl("^BA[0-9]+$", detailed_region)
        developmental_structure <- detailed_region %in%
          c("LGE", "CGE", "MGE", "GE")
        meta$brain_region_anatomical_level <- ifelse(
          brodmann_area,
          "Brodmann area; ontology records its parent cortical region",
          ifelse(
            developmental_structure,
            "developmental brain structure",
            paste(
              "source cortical region or subregion; ontology records",
              "the released parent region"
            )
          )
        )
        meta$brain_region_mapping_method <- paste(
          "source detailed region retained; released CELLxGENE",
          "parent-region ontology mapping retained"
        )
        meta$brain_region_mapping_confidence <- ifelse(
          developmental_structure,
          "source-reported ontology mapping",
          paste(
            "high for parent-region ontology; detailed source subregion",
            "is preserved separately"
          )
        )
        meta$Source_Original_Lineage <- meta$UCSC_Lineage
        meta$Source_Original_Seurat_Cluster <-
          meta$UCSC_Seurat_Clusters
        meta
      },
      field_map = list(
        cell_id = "Original_Cell_ID",
        donor_id = "donor_id",
        sample_id = "sample",
        specimen_id = "sample",
        library_id = "sample",
        source_record_id = "sample",
        technical_batch_id = "sample",
        species = NULL,
        species_raw = NULL,
        age = "Age_Harmonization_Input",
        sex = "sex",
        brain_region = "tissue",
        brain_region_source = "region",
        brain_region_ontology_label = "tissue",
        brain_region_ontology_id = "tissue_ontology_term_id",
        diagnosis = "disease",
        cell_type_original = "Source_Original_Lineage",
        cell_type = "Source_Original_Lineage",
        cell_type_level_1 = "Source_Original_Lineage",
        cell_type_level_2 = "cell_type",
        cell_type_level_3 = NULL,
        cell_type_cluster_id = "Source_Original_Seurat_Cluster",
        cell_type_ontology_label = "cell_type",
        cell_type_ontology_id = "cell_type_ontology_term_id",
        technology = NULL,
        modality = NULL,
        assay = "assay",
        sequencing_platform = NULL,
        chemistry = "UCSC_Chemistry"
      ),
      constants = list(
        species = "Homo sapiens",
        technology = "10x Genomics",
        modality = "snRNA-seq",
        donor_id_semantics = "source-reported donor",
        specimen_id_semantics = "source-reported tissue sample",
        specimen_id_verification_status = "source reported",
        library_id_semantics = paste(
          "single-nucleus 10x gene-expression library prepared from one",
          "source tissue sample"
        ),
        library_id_verification_status =
          "verified: Science 2023 methods and UCSC sample metadata",
        library_id_evidence = paste(
          "DOI 10.1126/science.adf0834 reports one 10x capture and library",
          "preparation per individual sample before pooling libraries for",
          "sequencing; UCSC metadata retains 108 exact source samples"
        ),
        source_record_id_semantics = "source sample record",
        technical_batch_semantics =
          "sample-specific 10x v2 capture and gene-expression library",
        technical_batch_verification_status =
          "verified: Science 2023 methods and UCSC sample metadata",
        technical_batch_evidence = paste(
          "the paper reports target capture per sample, unmodified 10x",
          "library preparation and pooling of individual-sample libraries",
          "for NovaSeq 6000 sequencing"
        ),
        technical_batch_use_status =
          "record only: ineligible for technical-batch correction",
        technical_batch_use_evidence = paste(
          "each of 108 source sample libraries, and each of 107 retained",
          "libraries after donor de-duplication, contains one donor; the",
          "cross-donor gate fails; Chemistry is V2 for the full Velmeshev",
          "source subset and therefore cannot define an internal batch"
        ),
        brain_region_ontology_source = "CELLxGENE/UBERON",
        sequencing_platform = "Illumina NovaSeq 6000",
        cell_type_source_column = paste(
          "UCSC Lineage;UCSC Seurat_clusters;CELLxGENE cell_type;",
          "CELLxGENE cell_type_ontology_term_id",
          sep = ""
        ),
        cell_type_source_file =
          "Velmeshev_2023_ucsc_metadata.tsv + CELLxGENE obs",
        cell_type_original_label_semantics =
          "original Velmeshev et al. UCSC Lineage annotation",
        cell_type_label_semantics =
          "preferred original-study label without atlas harmonization"
      ),
      required_filters = list(
        list(
          column = "dataset",
          values = "Velmeshev",
          reason = "source object row belongs to a previously integrated dataset",
          excluded_role = "reference_excluded"
        )
      )
    )
  )
  if (!dataset %in% names(configs)) {
    stop("No dataset preprocessing config for: ", dataset)
  }
  configs[[dataset]]
}

transform_raw_metadata <- function(meta, config) {
  transform <- config$raw_metadata_transform
  if (is.null(transform)) {
    return(meta)
  }
  transformed <- transform(meta)
  if (!is.data.frame(transformed) || nrow(transformed) != nrow(meta)) {
    stop("raw metadata transform must preserve all input rows")
  }
  transformed
}

transform_external_metadata <- function(meta, external_config) {
  transform <- external_config$transform
  if (is.null(transform)) {
    return(meta)
  }
  transformed <- transform(meta)
  if (!is.data.frame(transformed) || nrow(transformed) != nrow(meta)) {
    stop("external metadata transform must preserve all input rows")
  }
  transformed
}

apply_required_filters <- function(meta, config) {
  filters <- config$required_filters
  if (length(filters) == 0L) {
    return(meta)
  }
  for (filter in filters) {
    if (!filter$column %in% names(meta)) {
      stop("required filter column is missing: ", filter$column)
    }
    filter_value <- meta[[filter$column]]
    missing_value <- is.na(filter_value)
    keep <- !missing_value & filter_value %in% filter$values
    newly_excluded <- meta$Analysis_Include & !keep
    meta$Analysis_Include[newly_excluded] <- FALSE
    meta$Exclusion_Reason[newly_excluded] <- filter$reason
    if (!is.null(filter$missing_reason)) {
      meta$Exclusion_Reason[newly_excluded & missing_value] <-
        filter$missing_reason
    }
    if (!is.null(filter$excluded_role)) {
      meta$Analysis_Role[newly_excluded] <- filter$excluded_role
    }
  }
  meta
}
