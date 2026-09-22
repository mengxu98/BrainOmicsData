age_provenance <- function(
  metadata,
  age,
  source_raw,
  source_unit,
  source_basis,
  source_reference,
  formula,
  confidence,
  converted = FALSE
) {
  metadata$Age <- age
  metadata$Age_Source_Raw <- source_raw
  metadata$Age_Source_Unit <- source_unit
  metadata$Age_Source_Basis <- source_basis
  metadata$Age_Source_Reference <- source_reference
  metadata$Age_Harmonization_Input <- age
  metadata$Age_Conversion_Formula <- formula
  metadata$Age_Conversion_Confidence <- confidence
  metadata$Age_Conversion_Applied <- converted
  metadata
}

base_source_metadata <- function(
  cells,
  donor,
  specimen,
  library,
  source_record,
  age,
  sex,
  region,
  diagnosis,
  cell_type,
  technology = "10x Genomics",
  sequence = "snRNA-seq"
) {
  data.frame(
    Cells = as.character(cells),
    Dataset = dataset,
    Technology = technology,
    Sequence = sequence,
    Assay_Type = sequence,
    Original_Cell_ID = as.character(cells),
    Original_Donor_ID = as.character(donor),
    Original_Sample_ID = as.character(specimen),
    Original_Specimen_ID = as.character(specimen),
    Original_Library_ID = as.character(library),
    Original_Source_Sample_ID = as.character(specimen),
    Original_Source_Record_ID = as.character(source_record),
    Sample = as.character(specimen),
    Sample_ID = as.character(specimen),
    Age = as.character(age),
    Sex = as.character(sex),
    Brain_Region = as.character(region),
    Region = as.character(region),
    Diagnosis = as.character(diagnosis),
    CellType_raw = as.character(cell_type),
    Source_Preprocessing_Retained_Object_Scope = paste(
      "complete publicly released author-processed RNA count matrices",
      "selected for this transcriptomic resource; every input matrix cell",
      "is retained in the reconstructed object"
    ),
    Source_Preprocessing_Attrition_Status = paste(
      "no cell was dropped during source reconstruction; cells outside the",
      "main analysis remain in the complete object and are identified by",
      "Analysis_Include and Exclusion_Reason"
    ),
    Source_Preprocessing_Diagnosis_Audit_Status = paste(
      "source diagnosis or study status is retained in Diagnosis; analysis",
      "eligibility is recorded separately without deleting source cells"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}
