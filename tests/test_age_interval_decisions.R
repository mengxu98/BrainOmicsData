source("functions/dataset_metadata.R")

meta <- data.frame(
  Dataset = c("GSE144136", "GSE67835"),
  Age = c("38 years", "14-16 PCW"),
  Age_Harmonization_Input = c("38 years", "14-16 PCW"),
  BrainRegion = c(
    "Dorsolateral prefrontal cortex",
    "Cerebral cortex"
  ),
  stringsAsFactors = FALSE
)

corrected <- add_age_schema(meta)
stopifnot(identical(corrected$Stage, c("S13", "S4")))
stopifnot(identical(corrected$Age, c("38.71 years", "14-16 PCW")))
stopifnot(identical(corrected$Age_num, c(38.71, 15)))
stopifnot(identical(corrected$Age_Source_Raw, meta$Age))

groups <- add_age_schema(data.frame(
  Dataset = c("GSE144136", "GSE144136", "Other"),
  Age = c("38.71 years", "41.06 years", "38 years"),
  Age_Source_Raw = c("controls 38.71 +/- 4.32 years", "MDD 41.06 +/- 4.66 years", "38 years"),
  stringsAsFactors = FALSE
))
stopifnot(identical(groups$Age_num, c(38.71, 41.06, 38)))
stopifnot(identical(groups$Stage, c("S13", "S14", "S13")))
stopifnot(all(groups$Continuous_Age_Eligibility[1:2] ==
  "not eligible: exact donor age is not available"))
stopifnot(all(groups$Age_Representative_Method[1:2] == "published diagnosis-group mean"))
stopifnot(identical(add_age_schema(groups)$Age_num, groups$Age_num))
