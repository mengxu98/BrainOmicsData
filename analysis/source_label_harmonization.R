suppressPackageStartupMessages({
  library(data.table)
})

source("functions/data_paths.R")
source("functions/dataset_metadata.R")
source("functions/processed_object.R")
source("functions/integration.R")


integration_dir <- brainomics_data_path("integration_25")
annotation_dir <- file.path(integration_dir, "annotation")
input_file <- file.path(annotation_dir, "source_label_hierarchy.tsv.gz")
output_file <- file.path(annotation_dir, "source_label_harmonization.tsv")
summary_file <- file.path(
  annotation_dir,
  "source_label_harmonization_summary.tsv"
)
if (!file.exists(input_file)) {
  stop("Missing source-label hierarchy audit: ", input_file)
}

mapping <- as.data.table(read.delim(
  gzfile(input_file),
  sep = "\t",
  quote = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
))
required_fields <- c(
  "Dataset", "Source_Original_Label", "Source_Label", "Source_Level_1",
  "Source_Level_2", "Source_Level_3", "Source_Cluster_ID", "Cells"
)
missing_fields <- setdiff(required_fields, names(mapping))
if (length(missing_fields) > 0L) {
  stop("Source hierarchy is missing: ", paste(missing_fields, collapse = ", "))
}
source_hierarchy_fields <- c(
  "Source_Original_Label", "Source_Label", "Source_Level_1",
  "Source_Level_2", "Source_Level_3", "Source_Cluster_ID"
)
for (field in source_hierarchy_fields) {
  set(mapping, j = field, value = normalize_missing_metadata(mapping[[field]]))
}

normalise_text <- function(value) {
  value <- ifelse(is.na(value), "", as.character(value))
  tolower(trimws(value))
}

mapping[, Search_Text := paste(
  normalise_text(Source_Original_Label),
  normalise_text(Source_Label),
  normalise_text(Source_Level_1),
  normalise_text(Source_Level_2),
  normalise_text(Source_Level_3)
)]
mapping[, Primary_Text := paste(
  normalise_text(Source_Original_Label),
  normalise_text(Source_Label)
)]
# GSE186538 subclass is a separate classification, not an original_name parent.
# HYPOMAP's broad compound names (AstroEpendymal, Immune/Vascular) must not
# override the specific author label. Keep all source fields in the mapping table.
mapping[Dataset %chin% c("GSE186538", "HYPOMAP"), Search_Text := Primary_Text]
mapping[, `:=`(
  Harmonized_CellClass = NA_character_,
  Harmonized_CellType = NA_character_,
  Mapping_Status = NA_character_,
  Mapping_Rule = NA_character_
)]

assign_mapping <- function(index, cell_class, cell_type, status, rule) {
  index <- index & is.na(mapping$Mapping_Status)
  if (!any(index)) {
    return(invisible(NULL))
  }
  mapping[index, `:=`(
    Harmonized_CellClass = cell_class,
    Harmonized_CellType = cell_type,
    Mapping_Status = status,
    Mapping_Rule = rule
  )]
  invisible(NULL)
}

has <- function(pattern) {
  grepl(pattern, mapping$Search_Text, perl = TRUE)
}
has_primary <- function(pattern) {
  grepl(pattern, mapping$Primary_Text, perl = TRUE)
}

# Explicit author hierarchies take precedence over ambiguous subtype names.
assign_mapping(mapping$Dataset == "GSE186538" & has_primary("(^| )aendo "),
  "Vascular or mesenchymal", "Endothelial cell",
  "mapped_for_common_comparison", "author aEndo identity")
assign_mapping(mapping$Dataset == "GSE186538" & mapping$Source_Label == "EC L6b TLE4 CCN2",
  "Excitatory neuronal", "Excitatory neuron",
  "mapped_for_common_comparison", "author entorhinal L6b excitatory identity")
assign_mapping(mapping$Dataset == "GSE186538" & has_primary("(^| )macro "),
  "Myeloid or immune", "Macrophage",
  "mapped_for_common_comparison", "author Macro identity; anatomical subtype not inferred")
assign_mapping(mapping$Dataset == "HYPOMAP" & has_primary("oligo-precursor"),
  "Oligodendrocyte lineage", "Oligodendrocyte precursor cell",
  "mapped_for_common_comparison", "author Oligo-Precursor identity")
assign_mapping(mapping$Dataset == "HYPOMAP" & has_primary("ependymal tanycytes"),
  "Epithelial or barrier", "Tanycyte",
  "mapped_for_common_comparison", "author specific tanycyte identity")
assign_mapping(mapping$Dataset == "HYPOMAP" & has_primary("ependymal choroid"),
  "Epithelial or barrier", "Choroid plexus cell",
  "mapped_for_common_comparison", "author specific choroid identity")
assign_mapping(mapping$Dataset == "HYPOMAP" & has_primary("ependymal ependymocytes"),
  "Epithelial or barrier", "Ependymal cell",
  "mapped_for_common_comparison", "author specific ependymocyte identity")
assign_mapping(mapping$Dataset == "HYPOMAP" & has_primary("immune macrophages"),
  "Myeloid or immune", "Macrophage",
  "mapped_for_common_comparison", "author macrophage identity; anatomical subtype not inferred")
assign_mapping(mapping$Dataset == "HYPOMAP" & has_primary("vascular smcs"),
  "Vascular or mesenchymal", "Vascular smooth muscle cell",
  "mapped_for_common_comparison", "author vascular SMC identity")
assign_mapping(mapping$Dataset == "HYPOMAP" & has_primary("immune myocytes"),
  "Other non-neuronal", "Source-reported Immune Myocytes",
  "excluded_incomparable_source_resolution", "compound author identity retained pending semantic verification")
assign_mapping(mapping$Dataset == "Wang_2025" & mapping$Source_Label == "IN-Mix-LAMP5" &
  mapping$Source_Level_2 == "GABAergic neuron",
  "Inhibitory neuronal", "Inhibitory neuron", "mapped_for_common_comparison",
  "author GABAergic superclass; Mix does not indicate a doublet")
assign_mapping(mapping$Dataset == "Li_et_al_2018" & mapping$Source_Label == "Unassigned" &
  mapping$Source_Level_1 == "ExN",
  "Excitatory neuronal", "Excitatory neuron", "mapped_for_common_comparison", "author ExN broad identity")
assign_mapping(mapping$Dataset == "Li_et_al_2018" & mapping$Source_Label == "Unassigned" &
  mapping$Source_Level_1 == "InN",
  "Inhibitory neuronal", "Inhibitory neuron", "mapped_for_common_comparison", "author InN broad identity")
assign_mapping(mapping$Dataset == "GSE168408" & mapping$Source_Label == "Micro_out" &
  mapping$Source_Level_2 == "Micro",
  "Myeloid or immune", "Microglia", "mapped_for_common_comparison", "Herring methods explicitly assign subcluster 53 to Microglia")
assign_mapping(mapping$Dataset == "GSE168408" & mapping$Source_Label == "OPC_MBP" &
  mapping$Source_Level_2 == "Oligo",
  "Oligodendrocyte lineage", "Oligodendrocyte", "mapped_for_common_comparison", "Herring methods explicitly assign cluster 49 to Oligodendrocytes")
assign_mapping(mapping$Dataset == "GSE104276" & mapping$Source_Level_2 == "Excitatory neurons",
  "Excitatory neuronal", "Excitatory neuron", "mapped_for_common_comparison", "exact cell ID in author excitatory neuron worksheet")

# Source-reported uncertainty and mixed/low-quality groups remain explicit and
# are never forced into a biological class for integration evaluation.
assign_mapping(
  has_primary("unknown|unassigned|n-undef") |
    mapping$Source_Label %chin% c("OUT", "UD"),
  "Mixed or unknown", "Source-reported unknown",
  "excluded_source_unknown", "source explicitly reports unknown or OUT"
)
assign_mapping(
  has_primary("(^|[ _-])mix([ _-]|$)|hybrid|splatter|poor-quality") |
    mapping$Source_Label == "Micro_out",
  "Mixed or unknown", "Mixed or low-quality",
  "excluded_mixed_or_low_quality", "source mixed or low-quality label"
)

# Dataset-specific abbreviations are resolved before general text rules.
assign_mapping(
  mapping$Dataset == "EGAS00001006537" & has("-exn-"),
  "Excitatory neuronal", "Excitatory neuron",
  "mapped_for_common_comparison", "Cameron ExN abbreviation"
)
assign_mapping(
  mapping$Dataset == "EGAS00001006537" & has("-inn-"),
  "Inhibitory neuronal", "Inhibitory neuron",
  "mapped_for_common_comparison", "Cameron InN abbreviation"
)
assign_mapping(
  mapping$Dataset == "EGAS00001006537" & has("-rg-"),
  "Neural progenitor", "Radial glia",
  "mapped_for_common_comparison", "Cameron RG abbreviation"
)
assign_mapping(
  mapping$Dataset == "EGAS00001006537" & has("cycpro"),
  "Neural progenitor", "Cycling progenitor",
  "mapped_for_common_comparison", "Cameron CycPro abbreviation"
)
assign_mapping(
  mapping$Dataset == "EGAS00001006537" & has("-ip($| )"),
  "Neural progenitor", "Intermediate progenitor",
  "mapped_for_common_comparison", "Cameron IP abbreviation"
)
assign_mapping(
  mapping$Dataset == "EGAS00001006537" & has("-cr-"),
  "Excitatory neuronal", "Cajal-Retzius cell",
  "mapped_for_common_comparison", "Cameron CR abbreviation"
)
assign_mapping(
  mapping$Dataset == "EGAS00001006537" & has("-opc"),
  "Oligodendrocyte lineage", "Oligodendrocyte precursor cell",
  "mapped_for_common_comparison", "Cameron OPC abbreviation"
)
assign_mapping(
  mapping$Dataset == "EGAS00001006537" & has("-mg"),
  "Myeloid or immune", "Microglia",
  "mapped_for_common_comparison", "Cameron MG abbreviation"
)
assign_mapping(
  mapping$Dataset == "EGAS00001006537" & has("-endo"),
  "Vascular or mesenchymal", "Endothelial cell",
  "mapped_for_common_comparison", "Cameron Endo abbreviation"
)

assign_mapping(
  mapping$Dataset == "GSE217511" & has("cpn|spn"),
  "Excitatory neuronal", "Excitatory neuron",
  "mapped_for_common_comparison", "source cortical or subplate projection neuron"
)
assign_mapping(
  mapping$Dataset == "GSE217511" &
    mapping$Source_Label %chin% c("IN", "IN 5HTR3a", "IN SOM", "IN PV", "MSN"),
  "Inhibitory neuronal", "Inhibitory neuron",
  "mapped_for_common_comparison", "source interneuron or medium spiny neuron"
)
assign_mapping(
  mapping$Dataset == "GSE217511" & mapping$Source_Label == "CRN",
  "Excitatory neuronal", "Cajal-Retzius cell",
  "mapped_for_common_comparison", "source Cajal-Retzius neuron"
)
assign_mapping(
  mapping$Dataset == "GSE217511" &
    mapping$Source_Label %chin% c("nIPC", "gIPC", "TAC"),
  "Neural progenitor", "Intermediate progenitor",
  "mapped_for_common_comparison", "source IPC or transit-amplifying cell"
)
assign_mapping(
  mapping$Dataset == "GSE217511" &
    mapping$Source_Label %chin% c("AC", "AC-f", "AC-p"),
  "Astroglial", "Astrocyte",
  "mapped_for_common_comparison", "source astrocyte abbreviation"
)
assign_mapping(
  mapping$Dataset == "GSE217511" & mapping$Source_Label == "OL",
  "Oligodendrocyte lineage", "Oligodendrocyte",
  "mapped_for_common_comparison", "source oligodendrocyte abbreviation"
)
assign_mapping(
  mapping$Dataset == "GSE217511" & mapping$Source_Label == "preOL",
  "Oligodendrocyte lineage", "Committed oligodendrocyte precursor",
  "mapped_for_common_comparison", "source pre-oligodendrocyte abbreviation"
)
assign_mapping(
  mapping$Dataset == "GSE217511" & mapping$Source_Label == "MG",
  "Myeloid or immune", "Microglia",
  "mapped_for_common_comparison", "source microglia abbreviation"
)
assign_mapping(
  mapping$Dataset == "GSE217511" & mapping$Source_Label == "BVC",
  "Vascular or mesenchymal", "Vascular or mesenchymal cell",
  "mapped_for_common_comparison", "source blood-vessel-cell abbreviation"
)

assign_mapping(
  mapping$Dataset == "GSE296073" & mapping$Source_Label == "CEN",
  "Excitatory neuronal", "Excitatory neuron",
  "mapped_for_common_comparison", "source cortical excitatory neuron abbreviation"
)
assign_mapping(
  mapping$Dataset == "GSE296073" & mapping$Source_Label == "CIN",
  "Inhibitory neuronal", "Inhibitory neuron",
  "mapped_for_common_comparison", "source cortical inhibitory neuron abbreviation"
)
assign_mapping(
  mapping$Dataset == "GSE296073" & mapping$Source_Label == "GP",
  "Neural progenitor", "Glial progenitor",
  "mapped_for_common_comparison",
  paste(
    "Yu et al. 2025 Extended Data Fig. 2 defines GP as glia progenitors;",
    "DOI 10.1038/s41586-025-09362-8"
  )
)
assign_mapping(
  mapping$Dataset == "GSE296073" & mapping$Source_Label == "SN",
  "Other neuronal", "Subpallial neuron",
  "mapped_for_common_comparison",
  paste(
    "Yu et al. 2025 Extended Data Fig. 2 defines SN as subpallium neurons;",
    "DOI 10.1038/s41586-025-09362-8"
  )
)

assign_mapping(
  mapping$Dataset == "GSE168408" & mapping$Source_Level_1 == "PN",
  "Excitatory neuronal", "Excitatory neuron",
  "mapped_for_common_comparison", "Herring PN source class"
)
assign_mapping(
  mapping$Dataset == "GSE168408" & mapping$Source_Level_1 == "IN",
  "Inhibitory neuronal", "Inhibitory neuron",
  "mapped_for_common_comparison", "Herring IN source class"
)
assign_mapping(
  mapping$Dataset == "GSE204683" & has("(^| )in-(mge|cge|fetal)"),
  "Inhibitory neuronal", "Inhibitory neuron",
  "mapped_for_common_comparison", "source interneuron lineage"
)
assign_mapping(
  mapping$Dataset == "GSE207334" & mapping$Source_Label == "PC",
  "Vascular or mesenchymal", "Pericyte",
  "mapped_for_common_comparison", "source pericyte abbreviation"
)
assign_mapping(
  mapping$Dataset == "GSE207334" & mapping$Source_Label == "L6b",
  "Excitatory neuronal", "Excitatory neuron",
  "mapped_for_common_comparison",
  paste(
    "Ma et al. 2022 classifies L6b as a glutamatergic subclass;",
    "DOI 10.1126/science.abo7257"
  )
)
assign_mapping(
  mapping$Dataset == "GSE207334" & mapping$Source_Label == "TH",
  "Inhibitory neuronal", "Inhibitory neuron",
  "mapped_for_common_comparison",
  paste(
    "Ma et al. 2022 describes human TH-expressing inhibitory neurons;",
    "DOI 10.1126/science.abo7257"
  )
)
assign_mapping(
  mapping$Dataset == "GSE294786" & has("(^| )in-(sst|vip|pv|sv2c)"),
  "Inhibitory neuronal", "Inhibitory neuron",
  "mapped_for_common_comparison", "source inhibitory-neuron subtype"
)
assign_mapping(
  mapping$Dataset == "GSE294786" & has("(^| )l[2-6](/l?[2-6])?(-bs|-cc)?($| )"),
  "Excitatory neuronal", "Excitatory neuron",
  "mapped_for_common_comparison", "source cortical-layer excitatory subtype"
)
assign_mapping(
  mapping$Dataset == "GSE97942" & has("(^| )ex[0-9]"),
  "Excitatory neuronal", "Excitatory neuron",
  "mapped_for_common_comparison", "source Ex cluster abbreviation"
)
assign_mapping(
  mapping$Dataset == "GSE97942" & has("(^| )in[0-9]"),
  "Inhibitory neuronal", "Inhibitory neuron",
  "mapped_for_common_comparison", "source In cluster abbreviation"
)
assign_mapping(
  mapping$Dataset == "GSE97942" & mapping$Source_Label == "Per",
  "Vascular or mesenchymal", "Pericyte",
  "mapped_for_common_comparison", "source pericyte abbreviation"
)
assign_mapping(
  mapping$Dataset == "PRJCA015229" & mapping$Source_Level_1 == "IN",
  "Inhibitory neuronal", "Inhibitory neuron",
  "mapped_for_common_comparison", "source IN class"
)
assign_mapping(
  mapping$Dataset == "Li_et_al_2018" & has("(^| )ipc[0-9]?($| )"),
  "Neural progenitor", "Intermediate progenitor",
  "mapped_for_common_comparison", "source IPC abbreviation"
)
assign_mapping(
  mapping$Dataset == "Li_et_al_2018" & has("(^| )nasn[0-9]?($| )"),
  "Other neuronal", "Immature neuron",
  "mapped_for_common_comparison", "source nascent-neuron abbreviation"
)
assign_mapping(
  mapping$Dataset == "GSE186538" & has("(^| )dg mc "),
  "Excitatory neuronal", "Excitatory neuron",
  "mapped_for_common_comparison", "source dentate-gyrus mossy-cell identity"
)
assign_mapping(
  mapping$Dataset == "GSE186538" & has("(^| )t skap1 cd247"),
  "Myeloid or immune", "Other immune cell",
  "mapped_for_common_comparison", "source T-cell identity"
)
assign_mapping(
  mapping$Dataset == "GSE81475" & mapping$Source_Label == "Olig",
  "Oligodendrocyte lineage", "Oligodendrocyte",
  "mapped_for_common_comparison", "source oligodendrocyte abbreviation"
)
assign_mapping(
  mapping$Dataset == "Ma_et_al_2022" & has("(^| )rb hba1 hbb"),
  "Myeloid or immune", "Other hematopoietic cell",
  "mapped_for_common_comparison", "source erythroid-cell marker label"
)
assign_mapping(
  mapping$Dataset == "HYPOMAP" & has("(^| )gaba([ -]|$)"),
  "Inhibitory neuronal", "Inhibitory neuron",
  "mapped_for_common_comparison", "HYPOMAP author GABA hierarchy"
)
assign_mapping(
  mapping$Dataset == "HYPOMAP" & has("(^| )glu([ -]|$)"),
  "Excitatory neuronal", "Excitatory neuron",
  "mapped_for_common_comparison", "HYPOMAP author GLU hierarchy"
)
assign_mapping(
  mapping$Dataset == "HYPOMAP" & has("(^| )chol([ -]|$)"),
  "Other neuronal", "Cholinergic neuron",
  "mapped_for_common_comparison", "HYPOMAP author CHOL hierarchy"
)
assign_mapping(
  mapping$Dataset == "HYPOMAP" & has("(^| )hdc([ -]|$)"),
  "Other neuronal", "Histaminergic neuron",
  "mapped_for_common_comparison", "HYPOMAP author HDC hierarchy"
)
assign_mapping(
  mapping$Dataset == "HYPOMAP" & has("p2rx2[_ -]?otp"),
  "Other neuronal", "P2RX2-OTP neuron",
  "mapped_for_common_comparison", "HYPOMAP author P2RX2-OTP identity"
)
assign_mapping(
  mapping$Dataset == "Velmeshev_2023" & mapping$Source_Label == "VASC",
  "Vascular or mesenchymal", "Vascular or mesenchymal cell",
  "mapped_for_common_comparison", "Velmeshev VASC lineage"
)

assign_mapping(
  mapping$Dataset == "Velmeshev_2023" & mapping$Source_Label == "ExNeu",
  "Excitatory neuronal", "Excitatory neuron",
  "mapped_for_common_comparison", "Velmeshev ExNeu lineage"
)
assign_mapping(
  mapping$Dataset == "Velmeshev_2023" & mapping$Source_Label == "IN",
  "Inhibitory neuronal", "Inhibitory neuron",
  "mapped_for_common_comparison", "Velmeshev IN lineage"
)
assign_mapping(
  mapping$Dataset == "Velmeshev_2023" & mapping$Source_Label == "GLIALPROG",
  "Neural progenitor", "Glial progenitor",
  "mapped_for_common_comparison", "Velmeshev GLIALPROG lineage"
)

assign_mapping(
  mapping$Dataset == "GSE97942" & has("(^| )purk[12]?($| )"),
  "Inhibitory neuronal", "Purkinje neuron",
  "mapped_for_common_comparison", "source Purkinje label"
)
assign_mapping(
  mapping$Dataset == "GSE97942" & has("(^| )gran($| )"),
  "Excitatory neuronal", "Granule neuron",
  "mapped_for_common_comparison", "source granule-neuron label"
)

# Explicit non-neuronal identities must be resolved before interpreting layer
# tokens such as L2-6 or the generic word 'progenitor'. This prevents labels
# such as 'Oligo L2-6' and 'oligodendrocyte progenitor' from being assigned to
# neuronal or generic neural-progenitor groups.
assign_mapping(
  has("committed oligodendrocyte precursor|(^|[ _-])cop($|[ _-])|preol|gpr17[.]oligodendrocytes"),
  "Oligodendrocyte lineage", "Committed oligodendrocyte precursor",
  "mapped_for_common_comparison", "committed oligodendrocyte precursor identity"
)
assign_mapping(
  has("oligodendrocyte (precursor|progenitor)|(^|[ _-])opc(s)?([ _-]|$)"),
  "Oligodendrocyte lineage", "Oligodendrocyte precursor cell",
  "mapped_for_common_comparison", "OPC source identity"
)
assign_mapping(
  has("oligodendrocyte|(^|[ _-])oligo(s)?([ _-]|$)|(^|[ _-])oli($|[ _-])|(^|[ _-])odc($|[ _-])|(^|[ _-])ol($|[ _-])"),
  "Oligodendrocyte lineage", "Oligodendrocyte",
  "mapped_for_common_comparison", "oligodendrocyte source identity"
)
assign_mapping(
  has("astro|(^|[ _-])ast($|[ _-])|(^|[ _-])ac-f($|[ _-])|(^|[ _-])ac-p($|[ _-])"),
  "Astroglial", "Astrocyte",
  "mapped_for_common_comparison", "astrocyte source identity"
)
assign_mapping(
  has("microglia|(^|[ _-])micro($|[ _-])|(^|[ _-])mic($|[ _-])|(^|[ _-])mg($|[ _-])|tmem119[.]immune|p2ry12|micro/macro"),
  "Myeloid or immune", "Microglia",
  "mapped_for_common_comparison", "microglia source identity"
)
assign_mapping(
  has("(^|[ _-])macro([ _-]|$)"),
  "Myeloid or immune", "Border-associated macrophage",
  "mapped_for_common_comparison", "macrophage source identity"
)
assign_mapping(
  has("(^|[ _-])myeloid([ _-]|$)"),
  "Myeloid or immune", "Other immune cell",
  "mapped_for_common_comparison", "myeloid source identity"
)
assign_mapping(
  has("endothelial|(^|[ _-])endo($|[ _-])|(^|[ _-])end($|[ _-])"),
  "Vascular or mesenchymal", "Endothelial cell",
  "mapped_for_common_comparison", "endothelial source identity"
)
assign_mapping(
  has("pericyte|(^|[ _-])perc($|[ _-])|(^|[ _-])pc p2ry14|(^|[ _-])pc cldn5"),
  "Vascular or mesenchymal", "Pericyte",
  "mapped_for_common_comparison", "pericyte source identity"
)
assign_mapping(
  has("vascular smooth|(^|[ _-])vsmc($|[ _-])|(^|[ _-])asmc($|[ _-])|(^|[ _-])smc($|[ _-])"),
  "Vascular or mesenchymal", "Vascular smooth muscle cell",
  "mapped_for_common_comparison", "vascular smooth-muscle source identity"
)
assign_mapping(
  has("vlmc|fibroblast|mural|(^|[ _-])vasc?($|[ _-])|vascular|blood-vessel|(^|[ _-])bvc($|[ _-])"),
  "Vascular or mesenchymal", "Vascular or mesenchymal cell",
  "mapped_for_common_comparison", "vascular or mesenchymal source identity"
)

# Specific developmental and specialised identities.
assign_mapping(
  has("cajal|(^|[ _-])crn($|[ _-])|(^|[ _-])cr[12]?($|[ _-])"),
  "Excitatory neuronal", "Cajal-Retzius cell",
  "mapped_for_common_comparison", "Cajal-Retzius source identity"
)
assign_mapping(
  has("purkinje|(^|[ _-])purk"),
  "Inhibitory neuronal", "Purkinje neuron",
  "mapped_for_common_comparison", "Purkinje source identity"
)
assign_mapping(
  has("granule|(^|[ _-])gran($|[ _-])|dg gc"),
  "Excitatory neuronal", "Granule neuron",
  "mapped_for_common_comparison", "granule-neuron source identity"
)
assign_mapping(
  has("cycling progenitor|cycpro"),
  "Neural progenitor", "Cycling progenitor",
  "mapped_for_common_comparison", "cycling progenitor source identity"
)
assign_mapping(
  has("intermediate progenitor|(^|[ _-])[gn]?ipc([ _-]|$)|transit"),
  "Neural progenitor", "Intermediate progenitor",
  "mapped_for_common_comparison", "intermediate progenitor source identity"
)
assign_mapping(
  has("radial glia|(^|[ _-])rg([ _-]|$)|neprgc"),
  "Neural progenitor", "Radial glia",
  "mapped_for_common_comparison", "radial-glia source identity"
)
assign_mapping(
  has("glialprog|glial progenitor"),
  "Neural progenitor", "Glial progenitor",
  "mapped_for_common_comparison", "glial-progenitor source identity"
)
assign_mapping(
  has("neural progenitor|neuroepithelial|(^|[ _-])npc([ _-]|$)|stem cells|fetal_quiescent|fetal_replicating|(^|[ _-])progenitor($|[ _-])|upper rhombic lip|lower rhombic lip"),
  "Neural progenitor", "Neural progenitor",
  "mapped_for_common_comparison", "neural-progenitor source identity"
)

# Neuronal identities. These rules use biological terms or unambiguous source
# abbreviations; generic 'Neuron' remains an explicit other-neuron category.
# When an original study supplies an explicit broad hierarchy, that hierarchy
# takes precedence over subtype tokens embedded in a finer label. This avoids,
# for example, treating an Allen glutamatergic "Exc ... LAMP5" type as an
# inhibitory neuron merely because its fine label contains LAMP5.
source_level_1 <- normalise_text(mapping$Source_Level_1)
assign_mapping(
  source_level_1 %chin% c(
    "glutamatergic", "excitatory", "excitatory neuronal"
  ),
  "Excitatory neuronal", "Excitatory neuron",
  "mapped_for_common_comparison",
  "explicit source level-1 excitatory hierarchy"
)
assign_mapping(
  source_level_1 %chin% c(
    "gabaergic", "inhibitory", "inhibitory neuronal"
  ),
  "Inhibitory neuronal", "Inhibitory neuron",
  "mapped_for_common_comparison",
  "explicit source level-1 inhibitory hierarchy"
)
assign_mapping(
  has("gabaergic|inhibitory|interneuron|(^|[ _-])inh([ _-]|$)|(^|[ _-])inn[0-9]?([ _-]|$)|(^|[ _-])inhib([ _-]|$)|(^|[ _-])cin($|[ _-])|pvalb|(^|[ _-])pv($|[ _-])|(^|[ _-])sst($|[ _-])|vip|lamp5|sncg|adarb2|medium spiny|(^|[ _-])msn($|[ _-])"),
  "Inhibitory neuronal", "Inhibitory neuron",
  "mapped_for_common_comparison", "inhibitory-neuron source identity"
)
assign_mapping(
  has("glutamatergic|excitatory|(^|[ _-])exc([ _-]|$)|(^|[ _-])exn?[0-9]?([ _-]|$)|(^|[ _-])en-fetal|(^|[ _-])en($|[ _-])|(^|[ _-])cen($|[ _-])|projection neuron|intratelencephalic|corticothalamic|near-projecting|(^|[ _-])cpn($|[ _-])|(^|[ _-])spn($|[ _-])|(^|[ _-])l[2-6]/?[2-6]?($|[ _-])|(^|[ _-])ca[1-4]($|[ _-])|(^|[ _-])sub($|[ _-])|amygdala excitatory|thalamic excitatory"),
  "Excitatory neuronal", "Excitatory neuron",
  "mapped_for_common_comparison", "excitatory-neuron source identity"
)
assign_mapping(
  has("(^|[ _-])neuron(s)?($|[ _-])|(^|[ _-])sn($|[ _-])"),
  "Other neuronal", "Other neuron",
  "mapped_for_common_comparison", "source reports neuron without common subtype"
)

# Macroglial and oligodendrocyte-lineage identities.
assign_mapping(
  has("bergmann"),
  "Astroglial", "Bergmann glia",
  "mapped_for_common_comparison", "Bergmann-glia source identity"
)
assign_mapping(
  has("astro|(^|[ _-])ast($|[ _-])|(^|[ _-])ac-f($|[ _-])|(^|[ _-])ac-p($|[ _-])"),
  "Astroglial", "Astrocyte",
  "mapped_for_common_comparison", "astrocyte source identity"
)
assign_mapping(
  has("committed oligodendrocyte precursor|(^|[ _-])cop($|[ _-])|preol|gpr17[.]oligodendrocytes"),
  "Oligodendrocyte lineage", "Committed oligodendrocyte precursor",
  "mapped_for_common_comparison", "committed oligodendrocyte precursor identity"
)
assign_mapping(
  has("oligodendrocyte precursor|(^|[ _-])opc(s)?([ _-]|$)"),
  "Oligodendrocyte lineage", "Oligodendrocyte precursor cell",
  "mapped_for_common_comparison", "OPC source identity"
)
assign_mapping(
  has("oligodendrocyte|(^|[ _-])oligo(s)?([ _-]|$)|(^|[ _-])oli($|[ _-])|(^|[ _-])odc($|[ _-])|(^|[ _-])ol($|[ _-])"),
  "Oligodendrocyte lineage", "Oligodendrocyte",
  "mapped_for_common_comparison", "oligodendrocyte source identity"
)

# Myeloid, vascular, mesenchymal and epithelial identities.
assign_mapping(
  has("microglia|(^|[ _-])micro($|[ _-])|(^|[ _-])mic($|[ _-])|(^|[ _-])mg($|[ _-])|tmem119[.]immune|p2ry12|micro/macro"),
  "Myeloid or immune", "Microglia",
  "mapped_for_common_comparison", "microglia source identity"
)
assign_mapping(
  has("immune|macrophage|f13a1[.]immune|ms4a4b[.]immune"),
  "Myeloid or immune", "Other immune cell",
  "mapped_for_common_comparison", "other immune source identity"
)
assign_mapping(
  has("endothelial|(^|[ _-])endo($|[ _-])|(^|[ _-])end($|[ _-])"),
  "Vascular or mesenchymal", "Endothelial cell",
  "mapped_for_common_comparison", "endothelial source identity"
)
assign_mapping(
  has("pericyte|(^|[ _-])perc($|[ _-])|(^|[ _-])pc p2ry14|(^|[ _-])pc cldn5"),
  "Vascular or mesenchymal", "Pericyte",
  "mapped_for_common_comparison", "pericyte source identity"
)
assign_mapping(
  has("vascular smooth|(^|[ _-])vsmc($|[ _-])|(^|[ _-])smc($|[ _-])"),
  "Vascular or mesenchymal", "Vascular smooth muscle cell",
  "mapped_for_common_comparison", "vascular smooth-muscle source identity"
)
assign_mapping(
  has("vlmc|fibroblast|mural|(^|[ _-])vasc?($|[ _-])|vascular|blood-vessel|(^|[ _-])bvc($|[ _-])"),
  "Vascular or mesenchymal", "Vascular or mesenchymal cell",
  "mapped_for_common_comparison", "vascular or mesenchymal source identity"
)
assign_mapping(
  has("ependymal"),
  "Epithelial or barrier", "Ependymal cell",
  "mapped_for_common_comparison", "ependymal source identity"
)
assign_mapping(
  has("tanycyte"),
  "Epithelial or barrier", "Tanycyte",
  "mapped_for_common_comparison", "tanycyte source identity"
)
assign_mapping(
  has("choroid plexus"),
  "Epithelial or barrier", "Choroid plexus cell",
  "mapped_for_common_comparison", "choroid-plexus source identity"
)
assign_mapping(
  has("epith"),
  "Epithelial or barrier", "Other epithelial cell",
  "mapped_for_common_comparison", "epithelial source identity"
)

# Anything not supported by an explicit source-semantic rule remains visible.
assign_mapping(
  rep(TRUE, nrow(mapping)),
  NA_character_, NA_character_,
  "manual_review_required", "no reviewed common-label rule"
)

if (anyNA(mapping$Mapping_Status)) {
  stop("Source-label harmonization left rows without a status")
}
# These source labels do not distinguish the common neuronal identities.
# Retain them in the resource and correspondence table, but not in cLISI.
mapping[Mapping_Status == "mapped_for_common_comparison" &
  Harmonized_CellType %chin% c("Other neuron", "Immature neuron", "Subpallial neuron", "P2RX2-OTP neuron"),
  `:=`(Mapping_Status = "excluded_incomparable_source_resolution",
    Mapping_Rule = "source neuronal identity lacks comparable label resolution")]
mapping[, c("Search_Text", "Primary_Text") := NULL]
setorder(mapping, Dataset, -Cells, Source_Label)

summary <- mapping[, .(
  Source_Hierarchy_Rows = .N,
  Cells = sum(Cells),
  Cell_Fraction = sum(Cells) / sum(mapping$Cells),
  Datasets = uniqueN(Dataset)
), by = .(
  Mapping_Status, Harmonized_CellClass, Harmonized_CellType
)]
setorder(summary, Mapping_Status, Harmonized_CellClass, -Cells)

write_tsv(as.data.frame(mapping), output_file)
write_tsv(as.data.frame(summary), summary_file)

metadata_file <- file.path(annotation_dir, "source_label_evaluation_metadata.rds")
if (!file.exists(metadata_file)) stop("Run source_label_audit.R to refresh author labels first")
if (file.exists(metadata_file)) {
  thisutils::log_message("[source-label-harmonization] ", "Applying the reviewed hierarchy mapping to complete metadata")
  metadata <- readRDS(metadata_file)
  if (nrow(metadata) == 0L || anyDuplicated(metadata$Cells)) {
    stop("Complete metadata cells are empty or duplicated")
  }
  if (!identical(
    sort(unique(as.character(metadata$Dataset))),
    sort(reference_datasets())
  )) {
    stop("Complete metadata does not contain the reference datasets")
  }
  join_fields <- c(
    "Dataset", "Source_Original_Label", "Source_Label", "Source_Level_1",
    "Source_Level_2", "Source_Level_3", "Source_Cluster_ID"
  )
  metadata_source_fields <- c(
    Source_Original_Label = "source_cell_type_original_label",
    Source_Label = "source_cell_type_label",
    Source_Level_1 = "source_cell_type_level_1",
    Source_Level_2 = "source_cell_type_level_2",
    Source_Level_3 = "source_cell_type_level_3",
    Source_Cluster_ID = "source_cluster_id"
  )
  missing_metadata_fields <- setdiff(
    unname(metadata_source_fields),
    names(metadata)
  )
  if (length(missing_metadata_fields) > 0L) {
    stop(
      "Complete metadata is missing source hierarchy fields: ",
      paste(missing_metadata_fields, collapse = ", ")
    )
  }
  cell_source <- data.table(
    Row_Order = seq_len(nrow(metadata)),
    Cells = as.character(metadata$Cells),
    Dataset = as.character(metadata$Dataset)
  )
  for (output_name in names(metadata_source_fields)) {
    source_name <- metadata_source_fields[[output_name]]
    cell_source[[output_name]] <- normalize_missing_metadata(
      metadata[[source_name]]
    )
  }
  lookup <- mapping[, c(
    join_fields,
    "Harmonized_CellClass", "Harmonized_CellType", "Mapping_Status",
    "Mapping_Rule"
  ), with = FALSE]
  if (anyDuplicated(lookup[, ..join_fields])) {
    stop("Reviewed hierarchy mapping contains duplicated join keys")
  }
  cell_mapping <- lookup[cell_source, on = join_fields]
  setorder(cell_mapping, Row_Order)
  if (!identical(cell_mapping$Cells, as.character(metadata$Cells))) {
    stop("Per-cell source-label join changed cell order")
  }
  source_available <- rowSums(!is.na(as.data.frame(
    cell_mapping[, c(
      "Source_Original_Label", "Source_Label", "Source_Level_1",
      "Source_Level_2", "Source_Level_3", "Source_Cluster_ID"
    ), with = FALSE]
  ))) > 0L
  cell_mapping[
    is.na(Mapping_Status) & !source_available,
    `:=`(
      Mapping_Status = "source_label_unavailable",
      Mapping_Rule = "original study did not release a matching cell label"
    )
  ]
  if (any(is.na(cell_mapping$Mapping_Status))) {
    stop("Source-annotated cells failed the exact reviewed hierarchy join")
  }
  cell_mapping[, Source_Label_Evaluation_Eligible :=
    Mapping_Status == "mapped_for_common_comparison"]
  donor_field <- if ("Global_Donor_ID" %in% names(metadata)) {
    "Global_Donor_ID"
  } else {
    "Donor_ID"
  }
  cell_mapping[, Donor_ID := normalize_missing_metadata(
    metadata[[donor_field]]
  )]
  setcolorder(cell_mapping, c(
    "Cells", "Dataset", "Donor_ID", "Source_Original_Label",
    "Source_Label", "Source_Level_1", "Source_Level_2", "Source_Level_3",
    "Source_Cluster_ID", "Harmonized_CellClass", "Harmonized_CellType",
    "Mapping_Status", "Mapping_Rule", "Source_Label_Evaluation_Eligible",
    "Row_Order"
  ))
  cell_output_file <- file.path(
    annotation_dir,
    "source_labels_harmonized.rds"
  )
  temporary_cell_output <- paste0(cell_output_file, ".tmp.", Sys.getpid())
  saveRDS(as.data.frame(cell_mapping), temporary_cell_output, compress = TRUE)
  if (!file.rename(temporary_cell_output, cell_output_file)) {
    unlink(temporary_cell_output)
    stop("Could not publish per-cell source-label mapping")
  }
  cell_summary <- cell_mapping[, .(
    Cells = .N,
    Donors = uniqueN(Donor_ID[!is.na(Donor_ID)]),
    Cell_Fraction = .N / nrow(cell_mapping)
  ), by = .(
    Dataset, Mapping_Status, Harmonized_CellClass, Harmonized_CellType
  )]
  setorder(cell_summary, Dataset, Mapping_Status, -Cells)
  write_tsv(
    as.data.frame(cell_summary),
    file.path(annotation_dir, "source_label_cell_summary.tsv")
  )
  cell_contract <- data.frame(
    Field = c(
      "Cells", "Datasets", "Mapped_Evaluation_Cells",
      "Source_Unknown_Cells", "Mixed_Low_Quality_Cells",
      "Source_Label_Unavailable_Cells", "Incomparable_Source_Resolution_Cells", "Output_File"
    ),
    Value = c(
      nrow(cell_mapping),
      uniqueN(cell_mapping$Dataset),
      sum(cell_mapping$Mapping_Status == "mapped_for_common_comparison"),
      sum(cell_mapping$Mapping_Status == "excluded_source_unknown"),
      sum(cell_mapping$Mapping_Status == "excluded_mixed_or_low_quality"),
      sum(cell_mapping$Mapping_Status == "source_label_unavailable"),
      sum(cell_mapping$Mapping_Status == "excluded_incomparable_source_resolution"),
      cell_output_file
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(
    cell_contract,
    file.path(annotation_dir, "source_label_cell_contract.tsv")
  )
}

thisutils::log_message("[source-label-harmonization] ", paste(
  "Harmonized",
  nrow(mapping),
  "source hierarchy rows covering",
  format(sum(mapping$Cells), big.mark = ","),
  "source-annotated cells;",
  format(sum(mapping[Mapping_Status == "mapped_for_common_comparison"]$Cells),
    big.mark = ","
  ),
  "cells are eligible for common-label evaluation"
))
