#!/bin/bash

set -e

source "functions/utils.sh"

code_dir="processing"
overwrite="${1:-F}"
check_command Rscript
check_command python3

should_process() {
  local target_file="$1"
  if [[ "$overwrite" =~ ^([Tt]|[Tt][Rr][Uu][Ee]|1)$ ]]; then
    return 0
  fi
  if [ ! -f "$target_file" ]; then
    return 0
  fi
  return 1
}


# BCAtlas
# title: A brain cell atlas integrating single-cell transcriptomes across human brain regions
# paper: https://doi.org/10.1038/s41591-024-03150-z
# data: https://www.braincellatlas.org/dataSet
log_message "Processing BCAtlas data..."
run_r_script "$code_dir" "BCAtlas.R" "BCAtlas"

# BICCN
# data: https://brainscope.gersteinlab.org/integrative_files.html

# GSE103723
# paper: https://doi.org/10.1126/sciadv.adg3754
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/37824614
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103723
# code:
if should_process "../../data/BrainOmicsData/processed/GSE103723/GSE103723_processed.rds"; then
  log_message "Processing GSE103723 data..."
  run_r_script "$code_dir" "GSE103723.R" "GSE103723"
else
  log_message "GSE103723 data already processed!"
fi


# GSE104276
# paper: https://doi.org/10.1126/sciadv.adg3754
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/37824614
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104276
# code:
if should_process "../../data/BrainOmicsData/processed/GSE104276/GSE104276_processed.rds"; then
  log_message "Processing GSE104276 data..."
  run_r_script "$code_dir" "GSE104276.R" "GSE104276"
else
  log_message "GSE104276 data already processed!"
fi


# # GSE126836 (no age information)
# # paper: https://doi.org/10.1016/j.cell.2019.05.006
# # pmid: https://pubmed.ncbi.nlm.nih.gov/31178122/
# # data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126836
# # code:
# if should_process "../../data/BrainOmicsData/processed/GSE126836/GSE126836_processed.rds"; then
#   log_message "Processing GSE126836 data..."
#   Rscript $code_dir/GSE126836.R
#   log_message "GSE126836 data processed successfully!" --message-type success
# else
#   log_message "GSE126836 data already processed!"
# fi


# GSE186538
# paper: https://doi.org/10.1126/sciadv.adg3754
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/37824614
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186538
# code:
if should_process "../../data/BrainOmicsData/processed/GSE186538/GSE186538_processed.rds"; then
  log_message "Processing GSE186538 data..."
  run_r_script "$code_dir" "GSE186538.R" "GSE186538"
else
  log_message "GSE186538 data already processed!"
fi


# GSE199762
# paper: https://doi.org/10.1038/s41586-023-06981-x
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/38122823
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199762
# dbGaP, accession: https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs003509.v1.p1
# code: https://github.com/massisnascimento/ECstream
if should_process "../../data/BrainOmicsData/processed/GSE199762/GSE199762_processed.rds"; then
  log_message "Processing GSE199762 data..."
  run_r_script "$code_dir" "GSE199762.R" "GSE199762"
else
  log_message "GSE199762 data already processed!"
fi


# GSE204683 (multiome: snRNA-seq + snATAC-seq (GSE204684))
# paper: https://doi.org/10.1126/sciadv.adg3754
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/37824614
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE204683
# ATAC-seq: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE204684
# CELLxGENE (RRID: SCR_021059) data (h5ad):
# https://cellxgene.cziscience.com/collections/ceb895f4-ff9f-403a-b7c3-187a9657ac2c
# code: https://doi.org/10.5281/zenodo.7703253
if should_process "../../data/BrainOmicsData/processed/GSE204683/GSE204683_processed.rds"; then
  log_message "Processing GSE204683 data..."
  run_r_script "$code_dir" "GSE204683.R" "GSE204683"
else
  log_message "GSE204683 data already processed!"
fi


# GSE212606
# paper: https://doi.org/10.1126/sciadv.adg3754
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/37824614
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212606
# code:
if should_process "../../data/BrainOmicsData/processed/GSE212606/GSE212606_processed.rds"; then
  log_message "Processing GSE212606 data..."
  run_r_script "$code_dir" "GSE212606.R" "GSE212606"
else
  log_message "GSE212606 data already processed!"
fi


# GSE217511
# paper: https://doi.org/10.1038/s41467-022-34975-2
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/36509746
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217511
if should_process "../../data/BrainOmicsData/processed/GSE217511/GSE217511_processed.rds"; then
  log_message "Processing GSE217511 data..."
  run_r_script "$code_dir" "GSE217511.R" "GSE217511"
else
  log_message "GSE217511 data already processed!"
fi


# GSE67835
# paper: https://doi.org/10.1073/pnas.1507125112
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/26060301
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67835
if should_process "../../data/BrainOmicsData/processed/GSE67835/GSE67835_processed.rds"; then
  log_message "Processing GSE67835 data..."
  run_r_script "$code_dir" "GSE67835.R" "GSE67835"
else
  log_message "GSE67835 data already processed!"
fi


# GSE81475
# paper: https://doi.org/10.1016/j.celrep.2016.08.038
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/27568284
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81475
# code:
if should_process "../../data/BrainOmicsData/processed/GSE81475/GSE81475_processed.rds"; then
  log_message "Processing GSE81475 data..."
  run_r_script "$code_dir" "GSE81475.R" "GSE81475"
else
  log_message "GSE81475 data already processed!"
fi


# GSE97942 (contains GSE97887 + GSE97930)
# journal: Nature Biotechnology
# date: 2018
# paper: https://doi.org/10.1038/nbt.4038
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/29227469
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97942
# GSE97887 (scTHS-seq): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97887
# GSE97930 (snDrop-seq): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97930
# code:
if should_process "../../data/BrainOmicsData/processed/GSE97942/GSE97942_processed.rds"; then
  log_message "Processing GSE97942 data..."
  run_r_script "$code_dir" "GSE97942.R" "GSE97942"
else
  log_message "GSE97942 data already processed!"
fi


# Li et al. 2018
# paper: https://doi.org/10.1126/science.aat7615
# pmid: https://pubmed.ncbi.nlm.nih.gov/30545854/
# data:
# code:
if should_process "../../data/BrainOmicsData/processed/Li_et_al_2018/Li_et_al_2018_processed.rds"; then
  log_message "Processing Li_et_al_2018 data..."
  run_r_script "$code_dir" "Li_et_al_2018.R" "Li_et_al_2018"
else
  log_message "Li_et_al_2018 data already processed!"
fi


# Nowakowski_et_al_2017
# paper: https://doi.org/10.1126/science.aap8809
# pmid: https://pubmed.ncbi.nlm.nih.gov/29217575/
# data:
# code:
if should_process "../../data/BrainOmicsData/processed/Nowakowski_et_al_2017/Nowakowski_et_al_2017_processed.rds"; then
  log_message "Processing Nowakowski_et_al_2017 data..."
  run_r_script "$code_dir" "Nowakowski_et_al_2017.R" "Nowakowski_et_al_2017"
else
  log_message "Nowakowski_et_al_2017 data already processed!"
fi


# BTSatlas
# paper: https://doi.org/10.1038/s12276-024-01328-6
# data: https://zenodo.org/records/10939707
# code:
# Contains multiple datasets:
# "AllenM1", "EGAD00001006049", "EGAS00001006537",
# "GSE178175", "GSE168408",
# "GSE144136", "GSE202210"
# first run BTSatlas-1.py to create compatible h5ad file
# then run BTSatlas-2.R to split the data into multiple datasets
if [ ! -f "../../data/BrainOmicsData/processed/BTSatlas/BTS_atlas_compatible.h5ad" ]; then
  log_message "Processing BTSatlas data..."
  run_python_script "$code_dir" "BTSatlas_1.py" "BTSatlas"
else
  log_message "BTSatlas compatible h5ad file already exists!"
fi

log_message "Processing BTSatlas sub-datasets..."
run_r_script "$code_dir" "BTSatlas_2.R" "BTSatlas sub-datasets"

# AllenM1
# data: https://brain-map.org/our-research/cell-types-taxonomies/cell-types-database-rna-seq-data/human-m1-10x
if should_process "../../data/BrainOmicsData/processed/AllenM1/AllenM1_processed.rds"; then
  log_message "Processing AllenM1 data..."
  run_r_script "$code_dir" "AllenM1.R" "AllenM1"
else
  log_message "AllenM1 data already processed!"
fi


# EGAD00001006049
# paper: Comprehensive cell atlas of the first-trimester developing human brain
# doi: https://doi.org/10.1126/science.adf1226
# pmid: 37824650
# data: https://ega-archive.org/datasets/EGAD00001006049
# code:
if should_process "../../data/BrainOmicsData/processed/EGAD00001006049/EGAD00001006049_processed.rds"; then
  log_message "Processing EGAD00001006049 data..."
  run_r_script "$code_dir" "EGAD00001006049.R" "EGAD00001006049"
else
  log_message "EGAD00001006049 data already processed!"
fi


# EGAS00001006537
# paper: Single-Nuclei RNA Sequencing of 5 Regions of the Human Prenatal Brain Implicates Developing Neuron Populations in Genetic Risk for Schizophrenia
# pmid: 36150908
# data: https://ega-archive.org/studies/EGAS00001006537
# code:
if should_process "../../data/BrainOmicsData/processed/EGAS00001006537/EGAS00001006537_processed.rds"; then
  log_message "Processing EGAS00001006537 data..."
  run_r_script "$code_dir" "EGAS00001006537.R" "EGAS00001006537"
else
  log_message "EGAS00001006537 data already processed!"
fi


# GSE178175
# paper:
# pmid:
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178175
# code:
if should_process "../../data/BrainOmicsData/processed/GSE178175/GSE178175_processed.rds"; then
  log_message "Processing GSE178175 data..."
  run_r_script "$code_dir" "GSE178175.R" "GSE178175"
else
  log_message "GSE178175 data already processed!"
fi


# GSE168408
# paper:
# pmid:
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168408
# code:
if should_process "../../data/BrainOmicsData/processed/GSE168408/GSE168408_processed.rds"; then
  log_message "Processing GSE168408 data..."
  run_r_script "$code_dir" "GSE168408.R" "GSE168408"
else
  log_message "GSE168408 data already processed!"
fi


# GSE144136
# paper:
# pmid:
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144136
# code:
if should_process "../../data/BrainOmicsData/processed/GSE144136/GSE144136_processed.rds"; then
  log_message "Processing GSE144136 data..."
  run_r_script "$code_dir" "GSE144136.R" "GSE144136"
else
  log_message "GSE144136 data already processed!"
fi


# GSE202210
# paper:
# pmid:
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202210
# code:
if should_process "../../data/BrainOmicsData/processed/GSE202210/GSE202210_processed.rds"; then
  log_message "Processing GSE202210 data..."
  run_r_script "$code_dir" "GSE202210.R" "GSE202210"
else
  log_message "GSE202210 data already processed!"
fi


# SCR_016152
# SCR_016152 contains GSE207334 and Ma_et_al_2022,
# which are from the same study: https://www.ncbi.nlm.nih.gov/pubmed/36007006
# paper: https://doi.org/10.1126/science.abo7257
# https://pmc.ncbi.nlm.nih.gov/articles/PMC9614553/
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/36007006
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207334
# http://resources.sestanlab.org/PFC/
# https://brainscope.gersteinlab.org/
# code:
# all samples not include GSE207334
# object_all <- readRDS(
#   "../../data/BrainOmicsData/raw/GSE207334/PFC_snRNAseq_liftover.rds"
# )

# GSE207334 (multiome: snRNA-seq + snATAC-seq)
if should_process "../../data/BrainOmicsData/processed/GSE207334/GSE207334_processed.rds"; then
  log_message "Processing GSE207334 data..."
  run_r_script "$code_dir" "GSE207334.R" "GSE207334"
else
  log_message "GSE207334 data already processed!"
fi

# Ma_et_al_2022
if should_process "../../data/BrainOmicsData/processed/Ma_et_al_2022/Ma_et_al_2022_processed.rds"; then
  log_message "Processing Ma_et_al_2022 data..."
  run_r_script "$code_dir" "Ma_et_al_2022.R" "Ma_et_al_2022"
else
  log_message "Ma_et_al_2022 data already processed!"
fi


# GSE235493 (multiome: snRNA-seq + snATAC-seq, Macaque)
# paper: https://doi.org/10.1016/j.neuron.2025.04.025
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235493
# code: https://doi.org/10.5281/zenodo.15243470
# note: not found human data in GSE235493


# HYPOMAP
# paper: https://doi.org/10.1038/s41586-024-08504-8
# code:
#   https://github.com/lsteuernagel/HYPOMAP
#   https://github.com/georgiedowsett/HYPOMAP
#   https://github.com/lsteuernagel/scIntegration
#   https://github.com/mrcepid-rap
# data: https://cellxgene.cziscience.com/collections/d0941303-7ce3-4422-9249-cf31eb98c480
# data(spatial): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278848
if should_process "../../data/BrainOmicsData/processed/HYPOMAP/HYPOMAP_processed.rds"; then
  log_message "Processing HYPOMAP data..."
  run_r_script "$code_dir" "HYPOMAP.R" "HYPOMAP"
else
  log_message "HYPOMAP data already processed!"
fi


# SomaMut
# papaer: https://doi.org/10.1038/s41586-025-09435-8
# code: ~
# data: https://publications.wenglab.org/SomaMut/
if should_process "../../data/BrainOmicsData/processed/SomaMut/SomaMut_processed.rds"; then
  log_message "Processing SomaMut data..."
  run_r_script "$code_dir" "SomaMut.R" "SomaMut"
else
  log_message "SomaMut data already processed!"
fi


# PRJCA015229 (multiome: snRNA-seq + snATAC-seq, Human + Macaque)
# paper: https://doi.org/10.1016/j.xgen.2024.100703
# code: https://github.com/KIZ-SubLab/ACC-sn-Multiomes
# data: https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA015229
if should_process "../../data/BrainOmicsData/processed/PRJCA015229/PRJCA015229_processed.rds"; then
  log_message "Processing PRJCA015229 data..."
  run_r_script "$code_dir" "PRJCA015229.R" "PRJCA015229"
else
  log_message "PRJCA015229 data already processed!"
fi


# ROSMAP (Religious Order Study (ROS) or the Rush Memory and Aging Project (MAP))
# paper: https://doi.org/10.1016/j.cell.2023.08.039
# code: https://github.com/mathyslab7/ROSMAP_snRNAseq_PFC/
# data: https://compbio.mit.edu/ad_aging_brain/
if should_process "../../data/BrainOmicsData/processed/ROSMAP/ROSMAP_processed.rds"; then
  log_message "Processing ROSMAP data..."
  run_r_script "$code_dir" "ROSMAP.R" "ROSMAP"
else
  log_message "ROSMAP data already processed!"
fi


# GSE296073 (contains GSE274829 from PMID: 40770097)
# description for GSE274829 (organoids):
# Human embryonic stem cells-induced microglia (iMG) were transplanted to 4-week-old MGE organoids.
# We conducted scRNAseq to investigate the transcriptomics of 6-week-old MGE organoids with and without iMG.
# We also used Fluorescence-activated cell sorting (FACS) to enrich GFP-labelled iMG and condcuted scRNAseq.
# title: Microglia integration into organoids recapitulates human microglial biology
# journal: Nature
# date: 2025
# paper: https://doi.org/10.1038/s41586-025-09362-8
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/40770097
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE296073
# Seurat Objects of snRNAseq data of postmortem embryonic and perinatal human sampels 
# (h_pre_peri_DY for all, INS for interneurons) as well as MGE organoids(organoid6w_DY),
# induced microglia isolated from MGE organoids (img) in manuscript 
# Zenodo: https://zenodo.org/records/15299853
# code:https://github.com/DIANKUNYU/R-script-used-for-Yu-2025
# https://github.com/codycollier/mglia-nat25
if should_process "../../data/BrainOmicsData/processed/GSE296073/GSE296073_processed.rds"; then
  log_message "Processing GSE296073 data..."
  run_r_script "$code_dir" "GSE296073.R" "GSE296073"
else
  log_message "GSE296073 data already processed!"
fi


# GSE261983 (multiome: snRNA-seq + snATAC-seq, a part of PsychENCODE project and brainSCOPE (https://brainscope.gersteinlab.org/))
# title: Single-cell genomics and regulatory networks for 388 human brains
# paper: https://doi.org/10.1126/science.adi5199
# pmid: https://pubmed.ncbi.nlm.nih.gov/38781369/
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261983
# code:
if should_process "../../data/BrainOmicsData/processed/GSE261983/GSE261983_processed.rds"; then
  log_message "Processing GSE261983 data..."
  run_r_script "$code_dir" "GSE261983.R" "GSE261983"
else
  log_message "GSE261983 data already processed!"
fi
