#!/bin/bash

set -e

source "functions/utils.sh"

code_dir="download"
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


# BICCN
# data: https://brainscope.gersteinlab.org/integrative_files.html

# GSE103723
# paper: https://doi.org/10.1126/sciadv.adg3754
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/37824614
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103723
# code:


# GSE104276
# paper: https://doi.org/10.1126/sciadv.adg3754
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/37824614
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104276
# code:


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


# GSE199762
# paper: https://doi.org/10.1038/s41586-023-06981-x
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/38122823
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199762
# dbGaP, accession: https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs003509.v1.p1
# code: https://github.com/massisnascimento/ECstream

# GSE204684 (SCP1859, multiome: snRNA-seq + snATAC-seq)
# data: https://singlecell.broadinstitute.org/single_cell/study/SCP1859/multi-omic-profiling-of-the-developing-human-cerebral-cortex-at-the-single-cell-level#study-download

# GSE204683 (GSE204683 (multiome): snRNA-seq + snATAC-seq (GSE204682))
# title: Multi-omic profiling of the developing human cerebral cortex at the single-cell level
# doi: https://doi.org/10.1126/sciadv.adg3754
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/37824614
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE204683
# ATAC-seq: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE204682
# CELLxGENE (RRID: SCR_021059) data (h5ad):
# https://cellxgene.cziscience.com/collections/ceb895f4-ff9f-403a-b7c3-187a9657ac2c
# code: https://doi.org/10.5281/zenodo.7703253
bash $code_dir/GSE204683.sh
bash $code_dir/GSE204682.sh


# GSE212606
# paper: https://doi.org/10.1126/sciadv.adg3754
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/37824614
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212606
# code:


# GSE217511
# paper: https://doi.org/10.1038/s41467-022-34975-2
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/36509746
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217511


# GSE67835
# paper: https://doi.org/10.1073/pnas.1507125112
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/26060301
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67835


# GSE81475
# paper: https://doi.org/10.1016/j.celrep.2016.08.038
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/27568284
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81475
# code:


# GSE97942 (contains GSE97887 + GSE97930)
# journal: Nature Biotechnology
# date: 2018
# paper: https://doi.org/10.1038/nbt.4038
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/29227469
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97942
# GSE97887 (scTHS-seq): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97887
# GSE97930 (snDrop-seq): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97930
# code:


# Li et al. 2018
# paper: https://doi.org/10.1126/science.aat7615
# pmid: https://pubmed.ncbi.nlm.nih.gov/30545854/
# data:
# code:


# Nowakowski_et_al_2017
# paper: https://doi.org/10.1126/science.aap8809
# pmid: https://pubmed.ncbi.nlm.nih.gov/29217575/
# data:
# code:


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

# AllenM1
# data: https://brain-map.org/our-research/cell-types-taxonomies/cell-types-database-rna-seq-data/human-m1-10x


# EGAD00001006049
# paper: Comprehensive cell atlas of the first-trimester developing human brain
# doi: https://doi.org/10.1126/science.adf1226
# pmid: 37824650
# data: https://ega-archive.org/datasets/EGAD00001006049
# code:


# EGAS00001006537
# paper: Single-Nuclei RNA Sequencing of 5 Regions of the Human Prenatal Brain Implicates Developing Neuron Populations in Genetic Risk for Schizophrenia
# pmid: 36150908
# data: https://ega-archive.org/studies/EGAS00001006537
# code:


# GSE178175
# paper:
# pmid:
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178175
# code:


# GSE168408
# paper:
# pmid:
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168408
# code:


# GSE144136
# paper:
# pmid:
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144136
# code:


# GSE202210
# paper:
# pmid:
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202210
# code:


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

# Ma_et_al_2022


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


# SomaMut
# papaer: https://doi.org/10.1038/s41586-025-09435-8
# code: ~
# data: https://publications.wenglab.org/SomaMut/


# PRJCA015229 (multiome: snRNA-seq + snATAC-seq, Human + Macaque)
# paper: https://doi.org/10.1016/j.xgen.2024.100703
# code: https://github.com/KIZ-SubLab/ACC-sn-Multiomes
# data: https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA015229


# ROSMAP (Religious Order Study (ROS) or the Rush Memory and Aging Project (MAP))
# paper: https://doi.org/10.1016/j.cell.2023.08.039
# code: https://github.com/mathyslab7/ROSMAP_snRNAseq_PFC/
# data: https://compbio.mit.edu/ad_aging_brain/
bash $code_dir/rosmap_processed_data.sh
bash $code_dir/rosmap_ucsc_snRNAseq.sh
bash $code_dir/rosmap_ucsc_snRNAseqsnATACseq.sh
bash $code_dir/rosmap_ucsc_snATACseq.sh
bash $code_dir/rosmap_ucsc_snATACseq_Epigenomic.sh


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
bash $code_dir/GSE296073.sh


# GSE261983 (multiome: snRNA-seq + snATAC-seq, a part of PsychENCODE project and brainSCOPE (https://brainscope.gersteinlab.org/))
# title: Single-cell genomics and regulatory networks for 388 human brains
# paper: https://doi.org/10.1126/science.adi5199
# pmid: https://pubmed.ncbi.nlm.nih.gov/38781369/
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261983
# code:
