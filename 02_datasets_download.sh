#!/bin/bash
# Download one original source per script in download/ and keep the raw
# archives that the preprocessing stage reads.


set -e

source "functions/utils.sh"
BRAINOMICS_STAGE=02_source_download
source "functions/pipeline_lib.sh"

code_dir="download"
overwrite="${1:-F}"
check_command "$BRAINOMICS_RSCRIPT"
check_command "$BRAINOMICS_PYTHON"

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
# One download script per released source sits in download/ with concrete source
# URLs. Entries below without a download command record sources that are not part
# of the released cohort; citation and access details are in
# data/source_access_summary.tsv.

# BICCN
# data: https://brainscope.gersteinlab.org/integrative_files.html

# GSE103723
# paper: https://doi.org/10.1126/sciadv.adg3754
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/37824614
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103723
# code:
# GSE104276
# paper: https://doi.org/10.1038/nature25980
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/29539641
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104276
# code:
bash $code_dir/GSE104276.sh


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
# paper: https://doi.org/10.1016/j.neuron.2021.10.036
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/34798047
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186538
# code:
bash $code_dir/GSE186538.sh


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
# GSE212606
# paper: https://doi.org/10.1126/sciadv.adg3754
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/37824614
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212606
# code:
bash $code_dir/GSE212606.sh


# GSE217511
# paper: https://doi.org/10.1038/s41467-022-34975-2
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/36509746
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217511
bash $code_dir/GSE217511.sh


# GSE67835
# paper: https://doi.org/10.1073/pnas.1507125112
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/26060301
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67835
bash $code_dir/GSE67835.sh


# GSE81475
# paper: https://doi.org/10.1016/j.celrep.2016.08.038
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/27568284
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81475
# code:
bash $code_dir/GSE81475.sh


# GSE97942 (contains GSE97887 + GSE97930)
# journal: Nature Biotechnology
# date: 2018
# paper: https://doi.org/10.1038/nbt.4038
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/29227469
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97942
# GSE97887 (scTHS-seq): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97887
# GSE97930 (snDrop-seq): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97930
# code:
bash $code_dir/GSE97942.sh


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


# AllenM1
# data: https://brain-map.org/our-research/cell-types-taxonomies/cell-types-database-rna-seq-data/human-m1-10x
bash $code_dir/AllenM1.sh


# EGAS00001006537
# paper: Single-Nuclei RNA Sequencing of 5 Regions of the Human Prenatal Brain Implicates Developing Neuron Populations in Genetic Risk for Schizophrenia
# pmid: 36150908
# data: https://ega-archive.org/studies/EGAS00001006537
# code:
bash $code_dir/EGAS00001006537.sh


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
bash $code_dir/GSE168408.sh


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
bash $code_dir/GSE207334.sh

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
bash $code_dir/HYPOMAP.sh


# SomaMut
# paper: https://doi.org/10.1038/s41586-025-09435-8
# code: ~
# data: https://publications.wenglab.org/SomaMut/
bash $code_dir/SomaMut.sh


# PRJCA015229 (multiome: snRNA-seq + snATAC-seq, Human + Macaque)
# paper: https://doi.org/10.1016/j.xgen.2024.100703
# code: https://github.com/KIZ-SubLab/ACC-sn-Multiomes
# data: https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA015229
bash $code_dir/PRJCA015229.sh


# ROSMAP (Religious Order Study (ROS) or the Rush Memory and Aging Project (MAP))
# paper: https://doi.org/10.1016/j.cell.2023.08.039
# code: https://github.com/mathyslab7/ROSMAP_snRNAseq_PFC/
# data: https://compbio.mit.edu/ad_aging_brain/
bash $code_dir/rosmap_ucsc_snRNAseq.sh
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
# induced microglia isolated from MGE organoids (img)
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
# Catching et al. 2026
# title: Single-nucleus multiome analysis in the human prefrontal cortex identifies gene expression and cis-regulatory elements associated with aging
# modality: 10x Multiome (snRNA-seq + snATAC-seq)
# species: Homo sapiens
# paper: https://doi.org/10.1016/j.celrep.2026.117110
# pmid: https://pubmed.ncbi.nlm.nih.gov/41832957/
# data (latest processed release): https://zenodo.org/records/20834804
# data (paper accession): https://zenodo.org/records/18394349
# raw data (controlled): dbGaP phs004202.v1.p1 and NDA collection 3151
# code: https://github.com/NIH-CARD/scMAVERICS
# download: all public processed RNA, ATAC, and supplementary metadata files

# Wang et al. 2025
# title: Molecular and cellular dynamics of the developing human neocortex
# modality: 10x Multiome (snRNA-seq + snATAC-seq); MERFISH is also reported
# species: Homo sapiens
# paper: https://doi.org/10.1038/s41586-024-08351-7
# pmid: https://pubmed.ncbi.nlm.nih.gov/39779846/
# data: https://cellxgene.cziscience.com/collections/ad2149fc-19c5-41de-8cfe-44710fbada73
# raw data: https://assets.nemoarchive.org/dat-oiif74w
# processed data: https://doi.org/10.5061/dryad.2280gb612
# code: https://github.com/complexdisease/Human_Cortex_Dev_Multiome
# download: complete public CELLxGENE RNA h5ad plus ATAC fragments and index
# integration: RNA only; exclude GW27-2-7-18, NIH-5900, NIH-5554, and NIH-4341
# because they overlap the Velmeshev reference dataset.
bash $code_dir/Wang_2025.sh


# GSE294786 / Klavert et al. 2026
# title: Charting the human-specific properties of gene expression networks in the infant prefrontal cortex
# modality: snRNA-seq
# species: Homo sapiens, Pan troglodytes, and Macaca mulatta
# paper: https://doi.org/10.1126/sciadv.aea3316
# pmid: https://pubmed.ncbi.nlm.nih.gov/42234754/
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294786
# code:
# download: GEO provides one combined cross-species counts file and one combined
# metadata file, so both complete files must be downloaded; no separate
# non-human-only asset is downloaded.
# integration: human RNA only; retain donor 1385 in the reference and exclude
# the matching HBCC-1385 donor from the recorded but unused Catching dataset.

bash $code_dir/GSE294786.sh


# Velmeshev et al. 2023
# title: Single-cell analysis of prenatal and postnatal human cortical development
# modality: snRNA-seq; the public integrated CELLxGENE object also contains
# 10x Multiome-derived RNA from other source studies
# species: Homo sapiens
# paper: https://doi.org/10.1126/science.adf0834
# pmid: https://pubmed.ncbi.nlm.nih.gov/37824647/
# data: https://cellxgene.cziscience.com/collections/bacccb91-066d-4453-b70e-59de0b4598cd
# raw data: https://assets.nemoarchive.org/dat-3ah9h9x
# code: https://github.com/velmeshevlab/dev_hum_cortex
# download: complete public CELLxGENE object before source or donor filtering
# integration: use Dataset == "Velmeshev" in the reference and exclude donor
# 5936 because it is already present in GSE204683. Imported source datasets in
# the complete CELLxGENE object are retained in the full processed object but
# excluded from the analysis object as duplicates of existing inputs.
bash $code_dir/Velmeshev_2023.sh


# Clarence et al. 2025
# title: Multiomic single-cell profiling identifies critical regulators of postnatal brain
# modality: 10x Multiome (snRNA-seq + snATAC-seq)
# species: Homo sapiens
# paper: https://doi.org/10.1038/s41588-025-02083-8
# pmid: https://pubmed.ncbi.nlm.nih.gov/39962241/
# data (public RNA view): https://cellxgene.cziscience.com/collections/f406a653-c079-4bf9-aab6-85846c27571d
# raw multiome data (controlled): https://nda.nih.gov/edit_collection.html?id=5371
# code: https://github.com/DiseaseNeuroGenomics/snMultiome
# note: eight of ten donors overlap GSE204683; the two non-overlapping donors
# are adults. The public source does not provide a complete open RNA+ATAC
# processed pair, so no downloader is enabled.


# Palmer et al. 2026
# title: Single-cell multiomic human brain atlas reveals regulatory drivers of cortical regionality
# modality: SNARE-seq2 (paired snRNA-seq + snATAC-seq) and DART-FISH
# species: Homo sapiens
# paper: https://doi.org/10.1038/s41467-026-69368-2
# pmid: https://pubmed.ncbi.nlm.nih.gov/41723114/
# data (SNARE-seq2): https://data.nemoarchive.org/biccn/lab/zhang_kun/multimodal/sncell/
# data (DART-FISH): https://doi.org/10.35077/g.1179
# code: https://github.com/ypauling/human_brain_atlas_cortex_regionality
# note: only H19.30.004 and UW7118 are net-new donors after comparison with
# AllenM1 and HYPOMAP; UW7118 metadata remain incomplete. No downloader is
# enabled until the exact NeMO file manifest is fixed.


# GSE264624 / Thompson et al. 2025
# title: An integrated single-nucleus and spatial transcriptomics atlas reveals the molecular landscape of the human hippocampus
# modality: snRNA-seq, Visium spatial transcriptomics, and Visium-SPG
# species: Homo sapiens
# paper: https://doi.org/10.1038/s41593-025-02022-0
# pmid: https://pubmed.ncbi.nlm.nih.gov/40739059/
# data (snRNA-seq): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE264624
# data (spatial): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE230782
# code: https://github.com/Erik-D-Nelson/ARG_HPC_snRNAseq
# note: ten adult neurotypical donors; this does not address the paediatric
# coverage gap and is retained as a future regional expansion.


# GSE255968 / Kim et al. 2025
# title: An expanded subventricular zone supports postnatal cortical interneuron migration in gyrencephalic brains
# modality: snRNA-seq
# species: Homo sapiens and Sus scrofa; GSE255968 contains the human data
# paper: https://doi.org/10.1038/s41593-025-01987-2
# pmid: https://pubmed.ncbi.nlm.nih.gov/40659844/
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255968
# code:
# note: three prenatal human samples are public, but complete clinical-to-donor
# mapping is available for only two cases. Pig data are recorded but are not
# downloaded.


# Ament et al. 2023
# title: A single-cell genomic atlas for maturation of the human cerebellum during early childhood
# modality: snRNA-seq
# species: Homo sapiens
# paper: https://doi.org/10.1126/scitranslmed.ade1283
# pmid: https://pubmed.ncbi.nlm.nih.gov/37824600/
# data: https://assets.nemoarchive.org/collection/nemo%3Adat-wtqs68o
# code:
# note: the non-inflamed cohort includes substantial clinical conditions and
# the public sample table lacks stable donor identifiers; do not use in the
# primary reference.


# GSE280569 / Steyn et al. 2024
# title: A temporal cortex cell atlas highlights gene expression dynamics during human brain maturation
# modality: snRNA-seq and spatial transcriptomics
# species: Homo sapiens
# paper: https://doi.org/10.1038/s41588-024-01990-6
# pmid: https://pubmed.ncbi.nlm.nih.gov/39567748/
# data (snRNA-seq): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280569
# data (spatial): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280570
# code:
# note: all new donors are surgical or pathological cases, and the integrated
# object also contains all four Thrupp donors; do not use in the primary
# reference or count those donors twice.


# GSE153807 / Thrupp et al. 2020
# title: Single-Nucleus RNA-Seq Is Not Suitable for Detection of Microglial Activation Genes in Humans
# modality: matched scRNA-seq and snRNA-seq
# species: Homo sapiens
# paper: https://doi.org/10.1016/j.celrep.2020.108189
# pmid: https://pubmed.ncbi.nlm.nih.gov/32997994/
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153807
# code:
# note: four neurosurgical epilepsy donors; the same donors are already present
# in the Steyn integrated object. Do not use in the primary reference.


# GSE301953 / Lange et al. 2026
# title: A single-nucleus transcriptomic atlas of human basal ganglia during development forwarding diagnosis and therapy of pediatric movement disorders
# modality: snRNA-seq
# species: Homo sapiens
# paper (preprint): https://doi.org/10.64898/2026.06.04.26354648
# pmid: not available
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301953
# data (CELLxGENE): https://cellxgene.cziscience.com/collections/3332ad3e-8599-45af-bef8-17f16cf24245
# code:
# note: preprint cohort with several clinical conditions; do not use in the
# primary reference.


# Luquez et al. 2026
# title: Cell-type signatures of Alzheimer's disease shared across population groups
# modality: 10x Multiome (snRNA-seq + snATAC-seq)
# species: Homo sapiens
# paper: https://doi.org/10.1038/s41586-026-10793-0
# pmid: https://pubmed.ncbi.nlm.nih.gov/42457956/
# data: https://www.synapse.org/Synapse:syn53649093
# code:
# note: controlled access and Alzheimer disease-focused; record only and do not
# download.


# HRA006553 / Zhang et al. 2025
# title: Single-cell spatiotemporal transcriptomic and chromatin accessibility profiling in developing postnatal human and macaque prefrontal cortex
# modality: snRNA-seq, snATAC-seq, and spatial transcriptomics
# species: Homo sapiens and Macaca fascicularis
# paper: https://doi.org/10.1038/s41593-025-02150-7
# pmid: https://pubmed.ncbi.nlm.nih.gov/41381947/
# human data (controlled): https://ngdc.cncb.ac.cn/gsa-human/browse/HRA006553
# macaque data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305245
# code: https://github.com/Wu-lab-code-repository/PostnatalPFC
# note: HRA006553 has no open processed human matrix. GSE305245 contains only
# Macaca fascicularis RNA, ATAC, and spatial data, so neither source is
# downloaded for this human reference.
