#!/usr/bin/env Rscript
# Build Figure 5 source tables from all 22 extracted sources and the S15/PFC reuse analysis.
suppressPackageStartupMessages(library(data.table))
setDTthreads(2L)
source('functions/data_paths.R')
repo <- brainomics_repo_root()
input_dir <- Sys.getenv('BRAINOMICS_FIG5_PANEL_DIR',file.path(repo,'results/fig5_full_panel'))
out <- Sys.getenv('BRAINOMICS_FIG5_SOURCE_DIR',file.path(repo,'figures'))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
files <- list.files(input_dir, '^panel_.*[.]tsv[.]gz$',full.names=TRUE)
coverage_files <- list.files(input_dir, '^source_coverage_.*[.]tsv$',full.names=TRUE)
stopifnot(length(files)==22L,length(coverage_files)==22L)
a <- rbindlist(lapply(coverage_files,fread))
stopifnot(!anyDuplicated(a$Dataset),all(a$Input_Cells==a$Covered_Cells),
          sum(a$Input_Cells)==2602031L)
p <- rbindlist(lapply(files,fread))
if('Formal_CellType'%in%names(p) && !'CellType'%in%names(p)) setnames(p,'Formal_CellType','CellType')
genes <- c('PPP4R2','GXYLT2','KCNJ3','SMCHD1','CTSD','MRPL23',
           'FHIT','SLC10A7','KLHDC4','DACT1','DAAM1')
stopifnot(setequal(unique(p$Dataset),a$Dataset),setequal(unique(p$Gene),genes),
          uniqueN(p$CellType)==12L,
          !anyDuplicated(p,by=c('Dataset','Canonical_Donor_ID','AgeIntervalID',
                                'BrainRegion','CellType','Gene')),
          all(p$Measured_Cells<=p$Cells),all(p$Detected_Cells<=p$Measured_Cells))
fwrite(a,file.path(out,'full_panel_source_coverage.tsv'),sep='\t')
m <- p[,.(Cells=sum(Cells),Measured_Cells=sum(Measured_Cells),
          Groups=.N,Complete_Groups=sum(Measurement_Status=='complete'),
          Partial_Groups=sum(Measurement_Status=='partial')),
       by=.(Dataset,Gene)]
fwrite(m,file.path(out,'full_panel_measurement_by_study.tsv'),sep='\t')
d <- p[,.(Cells=sum(Cells),Library_Size=sum(Library_Size),
          Measured_Cells=sum(Measured_Cells),Gene_Count=sum(Gene_Count),
          Detected_Cells=sum(Detected_Cells),Age_Intervals=uniqueN(AgeIntervalID),
          Brain_Regions=uniqueN(BrainRegion)),
       by=.(Dataset,Canonical_Donor_ID,CellType,Gene)]
d[,`:=`(Measurement_Status=fcase(Measured_Cells==Cells,'complete',
                                 Measured_Cells==0L,'not_measured',default='partial'),
         Eligible=Cells>=20L & Library_Size>=1000 & Measured_Cells==Cells)]
d[,`:=`(CPM=fifelse(Eligible,1e6*Gene_Count/Library_Size,NA_real_),
         Detection=fifelse(Eligible,Detected_Cells/Cells,NA_real_))]
d[,Log1pCPM:=log1p(CPM)]
fwrite(d,file.path(out,'full_panel_donor_type_expression.tsv.gz'),sep='\t')
st <- d[Eligible==TRUE,.(Donors=.N,Cells=sum(Cells),
                           MeanLog1pCPM=mean(Log1pCPM),
                           MeanDetection=mean(Detection)),
        by=.(Dataset,CellType,Gene)]
fwrite(st,file.path(out,'full_panel_gene_type_by_study.tsv'),sep='\t')
eq <- st[,.(Studies=.N,Study_Equal_Log1pCPM=mean(MeanLog1pCPM),
             Study_Equal_Detection=mean(MeanDetection),
             Study_Donor_Groups=sum(Donors),Eligible_Cells=sum(Cells)),
         by=.(CellType,Gene)]
unq <- d[Eligible==TRUE,.(Unique_Donors=uniqueN(Canonical_Donor_ID)),
         by=.(CellType,Gene)]
eq <- merge(eq,unq,by=c('CellType','Gene'),all.x=TRUE)
eq[,Z:=if(.N>1L && sd(Study_Equal_Log1pCPM)>0)
            (Study_Equal_Log1pCPM-mean(Study_Equal_Log1pCPM))/sd(Study_Equal_Log1pCPM)
         else 0,by=Gene]
stopifnot(nrow(eq)==length(genes)*12L,all(is.finite(eq$Z)))
fwrite(eq,file.path(out,'full_panel_gene_type_study_equal.tsv'),sep='\t')

# Matched age/region/donor contrast across the complete resource.
base <- p[Eligible==TRUE & Measurement_Status=='complete' &
            CellType %in% c('Oligodendrocytes','Microglia')]
base[,Log1pCPM:=log1p(CPM)]
ol <- base[CellType=='Oligodendrocytes',
           .(Dataset,Canonical_Donor_ID,AgeIntervalID,BrainRegion,Gene,
             OL=Log1pCPM,OL_Cells=Cells)]
mg <- base[CellType=='Microglia',
           .(Dataset,Canonical_Donor_ID,AgeIntervalID,BrainRegion,Gene,
             Micro=Log1pCPM,Micro_Cells=Cells)]
pair <- merge(ol,mg,by=c('Dataset','Canonical_Donor_ID','AgeIntervalID','BrainRegion','Gene'))
pair[,Difference:=OL-Micro]
fwrite(pair,file.path(out,'full_panel_paired_donor_region_age.tsv.gz'),sep='\t')
donor <- pair[,.(Difference=mean(Difference),Matched_Regions=.N),
              by=.(Dataset,Canonical_Donor_ID,Gene)]
fwrite(donor,file.path(out,'full_panel_paired_donor.tsv.gz'),sep='\t')
st_pair <- donor[,.(Paired_Donors=.N,Mean_Difference=mean(Difference),
                   SD_Difference=if(.N>1L)sd(Difference) else NA_real_),
                 by=.(Dataset,Gene)]
fwrite(st_pair,file.path(out,'full_panel_paired_by_study.tsv'),sep='\t')
overall <- st_pair[Paired_Donors>=2L,
  .(Studies=.N,Paired_Study_Donor_Groups=sum(Paired_Donors),
    Study_Equal_Mean=mean(Mean_Difference),
    Study_Mean_SD=if(.N>1L)sd(Mean_Difference) else NA_real_),by=Gene]
fwrite(overall,file.path(out,'full_panel_paired_study_equal.tsv'),sep='\t')

# Two-sided exact sign tests on independent source-level mean directions.
x <- fread(file.path(out,'full_panel_paired_by_study.tsv'))[Paired_Donors>=2L]
stopifnot(setequal(x$Gene,genes),!anyDuplicated(x,by=c('Dataset','Gene')))
ans <- x[,{
  v <- Mean_Difference[Mean_Difference != 0]
  n <- length(v)
  positive <- sum(v>0)
  bt <- binom.test(positive,n,p=0.5,alternative='two.sided')
  .(Studies=.N,Nonzero_Studies=n,Positive_Studies=positive,
    Negative_Studies=sum(v<0),Paired_Study_Donor_Groups=sum(Paired_Donors),
    Study_Equal_Mean=mean(Mean_Difference),
    Study_Median=median(Mean_Difference),P_Sign=bt$p.value)
},by=Gene]
ans[,Q_BH:=p.adjust(P_Sign,method='BH')]
ans <- ans[match(genes,Gene)]
fwrite(ans,file.path(out,'full_panel_source_direction_sign_test.tsv'),sep='\t')
# Reconstruct fixed-context C from donor pseudobulk counts. This
# S15/PFC analysis is separate from the all-age/full-region extraction above.
analysis_dir <- Sys.getenv('BRAINOMICS_ANALYSIS_DIR',file.path(repo,'results/analysis_run'))
fixed_dir <- Sys.getenv('BRAINOMICS_FIG5_FIXED_DIR',unset='')
if(!nzchar(fixed_dir)) {
  downstream <- file.path(analysis_dir,'07_downstream')
  run_candidates <- c(analysis_dir,
    if(dir.exists(downstream)) list.dirs(downstream,recursive=FALSE,full.names=TRUE))
  run_candidates <- run_candidates[
    file.exists(file.path(run_candidates,'01_metadata','metadata_working.rds'))
  ]
  if(length(run_candidates)!=1L) stop('Set BRAINOMICS_FIG5_FIXED_DIR; found ',length(run_candidates),' analysis runs')
  candidates <- list.dirs(run_candidates[[1L]],recursive=FALSE,full.names=TRUE)
  candidates <- candidates[
    file.exists(file.path(candidates,'fixed_panel_donor_counts.tsv.gz')) &
    file.exists(file.path(candidates,'donor_type_eligibility.tsv'))
  ]
  if(length(candidates)!=1L) stop('Set BRAINOMICS_FIG5_FIXED_DIR; found ',length(candidates),' candidate directories')
  fixed_dir <- candidates[[1L]]
}
fixed_counts <- fread(file.path(fixed_dir,'fixed_panel_donor_counts.tsv.gz'))
fixed_eligibility <- fread(file.path(fixed_dir,'donor_type_eligibility.tsv'))
eligible_keys <- fixed_eligibility[Eligible==TRUE,
  .(Dataset,Global_Donor_ID,CellType)]
fixed_counts <- merge(fixed_counts,eligible_keys,
  by=c('Dataset','Global_Donor_ID','CellType'))[
    CellType %in% c('Oligodendrocytes','Microglia')]
fixed_counts[,Log1p_CPM:=log1p(CPM)]
ol_fixed <- fixed_counts[CellType=='Oligodendrocytes',
  .(Dataset,Global_Donor_ID,Symbol,OL_Log1p_CPM=Log1p_CPM,OL_Counts=Counts)]
micro_fixed <- fixed_counts[CellType=='Microglia',
  .(Dataset,Global_Donor_ID,Symbol,Micro_Log1p_CPM=Log1p_CPM,Micro_Counts=Counts)]
fixed <- merge(ol_fixed,micro_fixed,by=c('Dataset','Global_Donor_ID','Symbol'))
fixed[,Difference:=OL_Log1p_CPM-Micro_Log1p_CPM]
setcolorder(fixed,c('Dataset','Global_Donor_ID','Symbol','OL_Log1p_CPM',
                    'OL_Counts','Micro_Log1p_CPM','Micro_Counts','Difference'))
stopifnot(nrow(fixed)==440L,uniqueN(fixed$Global_Donor_ID)==40L,
          setequal(fixed$Symbol,genes),
          setequal(fixed$Dataset,c('ROSMAP','SomaMut')),
          fixed[,uniqueN(Global_Donor_ID),by=Dataset][Dataset=='ROSMAP',V1]==33L,
          fixed[,uniqueN(Global_Donor_ID),by=Dataset][Dataset=='SomaMut',V1]==7L)
frozen_pairs_path <- file.path(fixed_dir,'paired_ol_micro_donors.tsv')
if (file.exists(frozen_pairs_path)) {
  frozen_pairs <- fread(frozen_pairs_path)
  setorderv(fixed,c('Dataset','Global_Donor_ID','Symbol'))
  setorderv(frozen_pairs,c('Dataset','Global_Donor_ID','Symbol'))
  stopifnot(identical(as.data.frame(fixed[,.(Dataset,Global_Donor_ID,Symbol)]),
                      as.data.frame(frozen_pairs[,.(Dataset,Global_Donor_ID,Symbol)])),
            identical(fixed$OL_Counts,frozen_pairs$OL_Counts),
            identical(fixed$Micro_Counts,frozen_pairs$Micro_Counts),
            max(abs(fixed$Difference-frozen_pairs$Difference))<1e-12)
  # Preserve the source table's numeric formatting when available.
  stopifnot(file.copy(frozen_pairs_path,
                      file.path(out,'fig5_fixed_paired_donor_source.tsv'),overwrite=TRUE))
} else {
  fwrite(fixed,file.path(out,'fig5_fixed_paired_donor_source.tsv'),sep='\t')
}
# Four compact tables are the plotting inputs; the remaining tables document
# source coverage and intermediate summaries.
stopifnot(file.copy(file.path(out,'full_panel_gene_type_study_equal.tsv'),
                    file.path(out,'fig5_full_gene_type_source.tsv'),overwrite=TRUE),
          file.copy(file.path(out,'full_panel_paired_by_study.tsv'),
                    file.path(out,'fig5_full_paired_study_source.tsv'),overwrite=TRUE),
          file.copy(file.path(out,'full_panel_source_direction_sign_test.tsv'),
                    file.path(out,'fig5_full_direction_statistics.tsv'),overwrite=TRUE))
print(ans)
