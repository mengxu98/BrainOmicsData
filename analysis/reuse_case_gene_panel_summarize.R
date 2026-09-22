#!/usr/bin/env Rscript
# Summarize frozen-panel donor pseudobulk; retain inadequate coverage explicitly.
suppressPackageStartupMessages({library(data.table);library(Matrix)})
source("functions/data_paths.R")
source("functions/processed_object.R")
setDTthreads(2L)
args <- commandArgs(TRUE)
out <- if(length(args)) args[1] else "results/gene_reuse"
stopifnot(file.exists(file.path(out,"COUNTS_READY")))
unlink(file.path(out,"_SUCCESS"))
writeLines(paste0("PID=",Sys.getpid()),file.path(out,"SUMMARY_RUNNING"))
p <- fread(file.path(out,"reuse_panel_expression.tsv.gz"))
genes <- c("PPP4R2","GXYLT2","KCNJ3","SMCHD1","CTSD","MRPL23","FHIT","SLC10A7","KLHDC4","DACT1","DAAM1")
targets <- c("Oligodendrocytes","Microglia")
seed <- 20260730L; B <- 10000L
save_table <- function(x,name) fwrite(x,file.path(out,paste0(name,".tsv")),sep="\t",na="NA")
stopifnot(setequal(p$Gene,genes), !anyDuplicated(p,by=c("Scope","Dataset","Donor_ID","CellType","Gene")),
          all(p[Gene_Measured==TRUE, Gene_Detected_Cells<=Cells]),all(p[Included==TRUE,is.finite(Log1p_Gene_CPM)]))
chunks<-fread(file.path(out,"reuse_donor_pseudobulk_counts_manifest.tsv"))
count_check<-rbindlist(lapply(chunks$File,function(file) {
  chunk<-readRDS(file.path(out,file));g<-as.data.table(chunk$groups)
  stopifnot(nrow(g)==ncol(chunk$counts),nrow(chunk$cell_membership)==chunks[File==file,Cells])
  g[,Total_Count_Rebuilt:=as.numeric(Matrix::colSums(chunk$counts))]
  g[,.(Scope,Dataset,Donor_ID,CellType,Total_Count_Rebuilt)]
}))[,.(Total_Count_Rebuilt=sum(Total_Count_Rebuilt)),by=.(Scope,Dataset,Donor_ID,CellType)]
count_check<-merge(count_check,unique(p[,.(Scope,Dataset,Donor_ID,CellType,Total_Count)]),by=c("Scope","Dataset","Donor_ID","CellType"))
stopifnot(nrow(count_check)==uniqueN(p,by=c("Scope","Dataset","Donor_ID","CellType")),
  all(count_check$Total_Count_Rebuilt==count_check$Total_Count))
save_table(count_check,"reuse_full_gene_library_reconciliation")

# Reconcile every old donor row and all nine contrasts, including the historical
# rules and RNG order. The new minimum-two-pairs rule is deliberately separate.
old_dir <- brainomics_data_path("integration_25/evaluation/reuse_case_dact1_s15_pfc")
old <- fread(file.path(old_dir,"donor_celltype_pseudobulk.tsv.gz"))
current <- p[Scope=="atlas_full" & Gene=="DACT1",names(old),with=FALSE]
setorderv(old,c("Dataset","Donor_ID","CellType"));setorderv(current,c("Dataset","Donor_ID","CellType"))
stopifnot(nrow(old)==609L,nrow(current)==609L)
rowcheck <- rbindlist(lapply(names(old),function(nm) {
  if(is.numeric(old[[nm]])) {
    error <- max(abs(old[[nm]]-current[[nm]]),na.rm=TRUE)
    data.table(Field=nm,Equal=error<=1e-12,Maximum_Absolute_Difference=error)
  } else data.table(Field=nm,Equal=identical(old[[nm]],current[[nm]]),Maximum_Absolute_Difference=NA_real_)
}))
save_table(rowcheck,"legacy_DACT1_pseudobulk_reconciliation")
stopifnot(all(rowcheck$Equal))
legacy_comparators <- setdiff(sort(unique(old$CellType[old$Included])),"Oligodendrocytes")
legacy_boot <- function(paired,index) {
  set.seed(seed+index); datasets<-sort(unique(paired$Dataset)); estimates<-numeric(B)
  by_study<-lapply(datasets,function(dataset)paired[Dataset==dataset,Difference])
  for(iteration in seq_len(B)) estimates[iteration] <- mean(vapply(by_study,function(values) {
    mean(sample(values,length(values),replace=TRUE))
  },numeric(1)))
  quantile(estimates,c(.025,.975),names=FALSE)
}
legacy_pairs<-list();legacy_effects<-list();legacy_lod<-list()
for(i in seq_along(legacy_comparators)) {
  comparator<-legacy_comparators[i]
  d<-merge(old[Included==TRUE & CellType=="Oligodendrocytes",.(Dataset,Donor_ID,Target=Log1p_Gene_CPM)],
    old[Included==TRUE & CellType==comparator,.(Dataset,Donor_ID,Comparator=Log1p_Gene_CPM)],by=c("Dataset","Donor_ID"))
  d[,Difference:=Target-Comparator]
  if(nrow(d)<3L || uniqueN(d$Dataset)<2L) next
  ci<-legacy_boot(d,i);est<-mean(d[,.(Mean=mean(Difference)),by=Dataset]$Mean)
  legacy_pairs[[comparator]]<-copy(d)[,Comparator_CellType:=comparator]
  legacy_effects[[comparator]]<-data.table(Comparator_CellType=comparator,Paired_Donors=nrow(d),Datasets=uniqueN(d$Dataset),Estimate=est,CI95_Lower=ci[1],CI95_Upper=ci[2])
  legacy_lod[[comparator]]<-rbindlist(lapply(seq_len(nrow(d)),function(j) {
    remaining<-d[-j];v<-mean(remaining[,.(Mean=mean(Difference)),by=Dataset]$Mean)
    data.table(Comparator_CellType=comparator,Omitted_Dataset=d$Dataset[j],Omitted_Donor=d$Donor_ID[j],Paired_Donors=nrow(remaining),Datasets=uniqueN(remaining$Dataset),Estimate=v,Difference_From_Full=v-est)
  }))
}
legacy_new<-rbindlist(legacy_effects)
historical<-fread(file.path(old_dir,"paired_donor_effects.tsv"))
setnames(historical,"Dataset_Equal_Mean_Log1p_CPM_Difference","Estimate")
rec<-merge(legacy_new,historical,by="Comparator_CellType",suffixes=c("_Recomputed","_Historical"))
stopifnot(nrow(rec)==9L)
for(nm in c("Paired_Donors","Datasets","Estimate","CI95_Lower","CI95_Upper")) {
  rec[,(paste0(nm,"_Absolute_Difference")):=abs(get(paste0(nm,"_Recomputed"))-get(paste0(nm,"_Historical")))]
  stopifnot(all(rec[[paste0(nm,"_Absolute_Difference")]]<=1e-12))
}
save_table(rec,"legacy_DACT1_effect_reconciliation")
save_table(rbindlist(legacy_pairs),"legacy_DACT1_paired_donor_values")
lod<-rbindlist(legacy_lod)
save_table(lod,"legacy_DACT1_leave_one_donor")
save_table(lod[,.(Minimum_Estimate=min(Estimate),Maximum_Estimate=max(Estimate),
  Largest_Absolute_Influence=max(abs(Difference_From_Full)),Most_Influential_Donor=Omitted_Donor[which.max(abs(Difference_From_Full))],
  Crosses_Zero=min(Estimate)<=0 & max(Estimate)>=0),by=Comparator_CellType],"legacy_DACT1_leave_one_donor_summary")
message(Sys.time()," All 609 historical rows and nine DACT1 estimates/intervals reconciled")

# Donors, not cells, are resampled. Study membership is held fixed. sample.int
# also handles one value without R's sample(single-number) special case.
boot_means <- function(v) colMeans(matrix(v[sample.int(length(v),length(v)*B,replace=TRUE)],nrow=length(v)))
study_rows<-list();equal_rows<-list();los_rows<-list();paired_rows<-list();counter<-0L
for(scope in sort(unique(p$Scope))) for(gene in genes) for(target in targets) {
  base<-p[Scope==scope & Gene==gene]
  for(comparator in setdiff(sort(unique(base$CellType)),target)) {
    counter<-counter+1L;set.seed(seed+counter)
    d<-merge(base[Included==TRUE & CellType==target,.(Dataset,Donor_ID,Target=Log1p_Gene_CPM,Target_Cells=Cells)],
      base[Included==TRUE & CellType==comparator,.(Dataset,Donor_ID,Comparator=Log1p_Gene_CPM,Comparator_Cells=Cells)],by=c("Dataset","Donor_ID"))
    d[,`:=`(Difference=Target-Comparator,Scope=scope,Gene=gene,Target_CellType=target,Comparator_CellType=comparator)]
    paired_rows[[counter]]<-d
    ss<-list();draws<-list()
    for(study in sort(unique(base$Dataset))) {
      v<-d[Dataset==study,Difference];n<-length(v);ci<-c(NA_real_,NA_real_)
      if(n>=2L) {draws[[study]]<-boot_means(v);ci<-quantile(draws[[study]],c(.025,.975),names=FALSE)}
      ss[[study]]<-data.table(Scope=scope,Gene=gene,Target_CellType=target,Comparator_CellType=comparator,Dataset=study,
        Paired_Donors=n,Datasets=as.integer(n>0),Estimate=if(n)mean(v) else NA_real_,CI95_Lower=ci[1],CI95_Upper=ci[2],
        Target_Nonzero_Donors=sum(d[Dataset==study,Target]>0),Comparator_Nonzero_Donors=sum(d[Dataset==study,Comparator]>0),
        Both_Zero_Pairs=sum(d[Dataset==study,Target==0 & Comparator==0]),
        Expression_Information=if(!n)"no_eligible_pair" else if(all(d[Dataset==study,Target==0 & Comparator==0]))"no_detected_counts_in_either_type" else "detected_counts_present",
        Status=if(n>=2L)"eligible" else if(n==1L)"one_pair_no_interval" else "no_eligible_pair")
    }
    st<-rbindlist(ss);study_rows[[counter]]<-st;valid<-st[Status=="eligible"]
    ci<-c(NA_real_,NA_real_)
    if(nrow(valid)>=2L) ci<-quantile(Reduce(`+`,draws)/length(draws),c(.025,.975),names=FALSE)
    equal_rows[[counter]]<-data.table(Scope=scope,Gene=gene,Target_CellType=target,Comparator_CellType=comparator,
      Paired_Donors=sum(valid$Paired_Donors),Available_Paired_Donors=nrow(d),Datasets=nrow(valid),Study_Names=paste(valid$Dataset,collapse=";"),
      Estimate=if(nrow(valid))mean(valid$Estimate) else NA_real_,CI95_Lower=ci[1],CI95_Upper=ci[2],
      Target_Nonzero_Donors=sum(valid$Target_Nonzero_Donors),Comparator_Nonzero_Donors=sum(valid$Comparator_Nonzero_Donors),Both_Zero_Pairs=sum(valid$Both_Zero_Pairs),
      Studies_With_No_Detected_Counts=sum(valid$Expression_Information=="no_detected_counts_in_either_type"),
      Status=if(nrow(valid)>=2L)"eligible_multistudy" else if(nrow(valid)==1L)"single_eligible_study" else "no_eligible_study")
    if(nrow(valid)) los_rows[[counter]]<-rbindlist(lapply(valid$Dataset,function(omitted) {
      q<-valid[Dataset!=omitted]
      data.table(Scope=scope,Gene=gene,Target_CellType=target,Comparator_CellType=comparator,Omitted_Dataset=omitted,
        Paired_Donors=sum(q$Paired_Donors),Datasets=nrow(q),Estimate=if(nrow(q))mean(q$Estimate) else NA_real_,
        Difference_From_Full=if(nrow(q))mean(q$Estimate)-mean(valid$Estimate) else NA_real_,
        Status=if(nrow(q)>=2L)"multistudy_point_estimate" else if(nrow(q)==1L)"single_study_point_estimate" else "no_remaining_study")
    }))
  }
}
effects_study<-rbindlist(study_rows);effects_equal<-rbindlist(equal_rows)
save_table(effects_study,"reuse_panel_effects_by_study")
save_table(effects_equal,"reuse_panel_effects_study_equal")
save_table(rbindlist(los_rows),"reuse_panel_leave_one_study")
all_pairs<-rbindlist(paired_rows)
fwrite(all_pairs,file.path(out,"reuse_panel_paired_donor_values.tsv.gz"),sep="\t",compress="gzip")

coverage<-p[,.(All_Cells=sum(Cells),All_Donors=uniqueN(Donor_ID),Eligible_Cells=sum(Cells[Included]),Eligible_Donors=sum(Included),
  All_Gene_Counts=if(all(Gene_Measured))sum(Gene_Count) else NA_real_,Eligible_Gene_Counts=sum(Gene_Count[Included]),
  Measured_Cells=sum(Measured_Cells),Detected_Cells=sum(Gene_Detected_Cells),
  Eligible_Detected_Cells=sum(Gene_Detected_Cells[Included]),
  Mean_Donor_Log1p_CPM=if(any(Included))mean(Log1p_Gene_CPM[Included]) else NA_real_,
  Mean_Donor_Within_Type_Detection=if(any(Included))mean(Within_Type_Detection_Fraction[Included]) else NA_real_,
  Pooled_Eligible_Cell_Detection=if(any(Included))sum(Gene_Detected_Cells[Included])/sum(Cells[Included]) else NA_real_),
  by=.(Scope,Dataset,CellType,Gene)]
save_table(coverage,"reuse_panel_coverage_by_study")
save_table(p[Scope=="atlas_full",.(Cells=sum(Cells),Gene_Measured_In_All_Cells=all(Gene_Measured),
  Gene_Counts=if(all(Gene_Measured))sum(Gene_Count) else NA_real_,Detected_Cells=sum(Gene_Detected_Cells),
  Within_Cohort_Detection=if(all(Gene_Measured))sum(Gene_Detected_Cells)/sum(Cells) else NA_real_),by=.(Dataset,Gene)],"reuse_panel_gene_information_by_study")
save_table(coverage[Eligible_Donors>0,.(Datasets=.N,Donors=sum(Eligible_Donors),Cells=sum(Eligible_Cells),
  Study_Equal_Mean_Donor_Log1p_CPM=mean(Mean_Donor_Log1p_CPM),
  Study_Equal_Mean_Donor_Within_Type_Detection=mean(Mean_Donor_Within_Type_Detection)),by=.(Scope,CellType,Gene)],"reuse_panel_expression_study_equal")

# Compare label definitions in the same label-available cell pool and then
# matched donor/type groups. Native coverage remains in the table above.
selection_values<-list();selection_study<-list();selection_equal<-list()
scope_pairs<-list(source=c("source_original","atlas_source_available"),bts=c("bts_alternative","atlas_bts_available"))
for(kind in names(scope_pairs)) {
  scopes<-scope_pairs[[kind]]
  left<-p[Scope==scopes[1] & Included==TRUE,.(Dataset,Donor_ID,CellType,Gene,Alternative=Log1p_Gene_CPM,Alternative_Cells=Cells,Alternative_Detection=Within_Type_Detection_Fraction)]
  right<-p[Scope==scopes[2] & Included==TRUE,.(Dataset,Donor_ID,CellType,Gene,Atlas=Log1p_Gene_CPM,Atlas_Cells=Cells,Atlas_Detection=Within_Type_Detection_Fraction)]
  d<-merge(left,right,by=c("Dataset","Donor_ID","CellType","Gene"))
  d[,`:=`(Comparison=paste(scopes,collapse="_vs_"),Difference=Alternative-Atlas,Detection_Difference=Alternative_Detection-Atlas_Detection)]
  selection_values[[kind]]<-d
  st<-d[,{
    set.seed(seed+counter+match(.BY$Gene,genes));ci<-c(NA_real_,NA_real_)
    if(.N>=2L)ci<-quantile(boot_means(Difference),c(.025,.975),names=FALSE)
    .(Paired_Donors=.N,Estimate=mean(Difference),CI95_Lower=ci[1],CI95_Upper=ci[2],Mean_Detection_Difference=mean(Detection_Difference),
      Status=if(.N>=2L)"eligible" else "one_pair_no_interval")
  },by=.(Comparison,Dataset,CellType,Gene)]
  selection_study[[kind]]<-st
  selection_equal[[kind]]<-st[Status=="eligible",.(Paired_Donors=sum(Paired_Donors),Datasets=.N,Estimate=mean(Estimate),Mean_Detection_Difference=mean(Mean_Detection_Difference),
    Status=if(.N>=2L)"multistudy_point_estimate" else "single_study_point_estimate"),by=.(Comparison,CellType,Gene)]
}
save_table(rbindlist(selection_values),"reuse_selection_matched_donor_values")
save_table(rbindlist(selection_study),"reuse_selection_effects_by_study")
save_table(rbindlist(selection_equal),"reuse_selection_effects_study_equal")
# Compare the paired type contrast itself on shared donor sets, so a difference
# between native donor coverage is not mistaken for a label-definition effect.
selection_contrasts<-list();contrast_study<-list();contrast_equal<-list();contrast_los<-list()
for(kind in names(scope_pairs)) {
  scopes<-scope_pairs[[kind]]
  fields<-c("Dataset","Donor_ID","Gene","Target_CellType","Comparator_CellType")
  d<-merge(all_pairs[Scope==scopes[1],c(fields,"Difference"),with=FALSE],
    all_pairs[Scope==scopes[2],c(fields,"Difference"),with=FALSE],by=fields,suffixes=c("_Alternative","_Atlas"))
  d[,`:=`(Comparison=paste(scopes,collapse="_vs_"),Difference=Difference_Alternative-Difference_Atlas)]
  selection_contrasts[[kind]]<-d
  st<-d[,{
    set.seed(seed+counter+match(.BY$Gene,genes));ci<-c(NA_real_,NA_real_)
    if(.N>=2L)ci<-quantile(boot_means(Difference),c(.025,.975),names=FALSE)
    .(Paired_Donors=.N,Alternative_Estimate=mean(Difference_Alternative),Atlas_Estimate=mean(Difference_Atlas),
      Estimate=mean(Difference),CI95_Lower=ci[1],CI95_Upper=ci[2],Status=if(.N>=2L)"eligible" else "one_pair_no_interval")
  },by=.(Comparison,Dataset,Gene,Target_CellType,Comparator_CellType)]
  contrast_study[[kind]]<-st
  eligible_keys<-st[Status=="eligible",.(Comparison,Dataset,Gene,Target_CellType,Comparator_CellType)]
  eligible<-merge(d,eligible_keys,by=names(eligible_keys))
  contrast_equal[[kind]]<-eligible[,{
    set.seed(seed+counter+match(.BY$Gene,genes));ss<-split(Difference,Dataset)
    ci<-c(NA_real_,NA_real_);if(length(ss)>=2L)ci<-quantile(Reduce(`+`,lapply(ss,boot_means))/length(ss),c(.025,.975),names=FALSE)
    .(Paired_Donors=.N,Datasets=length(ss),Alternative_Estimate=mean(vapply(split(Difference_Alternative,Dataset),mean,numeric(1))),
      Atlas_Estimate=mean(vapply(split(Difference_Atlas,Dataset),mean,numeric(1))),Estimate=mean(vapply(ss,mean,numeric(1))),
      CI95_Lower=ci[1],CI95_Upper=ci[2],Status=if(length(ss)>=2L)"eligible_multistudy" else "single_eligible_study")
  },by=.(Comparison,Gene,Target_CellType,Comparator_CellType)]
  contrast_los[[kind]]<-st[Status=="eligible",{
    rbindlist(lapply(Dataset,function(omitted) {
      keep<-Dataset!=omitted
      data.table(Omitted_Dataset=omitted,Datasets=sum(keep),Paired_Donors=sum(Paired_Donors[keep]),
        Alternative_Estimate=if(any(keep))mean(Alternative_Estimate[keep]) else NA_real_,
        Atlas_Estimate=if(any(keep))mean(Atlas_Estimate[keep]) else NA_real_,
        Estimate=if(any(keep))mean(Estimate[keep]) else NA_real_)
    }))
  },by=.(Comparison,Gene,Target_CellType,Comparator_CellType)]
}
save_table(rbindlist(selection_contrasts),"reuse_selection_matched_contrast_values")
save_table(rbindlist(contrast_study),"reuse_selection_contrast_effects_by_study")
save_table(rbindlist(contrast_equal),"reuse_selection_contrast_effects_study_equal")
save_table(rbindlist(contrast_los),"reuse_selection_contrast_leave_one_study")
led<-fread(file.path(out,"reuse_cell_selection_ledger.tsv.gz"))
selection_counts<-list()
for(kind in names(scope_pairs)) for(target in targets) {
  zz<-copy(led)
  available<-if(kind=="source")zz$Source_Available else zz$BTS_Available
  type<-if(kind=="source")zz$Source_CellType else zz$BTS_CellType
  zz[,Selection_Status:=fcase(!available,"label_unavailable",Atlas_Common_Type==target & type==target,"both",
    Atlas_Common_Type==target,"atlas_only",type==target,"alternative_only",default="neither")]
  selection_counts[[paste(kind,target)]]<-zz[,.(Cells=.N),by=.(Dataset,Donor_ID,Selection_Status)][,`:=`(Comparison=kind,Target_CellType=target)]
}
save_table(rbindlist(selection_counts),"reuse_selection_cell_coverage")
save_table(data.table(Check=c("Legacy_Donor_Rows","Legacy_Contrasts","Full_Cohort_Cells","Frozen_Genes","All_Genes_All_Scopes_Complete","Minimum_Study_Pairs","Bootstrap_Seed","Bootstrap_Replicates"),
  Value=c(609,9,nrow(led),length(genes),all(p[,.N,by=.(Scope,Dataset,Donor_ID,CellType)]$N==11),2,seed,B)),"validation_checks")
writeLines(capture.output(sessionInfo()),file.path(out,"summary_sessionInfo.txt"))
writeLines(c("# Fixed S15/PFC multi-gene expression reuse analysis",
  "All 11 published Table S6 candidates are retained. DACT1 was examined previously; this is not a prospective validation panel.",
  "Cohort: current-atlas resource cells in S15 (>=60 years), prefrontal cortex, GSE202210/ROSMAP/SomaMut. Raw RNA counts; no reintegration or new cluster assignment.",
  "Pseudobulk counts are summed by study, donor and cell type. Counts RDS chunks retain ALL measured genes; sum matching groups across layers before normalization. Sample identifiers are traceability, not independent replicates.",
  "Gene CPM denominator is the sum over all measured RNA count genes. Not measured differs from a measured zero; partial measurement is excluded. Eligible groups have >=20 cells and >=1000 library counts.",
  "Within_Type_Detection_Fraction is the fraction of cells in the donor/type group with a positive raw count for that gene. It is not the fraction of donors with nonzero pseudobulk, nor Table S6's fraction of all brain cells.",
  "Effects are paired donor differences in log1p(CPM), target minus comparator. A study requires >=2 pairs; cross-study estimates weight eligible studies equally. All insufficient-coverage contrasts are retained with status; do not replace missing effects by zero.",
  "Pointwise 95% intervals use 10000 donor resamples within each study, conditional on observed studies. No new P values, multiplicity-adjusted inference, or gene-wise significance claims. Leave-one-study reports point estimates only.",
  "Expression-information fields separately flag contrasts with no detected counts in either type. A zero effect and zero-width interval for all-zero pairs is mathematically degenerate and is not evidence of biological equivalence. Sparse measured genes remain in the frozen panel; see reuse_panel_gene_information_by_study.tsv.",
  "Main panel: atlas_full; all eleven OL-minus-Microglia contrasts are the fixed primary visualization. Both targets versus all other types remain in source tables regardless of direction or coverage.",
  "source_original versus atlas_source_available uses original ROSMAP/SomaMut labels in the same available pool; broad atlas inhibitory counts are merged BEFORE CPM. GSE202210 BTSatlas is a separately named alternative-label sensitivity, not original source labels or a new independent reference.",
  "Selection tables retain native coverage, both/atlas_only/alternative_only/label_unavailable/neither cells, and expression in matched eligible donor/type groups. Expression differences describe label-selection dependence, not annotation accuracy.",
  "Historical DACT1 outputs are untouched. Legacy tables reproduce all 609 rows and all nine old effect estimates/intervals exactly (tolerance 1e-12), including the old study-inclusion rule. Their estimates need not equal the new >=2-pairs-per-study results.",
  "The case provides expression context, not GWAS replication, AD case-control association, mediation, causality, general age preservation, or validation of all 75 clusters.",
  "Input/provenance/contract/coverage/session tables record source identities and actual denominators. Table S6 does not assign gene-specific biomarker or genetic-evidence categories; missing provenance is not inferred from DACT1."
),file.path(out,"README.md"))
code_files<-c("analysis/reuse_case_gene_panel.R","analysis/reuse_case_gene_panel_summarize.R")
files<-list.files(out,full.names=TRUE);files<-files[!basename(files)%chin%c("manifest.tsv","_SUCCESS","SUMMARY_RUNNING","_INCOMPLETE","summary.log")]
writeLines(paste0("Completed=",Sys.time()),file.path(out,"_SUCCESS"))
unlink(file.path(out,c("SUMMARY_RUNNING","_INCOMPLETE")))
message(Sys.time()," Summary complete: ",out)
