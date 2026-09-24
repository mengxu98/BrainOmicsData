#!/usr/bin/env Rscript
# Extract the prespecified 11-gene panel from one processed source dataset.
# Run from the repository root, once per dataset with its processed object available.
suppressPackageStartupMessages({library(data.table);library(Matrix);library(SeuratObject)})
setDTthreads(2L)
source('functions/data_paths.R')
repo <- brainomics_repo_root()
source(file.path(repo,'functions/processed_object.R'))
source(file.path(repo,'functions/integration.R'))
source(file.path(repo,'functions/metadata_schema.R'))
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L)
dataset <- args[[1]]
analysis_dir <- Sys.getenv('BRAINOMICS_ANALYSIS_DIR',file.path(repo,'results/analysis_run'))
processed_dir <- Sys.getenv('BRAINOMICS_PROCESSED_DIR',brainomics_data_path('processed'))
out <- Sys.getenv('BRAINOMICS_FIG5_PANEL_DIR',file.path(repo,'results/fig5_full_panel'))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
find_unique_file <- function(root, filename) {
  candidates <- list.files(root, recursive=TRUE, full.names=TRUE)
  matches <- candidates[basename(candidates)==filename]
  if(length(matches)!=1L) stop('Expected one ',filename,' under ',root,'; found ',length(matches))
  matches[[1L]]
}
genes <- c('PPP4R2','GXYLT2','KCNJ3','SMCHD1','CTSD','MRPL23',
           'FHIT','SLC10A7','KLHDC4','DACT1','DAAM1')
metadata_file <- Sys.getenv('BRAINOMICS_METADATA_FILE',unset='')
if(!nzchar(metadata_file)) metadata_file <- find_unique_file(analysis_dir,'metadata_working.rds')
annotation_file <- Sys.getenv('BRAINOMICS_ANNOTATION_TABLE',
  unset=file.path(repo,'results','annotation','cluster_annotation.tsv'))
meta <- as.data.table(readRDS(metadata_file))
a <- fread(annotation_file)
if(!'CellType'%in%names(a) && 'Working_CellType'%in%names(a)) a[,CellType:=Working_CellType]
m <- meta[Dataset==dataset,.(Cells,Dataset,Global_Donor_ID,Canonical_Donor_ID,
                             AgeIntervalID,BrainRegion,Cluster)]
m[,CellType:=a$CellType[match(Cluster,a$Cluster)]]
stopifnot(nrow(m)>0L,!anyDuplicated(m$Cells),!anyNA(m$Canonical_Donor_ID),
          !anyNA(m$AgeIntervalID),!anyNA(m$BrainRegion),!anyNA(m$CellType))
folder <- file.path(processed_dir,dataset)
fmt <- fread(file.path(folder,'processed_object_format.tsv'))
cm <- fread(file.path(folder,fmt$Canonical_Metadata[[1L]]),
            select=c('Cells','Original_Cell_ID'))
cm[,Qualified_Cell:=paste(dataset,Original_Cell_ID,sep='::')]
stopifnot(!anyDuplicated(cm$Qualified_Cell))
obj <- load_processed_object(file.path(folder,fmt$Full_Object[[1L]]))
fm <- read_feature_metadata(file.path(folder,fmt$Feature_Metadata[[1L]]))
if (dataset == 'GSE168408') {
  # Seurat normalizes underscores in this source's feature names on import.
  normalized <- gsub('_','-',as.character(fm$Original_Feature_ID),fixed=TRUE)
  stopifnot(identical(normalized,rownames(obj[['RNA']])),!anyDuplicated(normalized))
  fm$Original_Feature_ID <- normalized
}
cw <- feature_crosswalk(fm,rownames(obj[['RNA']]),dataset)
out_parts <- list();covered <- character();k <- 0L
for (layer in Layers(obj,assay='RNA',search='^counts')) {
  ids <- Cells(obj[['RNA']],layer=layer)
  q <- cm$Qualified_Cell[match(ids,cm$Cells)]
  selected <- which(q %chin% m$Cells)
  if (!length(selected)) next
  x <- as(LayerData(obj,assay='RNA',layer=layer,cells=ids[selected]),'dgCMatrix')
  q <- cm$Qualified_Cell[match(colnames(x),cm$Cells)]
  context <- m[match(q,Cells)]
  stopifnot(!anyNA(context$Cells),identical(q,context$Cells))
  feat <- cw$Canonical_Feature_ID[match(rownames(x),cw$Original_Feature_ID)]
  stopifnot(!anyNA(feat),!anyDuplicated(feat))
  keyfields <- c('Dataset','Canonical_Donor_ID','AgeIntervalID','BrainRegion','CellType')
  groups <- unique(context[,..keyfields])
  context[,GI:=groups[context,on=keyfields,which=TRUE]]
  inc <- sparseMatrix(i=seq_len(nrow(context)),j=context$GI,x=1,
                      dims=c(nrow(context),nrow(groups)))
  groups[, `:=`(Cells=as.integer(tabulate(context$GI,nbins=nrow(groups))),
                Library_Size=as.numeric(crossprod(inc,Matrix::colSums(x))))]
  for (gene in genes) {
    idx <- match(gene,feat)
    if (is.na(idx)) {
      rows <- copy(groups)[,`:=`(Gene=gene,Measured_Cells=0L,
                                Gene_Count=0,Detected_Cells=0L)]
    } else {
      v <- as.numeric(x[idx,,drop=TRUE])
      rows <- copy(groups)[,`:=`(Gene=gene,Measured_Cells=Cells,
                                Gene_Count=as.numeric(rowsum(v,context$GI,reorder=TRUE)[as.character(seq_len(nrow(groups))),1]),
                                Detected_Cells=as.integer(rowsum(as.integer(v>0),context$GI,reorder=TRUE)[as.character(seq_len(nrow(groups))),1]))]
    }
    k <- k+1L;out_parts[[k]] <- rows
  }
  covered <- c(covered,q)
  rm(x,inc);gc(verbose=FALSE)
}
stopifnot(length(covered)==nrow(m),!anyDuplicated(covered),setequal(covered,m$Cells))
res <- rbindlist(out_parts)
agg <- res[,.(Cells=sum(Cells),Library_Size=sum(Library_Size),
             Measured_Cells=sum(Measured_Cells),Gene_Count=sum(Gene_Count),
             Detected_Cells=sum(Detected_Cells)),
           by=.(Dataset,Canonical_Donor_ID,AgeIntervalID,BrainRegion,CellType,Gene)]
agg[,`:=`(Measurement_Status=fcase(Measured_Cells==Cells,'complete',
                                   Measured_Cells==0L,'not_measured',default='partial'),
           Eligible=Cells>=20L & Library_Size>=1000)]
agg[,`:=`(CPM=fifelse(Measurement_Status=='complete' & Library_Size>0,
                      1e6*Gene_Count/Library_Size,NA_real_),
           Detection=fifelse(Measurement_Status=='complete',Detected_Cells/Cells,NA_real_))]
stopifnot(nrow(agg)==length(genes)*uniqueN(agg,by=c('Dataset','Canonical_Donor_ID','AgeIntervalID','BrainRegion','CellType')),
          all(agg$Measured_Cells<=agg$Cells),all(agg$Detected_Cells<=agg$Measured_Cells))
fwrite(agg,file.path(out,paste0('panel_',dataset,'.tsv.gz')),sep='\t')
fwrite(data.table(Dataset=dataset,Input_Cells=nrow(m),Covered_Cells=length(covered),
                  Groups=uniqueN(agg,by=c('Dataset','Canonical_Donor_ID','AgeIntervalID','BrainRegion','CellType')),
                  Complete_Gene_Groups=sum(agg$Measurement_Status=='complete'),
                  Partial_Gene_Groups=sum(agg$Measurement_Status=='partial'),
                  Missing_Gene_Groups=sum(agg$Measurement_Status=='not_measured')),
       file.path(out,paste0('source_coverage_',dataset,'.tsv')),sep='\t')
