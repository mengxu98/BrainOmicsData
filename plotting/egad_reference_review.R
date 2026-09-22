#!/usr/bin/env Rscript
# Exploratory held-out review. Does not overwrite the manuscript figures.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(ggrastr);library(ComplexHeatmap)})
source("plotting/config.R")
grDevices::pdfFonts(Arial=grDevices::pdfFonts('ArialMT')[[1]])
options(device=function(...) grDevices::cairo_pdf(file=file.path(tempdir(),'egad-layout.pdf'),family='Arial',...))
results_root<-file.path(doc,'tables/egad_mapping')
stopifnot(file.exists(file.path(results_root,'COMPLETE.json')),file.exists(file.path(results_root,'REPORT_COMPLETE.json')))
current_transfer<-TRUE;tsne<-FALSE;preview<-FALSE;reduction<-'umap';out<-'figures/egad_transfer';dir.create(out,recursive=TRUE,showWarnings=FALSE)
d<-fread(file.path(results_root,'query_predictions.tsv.gz'))
setnames(d,c('CellID','Source_CellClass','Working_CellType'),c('Cells','Source_Label','Predicted_CellType'))
stopifnot(nrow(d)==1661798L,!anyDuplicated(d$Cells),all(d$Mapping_Status=='mapped'),all(is.finite(d$UMAP_1)),all(is.finite(d$UMAP_2)))
stopifnot(all(d$Predicted_CellType %in% names(brainomics_celltype_colors)),uniqueN(d$Predicted_CellType)==12L)
# Source names are display-only harmonizations; fresh coordinates and votes unchanged.
rename<-c('Neuron'='Neurons','Neuroblast'='Neuroblasts','Neuronal IPC'='Neuronal intermediate progenitor cells','Glioblast'='Glioblasts','Oligo'='Oligodendrocyte lineage','Immune'='Immune cells','Vascular'='Vascular cells','Fibroblast'='Perivascular fibroblasts','Erythrocyte'='Erythrocytes')
d[Source_Label%in%names(rename),Source_Label:=unname(rename[Source_Label])]
# Reuse Fig. 2G's CellDimPlot rendering and final-size arrow-axis theme.
suppressPackageStartupMessages({library(SeuratObject);library(scop)})
qcolors<-brainomics_query_celltype_colors
rcolors<-c(brainomics_celltype_colors,'Not mapped (no counts)'='#999999',
 'Not mapped (no reference-feature counts)'='#999999')
counts<-Matrix::sparseMatrix(i=integer(),j=integer(),dims=c(2L,nrow(d)),
 dimnames=list(c('displayA','displayB'),d$Cells))
metadata<-data.frame(Source_Label=factor(d$Source_Label,levels=names(qcolors)),
 Predicted_CellType=factor(d$Predicted_CellType,levels=names(rcolors)),row.names=d$Cells)
metadata<-droplevels(metadata)
object<-CreateSeuratObject(counts=counts,meta.data=metadata)
embedding<-as.matrix(d[,.(UMAP_1,UMAP_2)]);rownames(embedding)<-d$Cells
object[['display']]<-CreateDimReducObject(embeddings=embedding,key='UMAP_',assay='RNA')
set.seed(2026);draw_order<-d$Cells[sample.int(nrow(d))]
plot_umap <- function(field,colors,title) {
 p<-scop::CellDimPlot(object,reduction='display',group.by=field,
  palcolor=colors,label=FALSE,seed=11,show_stat=FALSE,
  raster=TRUE,raster.dpi=c(1600,1600),pt.size=2,
  xlab=if(tsne)'t-SNE_1' else 'UMAP_1',ylab=if(tsne)'t-SNE_2' else 'UMAP_2',
  legend.title='Cell type',theme_use='theme_blank_axis',
  theme_args=list(text=element_text(family='Arial')),combine=FALSE)[[1]]
 p<-order_dim_plot_cells(p,draw_order)
 stopifnot(identical(rownames(p$data),draw_order))
 p$layers<-Filter(function(layer)!inherits(layer$geom,'GeomCustomAnn'),p$layers)
 p + scale_colour_manual(name='Cell type',values=colors,labels=identity) +
  theme_blank_axis(lab_size=5.5,axis_lwd=.6,
  xlab=if(tsne)'t-SNE_1' else 'UMAP_1',ylab=if(tsne)'t-SNE_2' else 'UMAP_2')+
  p$theme+p$coordinates+labs(title=title,subtitle=NULL)+
  guides(colour=guide_legend(ncol=1,override.aes=list(size=1.1)))+
  theme(text=element_text(family='Arial',size=6,face='plain'),
   plot.title=element_text(size=6.5,face='plain'),axis.title=element_blank(),
   legend.position='right',legend.text=element_text(size=5.5),legend.title=element_text(size=5.5),
   legend.key.height=grid::unit(35/length(colors),'mm'),legend.key.width=grid::unit(2,'mm'),
   legend.spacing.x=grid::unit(.5,'mm'),plot.margin=margin(2,2,7,5,unit='mm'))
}
ggsave(file.path(out,paste0('source_',reduction,'.pdf')),plot_umap('Source_Label',qcolors,'Original labels on reference projection'),
 device=cairo_pdf,width=90,height=56,units='mm')
ggsave(file.path(out,paste0('predicted_',reduction,'.pdf')),plot_umap('Predicted_CellType',rcolors[names(rcolors)%in%unique(d$Predicted_CellType)],'Transferred labels on reference projection'),
 device=cairo_pdf,width=90,height=56,units='mm')
if (preview) {
 preview_name<-if(current_transfer) {if(tsne)'egad_transfer_author_tsne' else 'egad_transfer_umap'} else {if(tsne)'egad_author_tsne_preview' else 'egad_umap_spca_preview'}
 assemble_pdf_figure(preview_name,c(188,62),list(
  pdf_panel('A',4,59,90,height=56,source=normalizePath(file.path(out,paste0('source_',reduction,'.pdf'))),label_x=1,label_y=61),
  pdf_panel('B',98,59,90,height=56,source=normalizePath(file.path(out,paste0('predicted_',reduction,'.pdf'))),label_x=95,label_y=61)),normalizePath('.'))
 quit(save='no')
}
c<-fread(file.path(results_root,'full_source_prediction_counts.tsv'))
setnames(c,c('Source_CellClass','Working_CellType','Row_Percent'),c('Source_Label','Predicted_CellType','Source_Fraction'))
c[Source_Label%in%names(rename),Source_Label:=unname(rename[Source_Label])];c[,Source_Fraction:=Source_Fraction/100]
mat<-matrix(0,length(qcolors),length(unique(d$Predicted_CellType)),dimnames=list(names(qcolors),names(rcolors)[names(rcolors)%in%unique(d$Predicted_CellType)]))
mat[cbind(match(c$Source_Label,rownames(mat)),match(c$Predicted_CellType,colnames(mat)))]<-c$Source_Fraction
stopifnot(all(abs(rowSums(mat)-1)<1e-8),sum(c$Cells)==1661798L)
# This is label co-occurrence, not a cosine similarity or expression correlation.
cell_mm<-3.8
ht<-Heatmap(mat,width=grid::unit(ncol(mat)*cell_mm,'mm'),height=grid::unit(nrow(mat)*cell_mm,'mm'),name='Within source (%)',col=circlize::colorRamp2(c(0,.5,1),c('white','#6BAED6','#08306B')),
 cluster_rows=FALSE,cluster_columns=FALSE,row_names_side='left',column_names_rot=55,
 row_names_gp=grid::gpar(fontsize=5.5),column_names_gp=grid::gpar(fontsize=5.5),
 row_names_max_width=grid::unit(44,'mm'),column_names_max_height=grid::unit(38,'mm'),
 row_labels=rownames(mat),
 left_annotation=rowAnnotation(simple_anno_size=grid::unit(2,'mm'),Source=rownames(mat),col=list(Source=qcolors),show_legend=FALSE,show_annotation_name=FALSE),
 top_annotation=HeatmapAnnotation(simple_anno_size=grid::unit(2,'mm'),Predicted=colnames(mat),col=list(Predicted=rcolors),show_legend=FALSE,show_annotation_name=FALSE),
 border=TRUE,rect_gp=grid::gpar(col='white',lwd=.3),
 heatmap_legend_param=list(at=c(0,.25,.5,.75,1),labels=c('0','25','50','75','100'),legend_height=grid::unit(20,'mm'),grid_width=grid::unit(2,'mm'),
   title_gp=grid::gpar(fontsize=5.5,fontface='plain'),labels_gp=grid::gpar(fontsize=5.5)),
 cell_fun=function(j,i,x,y,w,h,fill){if(mat[i,j]>=.05)grid::grid.text(sprintf('%.0f%%',100*mat[i,j]),x,y,
   gp=grid::gpar(fontsize=5,col=ifelse(mat[i,j]>.5,'white','black')))})
cairo_pdf(file.path(out,'label_agreement.pdf'),width=140/25.4,height=92/25.4,family='Arial')
draw(ht,heatmap_legend_side='right');dev.off()
a<-fread(file.path(results_root,'comparable_donor_source_agreement.tsv'));setnames(a,c('Source_Comparable','Agreement'),c('Source_Label','Concordance'))
short<-c('Neuron'='Neurons','Immune'='Immune cells','Vascular'='Vascular cells','Fibroblast'='Perivascular fibroblasts')
a[Source_Label%in%names(short),Source_Label:=unname(short[Source_Label])]
source_order<-names(qcolors)[names(qcolors)%in%a$Source_Label];a[,Source_Label:=factor(Source_Label,levels=source_order)]
summary<-a[,.(Mean=mean(Concordance)),by=Source_Label]
stopifnot(nrow(summary)==5L,sum(a$Cells)==605701L)
p<-ggplot(a,aes(Source_Label,Concordance,colour=Source_Label))+
 geom_hline(yintercept=c(0,.5,1),colour='#EEEEEE',linewidth=.25)+
 geom_point(position=position_jitter(width=.16,height=0,seed=2026),size=.7,alpha=.45)+
 geom_point(data=summary,aes(y=Mean),shape=18,size=1.4,colour='#333333')+
 scale_colour_manual(values=qcolors,guide='none')+
 scale_y_continuous(breaks=c(0,.5,1),labels=scales::label_percent(),limits=c(0,1),expand=expansion(mult=c(.02,.02)))+
 labs(x=NULL,y='Common-class concordance')+
 theme_classic(base_size=6,base_family='Arial')+
 theme(axis.line=element_line(linewidth=.25),axis.ticks=element_line(linewidth=.25),
  axis.text=element_text(size=5,colour='black'),axis.text.x=element_text(angle=65,hjust=1,vjust=1),
  plot.margin=margin(2,2,2,2,unit='mm'))
ggsave(file.path(out,'donor_concordance.pdf'),p,device=cairo_pdf,width=72,height=84,units='mm')
assemble_pdf_figure(if(current_transfer)'egad_transfer_review' else 'egad_reference_review',c(188,150),list(
 pdf_panel('A',4,147,90,height=56,source=normalizePath(file.path(out,paste0('source_',reduction,'.pdf'))),label_x=1,label_y=149),
 pdf_panel('B',98,147,90,height=56,source=normalizePath(file.path(out,paste0('predicted_',reduction,'.pdf'))),label_x=95,label_y=149),
 pdf_panel('C',4,86,110,height=84,source=normalizePath(file.path(out,'label_agreement.pdf')),label_x=1,label_y=88),
 pdf_panel('D',116,86,72,height=84,source=normalizePath(file.path(out,'donor_concordance.pdf')),label_x=113,label_y=88)),normalizePath('.'))
