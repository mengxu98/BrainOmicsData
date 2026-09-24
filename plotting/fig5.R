#!/usr/bin/env Rscript
# Figure 5: panel-first vector assembly from the full-atlas source tables.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
source('functions/utils.R')

genes <- c('PPP4R2','GXYLT2','KCNJ3','SMCHD1','CTSD','MRPL23',
           'FHIT','SLC10A7','KLHDC4','DACT1','DAAM1')
types <- c('Astrocytes','Endothelial cells','Mural cells','Fibroblasts',
           'Excitatory neurons','Inhibitory neurons','Microglia','Lymphocytes',
           'Neural progenitors','Differentiating oligodendrocytes',
           'Oligodendrocyte progenitor cells','Oligodendrocytes')
type_labels <- setNames(types,types)

source_dir <- Sys.getenv('BRAINOMICS_FIG5_SOURCE_DIR','figures')
out <- Sys.getenv('BRAINOMICS_FIG5_OUTPUT_DIR','figures')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
e <- fread(file.path(source_dir,'fig5_full_gene_type_source.tsv'))
if('Formal_CellType'%in%names(e) && !'CellType'%in%names(e)) setnames(e,'Formal_CellType','CellType')
s <- fread(file.path(source_dir,'fig5_full_paired_study_source.tsv'))[Paired_Donors >= 2L]
f <- fread(file.path(source_dir,'fig5_fixed_paired_donor_source.tsv'))
q <- fread(file.path(source_dir,'fig5_full_direction_statistics.tsv'))
stopifnot(nrow(e) == 132L, uniqueN(s$Dataset) == 16L,
          nrow(f) == 440L, uniqueN(f$Global_Donor_ID) == 40L,
          setequal(e$Gene,genes), setequal(s$Gene,genes),
          setequal(f$Symbol,genes), setequal(q$Gene,genes),
          setequal(f$Dataset,c('ROSMAP','SomaMut')))

e[,`:=`(Gene=factor(Gene,levels=genes),
        CellType=factor(CellType,levels=rev(types)),
        DetectionPct=100*Study_Equal_Detection)]
s[,`:=`(Gene=factor(Gene,levels=rev(genes)),
        Y=as.numeric(factor(Gene,levels=rev(genes))))]
s[,Yplot:=Y+((frank(Dataset,ties.method='first')-1)/max(1,.N-1)-0.5)*0.25,
  by=Gene]
summary_s <- s[,.(Studies=.N,Mean=mean(Mean_Difference),
                   SD=sd(Mean_Difference)),by=Gene]
summary_s[,`:=`(Y=as.numeric(Gene),
  Lower=Mean-qt(0.975,Studies-1)*SD/sqrt(Studies),
  Upper=Mean+qt(0.975,Studies-1)*SD/sqrt(Studies))]
stopifnot(nrow(q)==11L,all(q$Q_BH>=0 & q$Q_BH<=1),
          setequal(as.character(summary_s$Gene),q$Gene))
q[,`:=`(Y=match(Gene,rev(genes)),
       Q_Label=formatC(Q_BH,format='fg',digits=2,flag='#'))]

f[,`:=`(Gene=factor(Symbol,levels=rev(genes)),
        Dataset=factor(Dataset,levels=c('ROSMAP','SomaMut')))]
f[,Y:=as.numeric(Gene)]
f[,Ygroup:=Y+fifelse(Dataset=='ROSMAP',0.14,-0.14)]
set.seed(1)
f[,Yplot:=Ygroup+runif(.N,-0.04,0.04)]
summary_f <- f[,.(Mean=mean(Difference),Q1=unname(quantile(Difference,.25)),
                 Q3=unname(quantile(Difference,.75)),Donors=.N,
                 Ygroup=Ygroup[1]),by=.(Gene,Dataset)]
stopifnot(nrow(summary_f)==22L,
          all(summary_f[Dataset=='ROSMAP',Donors]==33L),
          all(summary_f[Dataset=='SomaMut',Donors]==7L))

common <- theme_minimal(base_size=9,base_family='Arial')+
  theme(panel.grid=element_blank(),axis.ticks=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(size=8.2,face='italic'),
        axis.text.x=element_text(size=7.6),
        plot.margin=margin(3,3,2,3,'mm'))

# Full cell-type names do not fit as single-line 45-degree labels in 12 columns.
# Put cell types on the vertical axis to keep every name legible without wrapping.
a <- ggplot(e,aes(x=Gene,y=CellType,size=DetectionPct,colour=Z))+
  geom_point(alpha=.95)+coord_fixed(ratio=1,clip='off')+
  scale_x_discrete(expand=expansion(add=.55))+
  scale_y_discrete(labels=type_labels,expand=expansion(add=.5))+
  scale_size_area(name='Detected (%)',max_size=4.0,limits=c(0,100),
                  breaks=c(25,50,75))+
  scale_colour_gradient2(name='Z',low='#3976A2',mid='#D9DEE0',
    high='#C75252',midpoint=0,limits=c(-3,3),oob=scales::squish,
    breaks=c(-3,0,3))+
  guides(size=guide_legend(order=1,override.aes=list(colour='#788C97')),
         colour=guide_colorbar(order=2,barheight=grid::unit(20,'mm')))+
  labs(x=NULL)+common+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5,size=7,
                                 face='italic'),
        axis.text.y=element_text(size=7,face='plain'),
        panel.border=element_rect(colour='#B9C1C4',fill=NA,linewidth=.45),
        legend.position='right',legend.title=element_text(size=7.3),
        legend.text=element_text(size=7),
        legend.key.height=grid::unit(4,'mm'),
        legend.spacing.y=grid::unit(1.2,'mm'))

b <- ggplot()+
  geom_segment(data=data.frame(Y=seq_along(genes)),
               aes(x=-3.1,xend=3.75,y=Y,yend=Y),
               colour='#E6E8E8',linewidth=.25)+
  geom_vline(xintercept=0,colour='#899498',linewidth=.4,linetype='dashed')+
  geom_vline(xintercept=4.0,colour='#D8DDDF',linewidth=.3)+
  geom_point(data=s,aes(x=Mean_Difference,y=Yplot),colour='#8AA6B2',
             size=1.15,alpha=.7)+
  geom_segment(data=summary_s,aes(x=Lower,xend=Upper,y=Y,yend=Y),
               colour='#334A55',linewidth=.75)+
  geom_point(data=summary_s,aes(x=Mean,y=Y),colour='#142B35',size=2.45)+
  geom_text(data=q,aes(x=4.4,y=Y,label=Q_Label),hjust=0,
            family='Arial',size=2.8,colour='#263238')+
  annotate('text',x=4.4,y=12.05,label='BH q',hjust=0,
           family='Arial',fontface='bold',size=2.7,colour='#263238')+
  scale_y_continuous(breaks=seq_along(genes),labels=rev(genes),
    limits=c(.45,length(genes)+1.35),expand=c(0,0))+
  scale_x_continuous(breaks=c(-2,0,2),limits=c(-3.1,7.2),expand=c(0,0))+
  labs(x=NULL)+common

dataset_colours <- c(ROSMAP='#386E9F',SomaMut='#CA7940')
c <- ggplot()+
  geom_hline(yintercept=seq_along(genes),colour='#ECEEEE',linewidth=.25)+
  geom_vline(xintercept=0,colour='#899498',linewidth=.4,linetype='dashed')+
  geom_point(data=f,aes(x=Difference,y=Yplot,colour=Dataset),
             alpha=.42,size=1.5)+
  geom_segment(data=summary_f,aes(x=Q1,xend=Q3,y=Ygroup,yend=Ygroup,
             colour=Dataset),linewidth=1.5,alpha=.9)+
  geom_point(data=summary_f,aes(x=Mean,y=Ygroup,colour=Dataset),
             shape=18,size=3)+
  scale_colour_manual(values=dataset_colours,name='S15 · PFC',
                      labels=c('ROSMAP  n=33','SomaMut  n=7'))+
  scale_y_continuous(breaks=seq_along(genes),labels=rev(genes),
    limits=c(.45,length(genes)+.55),expand=c(0,0))+
  scale_x_continuous(breaks=c(-4,-2,0,2,4),limits=c(-4.4,5.5),expand=c(0,0))+
  labs(x='OL − Microglia')+
  common+theme(legend.position='right',legend.title=element_text(size=7.5),
               legend.text=element_text(size=7.5),
               legend.key.height=grid::unit(4,'mm'))

panels <- list(a=a,b=b,c=c)
dims <- list(a=c(102,65),b=c(80,52),c=c(173,55))
for (id in names(panels)) {
  ggsave(file.path(out,paste0('fig5',id,'.pdf')),panels[[id]],
         width=dims[[id]][1],height=dims[[id]][2],units='mm',
         bg='white',device=cairo_pdf)
}
assembled <- assemble_pdf_figure(
  stem='fig5',page=c(196,145),repo=getwd(),output_dir=out,
  panels=list(
    pdf_panel('a',x=5,y=141,width=101,label='A',label_x=3,label_y=140,
              source='fig5a.pdf'),
    pdf_panel('b',x=110,y=141,width=81,label='B',label_x=108,label_y=140,
              source='fig5b.pdf'),
    pdf_panel('c',x=9,y=82,width=176,label='C',label_x=3,label_y=81,
              source='fig5c.pdf')
  )
)
export_checked <- function(command,args) {
  status <- system2(command,args)
  if (!identical(as.integer(status),0L)) stop(command,' failed with status ',status)
}
export_checked('pdftocairo',c('-svg',shQuote(assembled),
                              shQuote(file.path(out,'fig5.svg'))))
export_checked('pdftoppm',c('-f','1','-singlefile','-r','600','-png',
                           shQuote(assembled),shQuote(file.path(out,'fig5'))))
