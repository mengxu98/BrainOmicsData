source("plotting/config.R")
fig2_theme <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#2B2B2B"),
      axis.ticks = element_line(linewidth = 0.3, colour = "#2B2B2B"),
      axis.text = element_text(size = 7, colour = "#2B2B2B"),
      strip.text = element_text(size = 7),
      panel.grid.major.y = element_line(linewidth = 0.18, colour = "#E4E7EB"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = base_size + 0.8, face = "bold"),
      plot.subtitle = element_text(size = base_size + 0.3, colour = "#4A4A4A")
    )
}
fig2_row_theme <- theme(
  text = element_text(face = "plain"),
  axis.title = element_text(size = 6, face = "plain"),
  axis.text = element_text(size = 5.5, face = "plain"),
  strip.text = element_text(size = 5.5, face = "plain"),
  plot.title = element_text(size = 6.5, face = "plain")
)

# Same former Figure 2 aesthetics. No inherited estimates, significance stars or resampling intervals.
full_lisi_complete <- file.exists(file.path(doc,"tables/full_lisi/COMPLETE.json"))
if(full_lisi_complete) {
 all_lisi <- fread(file.path(doc,"tables/full_lisi_dataset_summary.tsv"))
 x <- all_lisi[Space=="latent50" & Label_Scheme=="Source_Full"]
 setnames(x,"Mean","cLISI")
} else x <- fread(file.path(doc,"tables/full_source_clisi_all_method_dataset_summary.tsv"))
x[,Method:=factor(Method,levels=method_levels)]
stopifnot(nrow(x)==22L*4L,!anyDuplicated(x[,.(Dataset,Method)]))
p_source_lisi<-ggplot(x,aes(as.integer(Method),cLISI,group=Method))+geom_boxplot(aes(fill=Method),width=.58,outlier.shape=NA,alpha=.25,linewidth=.35)+geom_point(aes(colour=Method),size=.65,alpha=.6,position=position_jitter(width=.12,height=0,seed=2026))+scale_x_continuous(breaks=1:4,labels=method_levels)+scale_colour_manual(values=method_colors,guide="none")+scale_fill_manual(values=method_colors,guide="none")+labs(x=NULL,y="Source-label cLISI",title=if(full_lisi_complete)"Source cLISI (50D)" else "Source cLISI (UMAP2)")+fig2_theme(7)+fig2_row_theme+theme(axis.text.x=element_text(angle=30,hjust=1))
p_effect<-ggplot()+annotate("text",x=0,y=0,label="Full-cohort iLISI\nPending",size=2.2,family="Arial",colour="#64748B")+xlim(-1,1)+ylim(-1,1)+theme_void()+labs(title="Dataset mixing")+theme(plot.title=element_text(family="Arial",size=6.5))
if(full_lisi_complete) {
 mixing <- all_lisi[Space=="latent50" & Label_Scheme=="Dataset"]
 stopifnot(nrow(mixing)==22L*4L,!anyDuplicated(mixing[,.(Dataset,Method)]),all(is.finite(mixing$Mean)))
 mixing[,Method:=factor(Method,levels=method_levels)]
 p_effect <- ggplot(mixing,aes(Method,Mean,group=Method)) +
  geom_boxplot(aes(fill=Method),width=.58,outlier.shape=NA,alpha=.25,linewidth=.35) +
  geom_point(aes(colour=Method),size=.65,alpha=.6,position=position_jitter(width=.12,height=0,seed=2026)) +
  scale_colour_manual(values=method_colors,guide="none") + scale_fill_manual(values=method_colors,guide="none") +
  labs(x=NULL,y="Dataset iLISI",title="Dataset mixing (50D)") + fig2_theme(7)+fig2_row_theme +
  theme(axis.text.x=element_text(angle=30,hjust=1))
}
v<-fread(file.path(doc,"tables/age_signal/donor_structure_by_dataset.tsv"));v[,Method:=factor(Method,levels=method_levels)];ds<-v[,.(Estimate=mean(Spearman_Donor_Distance_to_Raw)),by=Method]
p_donor<-ggplot(v,aes(Method,Spearman_Donor_Distance_to_Raw,colour=Method))+geom_point(size=.6,alpha=.4,position=position_jitter(width=.12,height=0,seed=2026))+geom_point(data=ds,aes(y=Estimate),size=1.8)+scale_colour_manual(values=method_colors,guide="none")+labs(x=NULL,y="Distance correlation with Raw",title="Donor-structure preservation")+fig2_theme(7)+fig2_row_theme+theme(axis.text.x=element_text(angle=30,hjust=1))
age<-fread(file.path(doc,"tables/age_signal/age_metrics_summary.tsv"))[Distance=="euclidean"];age[,Method:=factor(Method,levels=method_levels)]
folds<-fread(file.path(doc,"tables/age_signal/age_metrics_by_dataset.tsv"))[Distance=="euclidean"];folds[,Method:=factor(Method,levels=method_levels)]
p_age<-ggplot(age,aes(Method,Dataset_Equal_MAE,colour=Method))+geom_hline(yintercept=unique(age$Dataset_Equal_Baseline_MAE),linewidth=.3,linetype=2,colour="#777777")+geom_point(data=folds,aes(y=MAE),size=.7,alpha=.5)+geom_point(size=2)+scale_colour_manual(values=method_colors,guide="none")+labs(x=NULL,y="Age prediction MAE (years)",title="Age prediction")+fig2_theme(7)+fig2_row_theme+theme(axis.text.x=element_text(angle=30,hjust=1))
p_validation_row<-patchwork::wrap_plots(p_effect,p_source_lisi,p_donor,p_age,nrow=1,widths=c(1,1,2,1))+patchwork::plot_annotation(tag_levels=list(c("B","C","D","E")))
p_validation_row<-p_validation_row & theme(plot.margin=margin(2,2,2,2,unit="mm"),plot.tag=element_text(family="Arial",size=9,face="plain"),plot.tag.position=c(0,1))
ggsave("figures/fig2_validation_row.pdf",p_validation_row,width=168,height=48,units="mm",device=cairo_pdf,family="Arial",bg="white")
