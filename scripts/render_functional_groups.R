# Microcosm functional-group overview.
# Usage: Rscript plot_functional_groups.R DATA_DIRECTORY OUTPUT_DIRECTORY
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(dplyr);library(jsonlite)})
args <- commandArgs(TRUE);if(length(args)!=2)stop("Supply model data and figure output directories.")
data_dir <- normalizePath(args[1]);out <- args[2];dir.create(out,recursive=TRUE,showWarnings=FALSE)
read_table <- function(n)read.csv(file.path(data_dir,n),check.names=FALSE)
raw <- read_table("model_inputs_long.csv");means <- read_table("estimated_marginal_means.csv")
tests <- read_table("factorial_anova.csv");pairs <- read_table("substrate_pairwise_Tukey.csv")
groups <- unique(raw$functional_group)
origins <- c("Full sunlight","Reduced UV","Reduced UV and blue-green")
palette <- setNames(c("#E3AC18","#77BCE8","#555555"),origins)
shapes <- setNames(c(21,24,22),origins);offsets <- setNames(c(-.20,0,.20),origins)
short_titles <- c(
 "Cellulose degradation"="Cellulose degradation",
 "Hemicellulose degradation"="Hemicellulose degradation",
 "Pectin degradation"="Pectin degradation",
 "Alpha-glucan and fructan processing"="Alpha-glucan/fructan processing",
 "Oligosaccharide and disaccharide utilization"="Oligo-/disaccharide utilization",
 "Cutin degradation"="Cutin degradation",
 "Photoreactivation and nucleotide-excision repair"="Photoreactivation and NER",
 "General DNA-damage repair and tolerance"="DNA repair and damage tolerance",
 "Enzymatic ROS defense"="Enzymatic ROS defense",
 "Carotenoid biosynthesis and modification"="Carotenoid biosynthesis",
 "Sensory photoreception"="Sensory photoreception",
 "Aromatic-compound catabolism"="Aromatic-compound catabolism")
for(v in c("inoculum","substrate")) {
 raw[[v]] <- factor(raw[[v]],levels=origins);means[[v]] <- factor(means[[v]],levels=origins)
}
set.seed(42)
raw$x <- as.integer(raw$substrate)+unname(offsets[as.character(raw$inoculum)])+runif(nrow(raw),-.045,.045)
means$x <- as.integer(means$substrate)+unname(offsets[as.character(means$inoculum)])
cells <- filter(means,kind=="inoculum-substrate cells")
stopifnot(length(groups)==12,nrow(raw)==492,nrow(cells)==108,nrow(tests)==36,nrow(pairs)==144)
fmt_p <- function(p)ifelse(p<.001,formatC(p,format="e",digits=1),sprintf("%.4f",p))
term_order <- c("substrate","inoculum","inoculum:substrate")
term_labels <- c("Substrate","Inoculum","Interaction")
make_plot <- function(g,overview=TRUE) {
 d <- filter(raw,functional_group==g);e <- filter(cells,functional_group==g)
 aa <- filter(tests,functional_group==g);aa <- aa[match(term_order,aa$term),]
 yr <- range(c(d$value,e$lower_percent,e$upper_percent));span <- diff(yr)
 aa$label <- sprintf("%s: F=%.2f; p=%s; q=%s",term_labels,aa$F,fmt_p(aa$p),fmt_p(aa$q_BH))
 aa$x <- .55;aa$y <- yr[2]+c(.44,.30,.16)*span
 p <- ggplot()+
  geom_point(data=d,aes(x=x,y=value,fill=substrate,shape=inoculum),size=1.7,stroke=.35,alpha=.33,show.legend=FALSE)+
  geom_errorbar(data=e,aes(x=x,ymin=lower_percent,ymax=upper_percent),width=.06,linewidth=.45,colour="#333333")+
  geom_point(data=e,aes(x=x,y=estimate_percent,fill=substrate,shape=inoculum),size=3.0,stroke=.7,colour="black")+
  scale_fill_manual(values=palette,breaks=origins,name="Microcosm substrate",
    guide=guide_legend(order=1,nrow=1,title.position="left",override.aes=list(shape=21,size=3.2)))+
  scale_shape_manual(values=shapes,breaks=origins,name="Inoculum origin",
    guide=guide_legend(order=2,nrow=1,title.position="left",override.aes=list(fill="white",size=3.2)))+
  scale_x_continuous(breaks=1:3,labels=c("Full\nsunlight","Reduced\nUV","Reduced UV and\nblue-green"),limits=c(.5,3.5))+
  scale_y_continuous(breaks=scales::breaks_pretty(n=4),expand=expansion(mult=c(0,0)))+
  coord_cartesian(ylim=c(yr[1]-.09*span,yr[2]+if(overview).57*span else .10*span))+
  labs(x=NULL,y="Relative abundance (%)",title=paste0(LETTERS[match(g,groups)],"  ",short_titles[g]))+
  theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),
    panel.grid.major.y=element_line(colour="#EFEFEF",linewidth=.25),
    plot.title=element_text(size=11.5,face="bold",margin=margin(b=5)),
    axis.title.y=element_text(size=10),axis.text=element_text(size=9.5,colour="#222222"),
    legend.position="bottom",legend.box="vertical",legend.title=element_text(size=11),
    legend.text=element_text(size=10.5),legend.margin=margin(2,0,2,0),
    plot.margin=margin(6,8,6,6))
 if(overview)p <- p+geom_text(data=aa,aes(x=x,y=y,label=label),inherit.aes=FALSE,
   hjust=0,vjust=.5,size=9.3/ggplot2::.pt,colour="#333333",show.legend=FALSE)
 p
}
overview <- wrap_plots(lapply(groups,make_plot),ncol=3,guides="collect")+
  plot_annotation(caption=paste("Faint points: samples; large points: estimated inoculum-substrate means with pointwise 95% CIs.",
    "q: BH adjustment across 12 groups separately for each model term.",sep="\n"),
    theme=theme(plot.caption=element_text(size=10,hjust=0),plot.margin=margin(10,8,8,8)))
overview <- overview & theme(legend.position="bottom",legend.box="vertical")

ggsave(file.path(out,"supplementary_functional_groups.pdf"),overview,width=12.5,height=16,device=cairo_pdf)
ggsave(file.path(out,"supplementary_functional_groups.png"),overview,width=12.5,height=16,dpi=240,bg="white")
