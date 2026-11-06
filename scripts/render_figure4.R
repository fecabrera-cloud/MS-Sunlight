suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(jsonlite)})
script_file<-normalizePath(sub("^--file=","",commandArgs(trailingOnly=FALSE)[grep("^--file=",commandArgs(trailingOnly=FALSE))]))
base<-dirname(dirname(script_file))
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=2)stop("Usage: Rscript plot_figure4.R DATA_DIRECTORY OUTPUT_DIRECTORY")
data_dir<-normalizePath(args[1],mustWork=TRUE);out<-args[2]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
config<-read_json(file.path(base,"docs","figure4_configuration.json"))
origins<-c("Full sunlight","Reduced UV","Reduced UV and blue-green")
substrate_colours<-unlist(config$substrates);inoculum_shapes<-unlist(config$inoculum_shapes)
components<-c("Inoculum","Substrate","Interaction","Residual")
component_colours<-unlist(config$variance_partition$component_colours)
make_inset <- function(partition) {
  partition$xmax <- cumsum(partition$proportion)
  partition$xmin <- c(0,head(partition$xmax,-1))
  partition$x <- (partition$xmin+partition$xmax)/2
  partition$label <- sprintf("%.1f",100*partition$proportion)
  # Two legend rows: inoculum / substrate, interaction / residual.
  legend <- data.frame(component=components,x=c(.045,.55,.045,.55),y=c(.30,.30,.10,.10))
  g <- ggplot()+
    annotate("rect",xmin=.002,xmax=.998,ymin=.002,ymax=.998,fill="white",colour="#BDBDBD",linewidth=.3)+
    annotate("text",x=.025,y=.91,label="Variance explained (%)",hjust=0,size=3.2,fontface="bold")+
    geom_rect(data=partition,aes(xmin=.025+.95*xmin,xmax=.025+.95*xmax,ymin=.51,ymax=.73,fill=component),
      colour="white",linewidth=.35)+
    geom_text(data=partition[partition$component!="Residual",],aes(x=.025+.95*x,y=.62,label=label),
      colour="white",size=2.95)+
    geom_text(data=partition[partition$component=="Residual",],aes(x=.025+.95*x,y=.62,label=label),
      colour="#222222",size=2.95)+
    geom_rect(data=legend,aes(xmin=x,xmax=x+.055,ymin=y-.045,ymax=y+.045,fill=component),colour=NA)+
    geom_text(data=legend,aes(x=x+.077,y=y,label=component),hjust=0,size=2.95)+
    scale_fill_manual(values=component_colours,guide="none")+
    coord_cartesian(xlim=c(0,1),ylim=c(0,1),expand=FALSE,clip="off")+theme_void()+
    theme(plot.margin=margin(0,0,0,0),plot.background=element_rect(fill="white",colour=NA))
  g
}
scores<-read.csv(file.path(data_dir,"plotted_ordination_scores.csv"))
ellipses<-read.csv(file.path(data_dir,"pooled_substrate_ellipses.csv"))
partitions<-read.csv(file.path(data_dir,"sequential_variance_partition.csv"))
statistics<-read.csv(file.path(data_dir,"printed_statistics.csv"))
panels<-list();footer_rows<-list()
summaries<-read.csv(file.path(data_dir,"panel_summary_labels.csv"))
for(marker in c("ITS","16S")) {
  xy<-scores[scores$marker==marker,];e<-ellipses[ellipses$marker==marker,]
  p<-partitions[partitions$marker==marker,];a<-statistics[statistics$marker==marker,]
  xy$substrate<-factor(xy$substrate,levels=origins);xy$inoculum<-factor(xy$inoculum,levels=origins)
  stress<-summaries$NMDS_stress[summaries$marker==marker]
  stopifnot(nrow(xy)==if(marker=="ITS")40 else 36,nrow(a)==3)
  xr<-range(c(xy$Axis1,e$Axis1));yr<-range(c(xy$Axis2,e$Axis2))
  xl<-xr+c(-1,1)*.075*diff(xr);yl<-c(yr[1]-.075*diff(yr),yr[2]+.48*diff(yr))
  footer<-sprintf("N = %d; NMDS stress = %.3f",nrow(xy),stress)
  footer_rows[[marker]]<-data.frame(marker,N=nrow(xy),NMDS_stress=stress,label=footer,
    x=xl[2]-.025*diff(xl),y=yl[1]+.025*diff(yl),text_size_mm=2.75)
  g<-ggplot(xy,aes(Axis1,Axis2,fill=substrate,shape=inoculum))+
    geom_polygon(data=e,aes(Axis1,Axis2,group=substrate,fill=substrate),inherit.aes=FALSE,
      colour=NA,alpha=.10,show.legend=FALSE)+
    geom_path(data=e,aes(Axis1,Axis2,group=substrate,colour=substrate),inherit.aes=FALSE,
      linewidth=.65,linetype="dotted",show.legend=FALSE)+
    geom_point(size=3.1,stroke=.55,colour="black")+
    annotate("text",x=rep(xl[1]+.025*diff(xl),3),y=yl[1]+c(.945,.875,.805)*diff(yl),
      label=a$label,hjust=0,vjust=.5,size=3.15)+
    annotate("text",x=xl[2]-.025*diff(xl),y=yl[1]+.025*diff(yl),label=footer,
      hjust=1,vjust=0,size=2.75)+
    scale_fill_manual(values=substrate_colours,breaks=origins,drop=FALSE,name="Microcosm substrate",
      guide=guide_legend(order=1,override.aes=list(shape=21,size=3.1,colour="black")))+
    scale_shape_manual(values=inoculum_shapes,breaks=origins,drop=FALSE,name="Inoculum origin",
      guide=guide_legend(order=2,override.aes=list(fill="#EEEEEE",colour="black",size=3.1)))+
    scale_colour_manual(values=substrate_colours,guide="none")+
    coord_cartesian(xlim=xl,ylim=yl,expand=FALSE)+theme_bw(base_size=13)+
    theme(panel.grid=element_blank(),plot.title=element_text(size=14,face="bold"),
      plot.margin=margin(7,9,5,7),legend.position="bottom",legend.title=element_text(size=11),
      legend.text=element_text(size=10.5),legend.box="vertical",legend.spacing.y=grid::unit(.05,"cm"),
      legend.margin=margin(0,0,0,0))+
    labs(x="NMDS1",y="NMDS2",title=if(marker=="ITS")"A  Fungi" else "B  Bacteria")
  panels[[marker]]<-g+inset_element(make_inset(p),left=.43,bottom=.76,right=.985,top=.985,
    align_to="panel",on_top=TRUE,clip=TRUE)
}
figure<-wrap_plots(panels,nrow=1,guides="collect")+
  plot_annotation(theme=theme(plot.margin=margin(7,9,6,7)))
figure<-figure & theme(legend.position="bottom",legend.box="vertical")
ggsave(file.path(out,"Figure4.pdf"),figure,width=11.8,height=6.6,device=cairo_pdf)
ggsave(file.path(out,"Figure4.png"),figure,width=11.8,height=6.6,dpi=240)
writeLines(capture.output(sessionInfo()),file.path(out,"plot_R_session_info.txt"))
cat("Rendered Figure4.pdf and Figure4.png from final model results.\n")
