suppressPackageStartupMessages({library(ggplot2);library(ggrepel);library(patchwork);library(jsonlite)})
script<-normalizePath(sub("^--file=","",commandArgs(FALSE)[grep("^--file=",commandArgs(FALSE))]))
base<-dirname(dirname(script));args<-commandArgs(TRUE);stopifnot(length(args)==1);out<-normalizePath(args[1])
data_out<-out
source(file.path(base,"scripts","functional_ordination_style.R"))
xy<-read.csv(file.path(out,"ordination_scores.csv"));loadings<-read.csv(file.path(out,"pcoa_loadings.csv"))
variance<-read.csv(file.path(out,"ordination_variance.csv"))
fig<-make_microcosm_pcoa(xy,loadings,variance)
legend_labels<-c("Full sunlight","Reduced UV","Reduced UV and blue-green")
for(aesthetic in c("fill","shape")) {
  legend_scale<-fig$scales$get_scales(aesthetic);legend_scale$labels<-legend_labels
}
original<-ggplot_build(fig)
stats<-read.csv(file.path(out,"sequential_permanova.csv"));config<-read_json(file.path(base,"docs","figure4_configuration.json"))
components <- c("Inoculum","Substrate","Interaction","Residual")
component_colours <- unlist(config$variance_partition$component_colours)
terms <- c("inoculum","substrate","inoculum:substrate","Residual")
idx <- match(terms,stats$term)
stopifnot(!anyNA(idx),identical(names(component_colours),components))
partition <- data.frame(component=factor(components,levels=components),term=terms,
  df=stats$df[idx],inertia=stats$SumOfSqs[idx],proportion=stats$R2[idx],percent=100*stats$R2[idx])
partition$xmax <- cumsum(partition$proportion)
partition$xmin <- c(0,head(partition$xmax,-1))
partition$x <- (partition$xmin+partition$xmax)/2
partition$label <- sprintf("%.1f",partition$percent)
annotations <- data.frame(factor=components[1:3],term=terms[1:3],p=stats$p[idx[1:3]])
annotations$label <- sprintf("%s: p = %.4f",annotations$factor,annotations$p)
annotations$x <- -9.25;annotations$y <- c(7.65,7.05,6.45)
stopifnot(abs(sum(partition$proportion)-1)<1e-12,
  max(abs(annotations$p-c(.0064,.0008,.1848)))<1e-12,
  all(partition$xmax-partition$xmin>0),
  max(abs((partition$xmax-partition$xmin)-partition$proportion))<1e-12)

# Same layout, ordering, and component colors as the amplicon Figure 4 inset.
legend_data <- data.frame(component=factor(components,levels=components),
  x=c(.045,.55,.045,.55),y=c(.30,.30,.10,.10))
inset <- ggplot()+
  annotate("rect",xmin=.002,xmax=.998,ymin=.002,ymax=.998,fill="white",colour="#BDBDBD",linewidth=.3)+
  annotate("text",x=.025,y=.91,label="Variance explained (%)",hjust=0,size=12.5/ggplot2::.pt,fontface="bold")+
  geom_rect(data=partition,aes(xmin=.025+.95*xmin,xmax=.025+.95*xmax,ymin=.51,ymax=.73,fill=component),
    colour="white",linewidth=.35)+
  geom_text(data=partition[partition$component!="Residual",],aes(x=.025+.95*x,y=.62,label=label),
    colour="white",size=10.5/ggplot2::.pt)+
  geom_text(data=partition[partition$component=="Residual",],aes(x=.025+.95*x,y=.62,label=label),
    colour="#222222",size=10.5/ggplot2::.pt)+
  geom_rect(data=legend_data,aes(xmin=x,xmax=x+.055,ymin=y-.045,ymax=y+.045,fill=component),colour=NA)+
  geom_text(data=legend_data,aes(x=x+.077,y=y,label=component),hjust=0,size=11.5/ggplot2::.pt)+
  scale_fill_manual(values=component_colours,guide="none")+
  coord_cartesian(xlim=c(0,1),ylim=c(0,1),expand=FALSE,clip="off")+theme_void()+
  theme(plot.margin=margin(0,0,0,0),plot.background=element_rect(fill="white",colour=NA))
fig <- fig+geom_text(data=annotations,aes(x=x,y=y,label=label),inherit.aes=FALSE,
  hjust=0,vjust=.5,size=12.5/ggplot2::.pt,colour="#333333",show.legend=FALSE)
built <- ggplot_build(fig)
stopifnot(is.null(fig$labels$title),is.null(fig$labels$subtitle),length(fig$layers)==7,
  identical(original$data,built$data[1:6]),
  nrow(built$data[[3]])==41,nrow(built$data[[5]])==12,
  all(built$data[[5]]$colour=="#3F3F3F"),all(built$data[[6]]$colour=="#3F3F3F"),
  all(built$data[[6]]$size==13/ggplot2::.pt),all(built$data[[5]]$alpha==1),
  sum(!is.na(built$data[[6]]$segment.colour))==3,
  identical(fig$scales$get_scales("colour")$guide,"none"),
  identical(as.character(built$plot$scales$get_scales("fill")$get_labels()),legend_labels),
  identical(as.character(built$plot$scales$get_scales("shape")$get_labels()),legend_labels))
inset_built <- ggplot_build(inset)
stopifnot(nrow(inset_built$data[[3]])==4,
  max(abs((inset_built$data[[3]]$xmax-inset_built$data[[3]]$xmin)/.95-partition$proportion))<1e-12)
figure <- fig+inset_element(inset,left=.69,bottom=.76,right=.985,top=.985,
  align_to="panel",on_top=TRUE,clip=TRUE)
for(ext in c("pdf","png"))ggsave(file.path(out,paste0("pcoa_publication.",ext)),figure,
  width=11.6,height=9.3,dpi=300,bg="white")
ggsave(file.path(out,"pcoa_preview.png"),figure,width=11.6,height=9.3,
  dpi=180*7.2/11.6,bg="white")
write.csv(partition,file.path(data_out,"sequential_variance_partition.csv"),row.names=FALSE)
write.csv(annotations,file.path(data_out,"displayed_factor_p_values.csv"),row.names=FALSE)
write.csv(stats,file.path(data_out,"sequential_permanova.csv"),row.names=FALSE,na="")
write.csv(data.frame(component=names(component_colours),colour=unname(component_colours)),
  file.path(data_out,"variance_component_colours.csv"),row.names=FALSE)
saveRDS(figure,file.path(data_out,"fig5v3_plot.rds"))
writeLines(capture.output(sessionInfo()),file.path(out,"R_session_info.txt"))
write_json(list(N=41,plotted_functional_groups=12,arrow_colour="#3F3F3F",loading_label_size_pt=13,
  functional_process_categories=FALSE,functional_category_legend=FALSE,
  coordinates_loadings_points_centroids_and_ellipses_unchanged=TRUE,
  model="Saved sequential PERMANOVA; inoculum, substrate, interaction; block restricts permutations",
  variance_bar_order=components,variance_bar_percentages=partition$percent,
  factor_p_values=annotations$p,inset_style="Matches accepted amplicon Figure 4",
  substrate_and_inoculum_legend_labels=legend_labels,
  title_present=FALSE),file.path(out,"validation.json"),pretty=TRUE,auto_unbox=TRUE)
cat("Rendered microcosm functional ordination.\n")
