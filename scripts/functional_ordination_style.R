# Monochrome microcosm correlation biplot.
make_microcosm_pcoa <- function(xy,loadings,variance) {
origins<-c("Full sunlight","Reduced UV","Reduced UV and blue-green")
legend_labels<-c("Full sunlight","Reduced UV","Reduced UV and blue-green")
colours<-setNames(c("#E3AC18","#77BCE8","#555555"),origins);shapes<-setNames(c(21,24,22),origins);grey<-"#3F3F3F"
for(v in c("substrate","inoculum"))xy[[v]]<-factor(xy[[v]],levels=origins)
short_labels <- c(
  "Sensory photoreception"="Sensory\nphotoreception",
  "Carotenoid biosynthesis and modification"="Carotenoid\nbiosynthesis",
  "Enzymatic ROS defense"="Enzymatic ROS\ndefense",
  "Photoreactivation and nucleotide-excision repair"="Photoreactivation\nand NER",
  "General DNA-damage repair and tolerance"="DNA repair and\ndamage tolerance",
  "Cellulose degradation"="Cellulose\ndegradation",
  "Hemicellulose degradation"="Hemicellulose\ndegradation",
  "Pectin degradation"="Pectin\ndegradation",
  "Cutin degradation"="Cutin\ndegradation",
  "Alpha-glucan and fructan processing"="Alpha-glucan and\nfructan processing",
  "Oligosaccharide and disaccharide utilization"="Oligo-/disaccharide\nutilization",
  "Aromatic-compound catabolism"="Aromatic-compound\ncatabolism")
loadings$display_label <- unname(short_labels[loadings$functional_group])
stopifnot(!anyNA(loadings$display_label),nrow(loadings)==12,!anyDuplicated(loadings$functional_group))
# One common 13% increase preserves vector directions and relative lengths.
mul <- min(diff(range(xy$Axis1))/diff(range(loadings$PCoA1)),
           diff(range(xy$Axis2))/diff(range(loadings$PCoA2))) * 1.20 * 1.13
loadings$xend <- loadings$PCoA1*mul;loadings$yend <- loadings$PCoA2*mul
cent <- aggregate(cbind(Axis1,Axis2)~substrate,xy,mean)
cent$substrate <- factor(cent$substrate,levels=origins)
cent$display_label <- "";xy$display_label <- ""
obstacles <- rbind(data.frame(xend=xy$Axis1,yend=xy$Axis2,display_label="",functional_group=NA_character_),
  data.frame(xend=cent$Axis1,yend=cent$Axis2,display_label="",functional_group=NA_character_))
arrow_obstacles <- do.call(rbind,lapply(seq_len(nrow(loadings)),function(i) {
  f <- seq(.08,.98,length.out=24)
  data.frame(xend=f*loadings$xend[i],yend=f*loadings$yend[i],display_label="",functional_group=NA_character_)
}))
label_data <- rbind(loadings[,c("xend","yend","display_label","functional_group")],obstacles,arrow_obstacles)
# Initial anchors are retuned for the microcosm geometry; repulsion includes
# every sample and substrate centroid as a blank-label obstacle.
anchors <- data.frame(functional_group=names(short_labels),
  x=c(-6.2,-3.8,4.2,.25,-5.7,-7.8,-4.7,-7.4,-7.2,-4.1,-3.8,-6.0),
  y=c(-1.7,-2.6,3.0,6.9,2.9,1.55,-3.8,-.7,.55,5.65,2.5,-2.8))
idx <- match(label_data$functional_group,anchors$functional_group)
label_data$dx <- ifelse(is.na(idx),0,anchors$x[idx]-label_data$xend)
label_data$dy <- ifelse(is.na(idx),0,anchors$y[idx]-label_data$yend)
leaders <- c("Alpha-glucan and fructan processing","Hemicellulose degradation","Enzymatic ROS defense")
has_leader <- label_data$functional_group %in% leaders
point_sizes <- ifelse(nzchar(label_data$display_label),4.25,6.5);point_sizes[has_leader] <- 0
point_sizes[seq.int(nrow(loadings)+nrow(obstacles)+1,nrow(label_data))] <- 1.3
xlim <- c(-9.6,7.5);ylim <- c(-4.6,7.9)
fig <- ggplot(xy,aes(Axis1,Axis2))+
  stat_ellipse(aes(fill=substrate,group=substrate),geom="polygon",type="t",level=.75,
    alpha=.11,colour=NA,show.legend=FALSE)+
  stat_ellipse(aes(colour=substrate,group=substrate),type="t",level=.75,
    linetype="dotted",linewidth=.55,show.legend=FALSE)+
  geom_point(aes(fill=substrate,shape=inoculum),size=4.25,stroke=.9,colour="black")+
  geom_point(data=cent,aes(fill=substrate),shape=23,size=5.6,stroke=1,colour="black",show.legend=FALSE)+
  geom_segment(data=loadings,aes(x=0,y=0,xend=xend,yend=yend),inherit.aes=FALSE,
    colour=grey,alpha=1,linewidth=.62,arrow=arrow(length=grid::unit(.23,"cm")))+
  geom_text_repel(data=label_data,aes(x=xend,y=yend,label=display_label),inherit.aes=FALSE,
    colour=grey,alpha=1,size=13/ggplot2::.pt,fontface="bold",lineheight=.95,
    segment.colour=ifelse(has_leader,grey,NA_character_),segment.size=.2,segment.alpha=.5,
    min.segment.length=0,point.padding=grid::unit(1e-6,"lines"),point.size=point_sizes,
    box.padding=.45,max.overlaps=Inf,force=1.2,force_pull=1,
    nudge_x=label_data$dx,nudge_y=label_data$dy,max.time=Inf,max.iter=30000,seed=42,
    xlim=xlim,ylim=ylim,show.legend=FALSE)+
  scale_colour_manual(values=colours,guide="none")+
  scale_fill_manual(values=colours,breaks=origins,labels=legend_labels,name="Microcosm substrate")+
  scale_shape_manual(values=shapes,breaks=origins,labels=legend_labels,name="Inoculum origin")+
  guides(fill=guide_legend(nrow=1,order=1,title.position="left",override.aes=list(shape=21,size=4.5)),
    shape=guide_legend(nrow=1,order=2,title.position="left",override.aes=list(fill="white",size=4.5)))+
  coord_cartesian(xlim=xlim,ylim=ylim,clip="on")+
  theme_bw(base_size=13)+theme(
    panel.grid.major=element_line(colour="#EEEEEE",linewidth=.25),panel.grid.minor=element_blank(),
    axis.title=element_text(size=17),axis.text=element_text(size=14,colour="#333333"),
    legend.position="bottom",legend.box="vertical",legend.box.just="left",
    legend.title=element_text(size=14),legend.text=element_text(size=13),
    legend.key.width=grid::unit(.65,"cm"),legend.key.height=grid::unit(.45,"cm"),
    legend.box.margin=margin(3,0,3,0),legend.margin=margin(2,4,2,4),
    legend.spacing.y=grid::unit(.08,"cm"),plot.margin=margin(10,18,8,10))+
  labs(x=sprintf("PCoA1 (%.1f%%)",100*variance$variance_fraction[1]),
       y=sprintf("PCoA2 (%.1f%%)",100*variance$variance_fraction[2]))

fig
}
