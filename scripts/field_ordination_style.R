# Field functional correlation biplot style.
make_pcoa_figure <- function(scores,loadings,variance,global,experiment) {
  tvar<-if(experiment=="field")"Sunlight exposure" else "Litter photodegradation"
  ptitle<-if(experiment=="field")"Functional profiles under field sunlight exposure" else
    "Functional profiles after litter photodegradation"
  legend_labels<-setNames(levels(scores$Treatment),levels(scores$Treatment))
  xlab<-sprintf("PCoA1 (%.1f%%)",100*variance[1]);ylab<-sprintf("PCoA2 (%.1f%%)",100*variance[2])
  short_labels<-c(
    "Sensory photoreception"="Sensory\nphotoreception",
    "Carotenoid biosynthesis and modification"="Carotenoid\nbiosynthesis",
    "Catalase, superoxide and peroxide defense"="Enzymatic ROS\ndefense",
    "Photoreactivation and nucleotide-excision repair"="Photoreactivation\nand NER",
    "General DNA-damage repair and tolerance"="DNA repair and\ndamage tolerance",
    "Cellulose degradation"="Cellulose\ndegradation",
    "Hemicellulose degradation"="Hemicellulose\ndegradation",
    "Pectin degradation"="Pectin\ndegradation",
    "Cutin degradation"="Cutin\ndegradation",
    "Alpha-glucan and fructan processing"="Alpha-glucan and\nfructan processing",
    "Oligosaccharide and disaccharide utilization"="Oligo-/disaccharide\nutilization",
    "Lignin-associated oxidation"="Lignin-associated\noxidation",
    "Aromatic-compound catabolism"="Aromatic-compound\ncatabolism")
  loadings$display_label<-unname(short_labels[loadings$functional_group])
  loadings$display_label[is.na(loadings$display_label)]<-loadings$functional_group[is.na(loadings$display_label)]
  if(experiment=="field") {
    # Interpretive display categories only: no regrouping of the input data,
    # loadings, ordination, ellipses, or inferential tests.
    process_colours<-c(
      "Light sensing and stress protection"="#7B3294",
      "Plant-polymer and carbohydrate processing"="#008837",
      "Lignin and aromatic transformations"="#A6611A")
    process_groups<-c(
      "Sensory photoreception"="Light sensing and stress protection",
      "Carotenoid biosynthesis and modification"="Light sensing and stress protection",
      "Enzymatic ROS defense"="Light sensing and stress protection",
      "Photoreactivation and nucleotide-excision repair"="Light sensing and stress protection",
      "General DNA-damage repair and tolerance"="Light sensing and stress protection",
      "Cellulose degradation"="Plant-polymer and carbohydrate processing",
      "Hemicellulose degradation"="Plant-polymer and carbohydrate processing",
      "Pectin degradation"="Plant-polymer and carbohydrate processing",
      "Cutin degradation"="Plant-polymer and carbohydrate processing",
      "Alpha-glucan and fructan processing"="Plant-polymer and carbohydrate processing",
      "Oligosaccharide and disaccharide utilization"="Plant-polymer and carbohydrate processing",
      "Lignin-associated oxidation"="Lignin and aromatic transformations",
      "Aromatic-compound catabolism"="Lignin and aromatic transformations")
    stopifnot(length(process_groups)==13,!anyDuplicated(names(process_groups)),
              !anyDuplicated(loadings$functional_group),
              setequal(loadings$functional_group,names(process_groups)))
    loadings$process_category<-factor(unname(process_groups[loadings$functional_group]),
                                     levels=names(process_colours))
    loadings$display_label[loadings$functional_group=="Enzymatic ROS defense"]<-"Enzymatic ROS\ndefense"
    # The common multiplier preserves all loading directions and relative lengths.
    arrow_scale<-.84*1.13
    mul<-min(diff(range(scores$PCoA1))/max(diff(range(loadings$PCoA1)),1e-12),
             diff(range(scores$PCoA2))/max(diff(range(loadings$PCoA2)),1e-12))*arrow_scale
    loadings$xend<-loadings$PCoA1*mul;loadings$yend<-loadings$PCoA2*mul
    cent<-scores%>%group_by(Treatment)%>%summarise(PCoA1=mean(PCoA1),PCoA2=mean(PCoA2),.groups="drop")
    xr<-range(c(scores$PCoA1,loadings$xend));yr<-range(c(scores$PCoA2,loadings$yend))
    xpad<-max(diff(xr)*.23,.01);ypad<-max(diff(yr),.01)
    xlim<-c(xr[1]-xpad,xr[2]+xpad);ylim<-c(yr[1]-.13*ypad,yr[2]+.25*ypad)
    # Empty labels make samples and centroids obstacles for label repulsion.
    label_data<-bind_rows(loadings,
      transmute(scores,xend=PCoA1,yend=PCoA2,display_label="",process_category=names(process_colours)[1]),
      transmute(cent,xend=PCoA1,yend=PCoA2,display_label="",process_category=names(process_colours)[1]))
    label_nudges<-data.frame(
      functional_group=names(process_groups),
      dx=c(-.25,-.3,.7,-.9,.9,-.9,-1.4,-.8,-1.7,.7,.4,-.5,1.3),
      dy=c(1.0,1.1,2.5,1.8,.8,-.9,1.0,2.2,-.2,1.2,-.8,-1.0,-.8))
    label_data<-left_join(label_data,label_nudges,by="functional_group")
    label_data$dx[is.na(label_data$dx)]<-0;label_data$dy[is.na(label_data$dy)]<-0
    main<-ggplot(scores,aes(PCoA1,PCoA2,fill=Treatment))+
      stat_ellipse(geom="polygon",type="t",level=.75,alpha=.11,colour=NA,show.legend=FALSE)+
      stat_ellipse(aes(colour=Treatment),type="t",level=.75,linetype="dotted",linewidth=.55,show.legend=FALSE)+
      geom_point(shape=21,size=4.25,stroke=.9,colour="black")+
      geom_point(data=cent,shape=23,size=5.6,stroke=1,colour="black",show.legend=FALSE)+
      geom_segment(data=loadings,aes(x=0,y=0,xend=xend,yend=yend,colour=process_category),
                   inherit.aes=FALSE,arrow=arrow(length=unit(.23,"cm")),alpha=1,linewidth=.62)+
      geom_text_repel(data=label_data,aes(x=xend,y=yend,label=display_label,colour=process_category),
                      inherit.aes=FALSE,show.legend=FALSE,size=13/ggplot2::.pt,fontface="bold",
                      alpha=1,segment.size=.2,segment.alpha=.5,min.segment.length=.2,
                      box.padding=.45,point.padding=.65,
                      point.size=4.25,max.overlaps=Inf,force=1.2,force_pull=1,
                      nudge_x=label_data$dx,nudge_y=label_data$dy,
                      max.time=Inf,max.iter=30000,seed=42,xlim=xlim,ylim=ylim)+
      scale_colour_manual(values=c(colours,process_colours),breaks=names(process_colours),
                          name="Functional-process category",na.translate=FALSE)+
      scale_fill_manual(values=colours,breaks=levels(scores$Treatment),labels=legend_labels,
                        name=tvar)+
      guides(fill=guide_legend(position="top",nrow=1,order=1,
                              override.aes=list(size=4.5)),
             colour=guide_legend(position="bottom",ncol=1,order=2,
                                title.position="left",override.aes=list(linewidth=.8,alpha=1)))+
      coord_cartesian(xlim=xlim,ylim=ylim,clip="off")+
      theme_bw(base_size=13)+theme(
        panel.grid.major=element_line(colour="#EEEEEE",linewidth=.25),
        panel.grid.minor=element_blank(),
        axis.title=element_text(size=17),axis.text=element_text(size=14,colour="#333333"),
        legend.title=element_text(size=14),legend.text=element_text(size=13),
        legend.key.width=unit(.65,"cm"),legend.key.height=unit(.45,"cm"),
        legend.box.margin=margin(3,0,3,0),legend.margin=margin(2,4,2,4),
        legend.spacing.y=unit(.08,"cm"),plot.margin=margin(10,18,8,10))+
      labs(x=xlab,y=ylab)
    return(main)
  }
}
