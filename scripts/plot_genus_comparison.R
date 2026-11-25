suppressPackageStartupMessages({
  library(phyloseq); library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(jsonlite)
})
script_file<-normalizePath(sub("^--file=","",commandArgs(trailingOnly=FALSE)[grep("^--file=",commandArgs(trailingOnly=FALSE))]))
root<-dirname(dirname(script_file));source_root<-root
out<-file.path(root,"results","reproduced","amplicon_comparison");dir.create(out,recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(out,"data"),showWarnings=FALSE);dir.create(file.path(out,"provenance"),showWarnings=FALSE)
stem<-"field_microcosm_genus_comparison"
data_root<-file.path(root,"data","amplicons")
for(ex in parse(file.path(root,"scripts","run_amplicon_audit.R"))) {
  if(is.call(ex)&&identical(ex[[1]],as.name("<-"))&&is.symbol(ex[[2]])&&as.character(ex[[2]])%in%
     c("split_taxonomy","load_16s","load_its","standardize_subset","clean_genus"))eval(ex)
}
write_table<-function(x,name)write.csv(x,file.path(out,"data",name),row.names=FALSE,na="")
colour<-c(Field="#0072B2",Microcosm="#D55E00")
summary_rows<-list();sample_rows<-list();taxon_rows<-list();qc<-list()
sets<-list();selections<-list();plots<-list();venns<-list();venn_counts<-list()
invalid_genus<-function(x)is.na(x)|!nzchar(trimws(x))|grepl(
  "uncultured|unidentified|unknown|unclassified|norank|ambiguous|metagenome|unassigned|incertae_sedis",x,ignore.case=TRUE)

axis_labels<-function(x)as.expression(lapply(x,function(k) {
  if(startsWith(k,"Unclassified_")) {
    family<-sub("^Unclassified_","",k)
    if(family=="family") return(bquote("Family unresolved"))
    return(bquote(paste("Uncl. ",.(family))))
  }
  if(grepl("_gen_Incertae_sedis$",k)) return(bquote(paste("Uncl. ",.(sub("_gen_Incertae_sedis$","",k)))))
  if(k=="Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium")
    return(bquote(atop(italic("Allorhizobium–Neorhizobium–"),italic("Pararhizobium–Rhizobium"))))
  if(k=="Methylobacterium-Methylorubrum") return(bquote(atop(italic("Methylobacterium–"),italic("Methylorubrum"))))
  bquote(italic(.(k)))
}))

for(marker in c("ITS","16S")) {
  ps0<-if(marker=="ITS")load_its() else load_16s()
  sets[[marker]]<-list()
  for(experiment in c("field","microcosm")) {
    prepared<-standardize_subset(ps0,marker,experiment);ps<-prepared$ps
    m<-as(otu_table(ps),"matrix");if(!taxa_are_rows(ps))m<-t(m)
    rel<-sweep(m,2,colSums(m),"/")
    tx<-as(tax_table(ps),"matrix");category<-clean_genus(ps)
    aggregated<-rowsum(rel,group=category,reorder=FALSE)
    stopifnot(max(abs(colSums(aggregated)-1))<1e-12)
    means<-rowMeans(aggregated);sds<-apply(aggregated,1,sd)
    experiment_label<-if(experiment=="field")"Field" else "Microcosm"
    sr<-data.frame(marker,experiment=experiment_label,category=names(means),mean_relative_abundance=means,
      sd_relative_abundance=sds,N_samples=nsamples(ps))
    summary_rows[[length(summary_rows)+1]]<-sr
    sample_rows[[length(sample_rows)+1]]<-as.data.frame(t(aggregated),check.names=FALSE)%>%
      tibble::rownames_to_column("sample_id")%>%pivot_longer(-sample_id,names_to="category",values_to="relative_abundance")%>%
      mutate(marker=marker,experiment=experiment_label)
    genus<-trimws(as.character(tx[,"Genus"]))
    is_named<-!invalid_genus(genus)
    sets[[marker]][[experiment_label]]<-sort(unique(genus[is_named]))
    taxon_rows[[length(taxon_rows)+1]]<-data.frame(marker,experiment=experiment_label,asv_id=rownames(tx),
      raw_genus=genus,composition_category=category,in_named_genus_venn=is_named)
    qc[[length(qc)+1]]<-data.frame(marker,experiment=experiment_label,N=nsamples(ps),ASVs_before_rarefaction=ntaxa(ps),
      minimum_prevalence_samples=max(2,ceiling(.05*nsamples(ps))),named_genus_labels=length(sets[[marker]][[experiment_label]]))
    if(experiment=="microcosm") {
      rare<-readRDS(file.path(root,"results","reproduced","amplicons","microcosm",marker,"rarefied_asv_matrix.rds"))
      stopifnot(setequal(rownames(m),colnames(rare)),setequal(colnames(m),rownames(rare)))
    }
  }
}
summaries<-bind_rows(summary_rows);qc<-bind_rows(qc)
stopifnot(identical(qc$N,c(18L,40L,19L,36L)),identical(qc$ASVs_before_rarefaction,c(416L,240L,961L,634L)))
write_table(summaries,"genus_category_means.csv");write_table(bind_rows(sample_rows),"sample_genus_relative_abundance.csv")
write_table(bind_rows(taxon_rows),"asv_taxonomy_and_venn_eligibility.csv");write_table(qc,"sample_and_filter_summary.csv")
membership_rows<-list()
for(marker in c("ITS","16S")) {
  mm<-summaries%>%filter(.data$marker==.env$marker)
  top<-mm%>%group_by(experiment)%>%arrange(desc(mean_relative_abundance),category,.by_group=TRUE)%>%slice_head(n=10)%>%ungroup()
  selected<-unique(top$category)
  wide<-mm%>%filter(category%in%selected)%>%select(category,experiment,mean_relative_abundance)%>%
    pivot_wider(names_from=experiment,values_from=mean_relative_abundance,values_fill=0)
  # Retain explicit zeroes for selected categories not detected in one experiment.
  for(k in c("Field","Microcosm")) if(!k%in%names(wide))wide[[k]]<-0
  wide<-wide%>%arrange(desc(pmax(Field,Microcosm)),category)
  selections[[marker]]<-data.frame(marker,wide)
  wide$category<-factor(wide$category,levels=rev(wide$category))
  wide$y<-as.numeric(wide$category)
  long<-wide%>%pivot_longer(c(Field,Microcosm),names_to="experiment",values_to="mean_relative_abundance")
  long$experiment<-factor(long$experiment,levels=c("Field","Microcosm"))
  long$y<-long$y+ifelse(long$experiment=="Field",.10,-.10)
  plots[[marker]]<-ggplot(wide,aes(y=y))+
    geom_segment(aes(x=100*Field,xend=100*Microcosm,y=y+.10,yend=y-.10),colour="#B6B6B6",linewidth=.65)+
    geom_point(data=long,aes(x=100*mean_relative_abundance,fill=experiment,shape=experiment),size=2.8,colour="white",stroke=.3)+
    scale_fill_manual(values=colour,name=NULL)+scale_shape_manual(values=c(Field=21,Microcosm=24),name=NULL)+
    scale_y_continuous(breaks=seq_len(nrow(wide)),labels=axis_labels(levels(wide$category)),limits=c(.5,nrow(wide)+.5),expand=expansion(mult=0))+
    scale_x_continuous(limits=c(0,32),breaks=c(0,10,20,30),expand=expansion(mult=c(.02,.03)))+
    theme_classic(base_size=12)+theme(axis.title.y=element_blank(),axis.text.y=element_text(size=10.5,colour="#222222"),
      axis.text.x=element_text(colour="#222222"),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
      panel.grid.major.x=element_line(colour="#EAEAEA",linewidth=.35),plot.title=element_text(face="bold",size=15),
      plot.margin=margin(9,15,8,8),legend.position="top",plot.title.position="plot")+
    labs(x="Mean relative abundance (%)",title=if(marker=="ITS")"A  Fungi (ITS)" else "B  Bacteria (16S)")

  field<-sets[[marker]]$Field;micro<-sets[[marker]]$Microcosm
  only_field<-setdiff(field,micro);shared<-intersect(field,micro);only_micro<-setdiff(micro,field)
  stopifnot(length(only_field)+length(shared)==length(field),length(only_micro)+length(shared)==length(micro),
            !length(intersect(only_field,only_micro)))
  all_genera<-sort(union(field,micro))
  membership_rows[[marker]]<-data.frame(marker,genus=all_genera,detected_field=all_genera%in%field,
    detected_microcosm=all_genera%in%micro,
    membership=ifelse(all_genera%in%shared,"Shared",ifelse(all_genera%in%only_field,"Field only","Microcosm only")))
  venn_counts[[marker]]<-data.frame(marker,field_only=length(only_field),shared=length(shared),microcosm_only=length(only_micro),
    field_total=length(field),microcosm_total=length(micro),union_total=length(all_genera))
  theta<-seq(0,2*pi,length.out=361)
  circles<-bind_rows(data.frame(x=-.55+cos(theta),y=sin(theta),experiment="Field"),
                     data.frame(x=.55+cos(theta),y=sin(theta),experiment="Microcosm"))
  venns[[marker]]<-ggplot(circles,aes(x,y,group=experiment))+
    geom_polygon(aes(fill=experiment),alpha=.15,colour=NA)+
    geom_path(aes(colour=experiment),linewidth=.65)+
    annotate("text",x=c(-1.02,0,1.02),y=0,label=c(length(only_field),length(shared),length(only_micro)),size=6.1,fontface="bold",colour="#222222")+
    annotate("text",x=c(-.91,.91),y=1.23,label=c(paste0("Field (",length(field),")"),paste0("Microcosm (",length(micro),")")),
      size=4.2,colour=unname(colour))+
    scale_fill_manual(values=colour,guide="none")+scale_colour_manual(values=colour,guide="none")+
    coord_equal(xlim=c(-1.75,1.75),ylim=c(-1.08,1.42),expand=FALSE)+theme_void(base_size=12)+
    theme(plot.title=element_text(face="bold",size=15,hjust=0),plot.margin=margin(15,15,6,8))+
    labs(title=if(marker=="ITS")"C  Shared fungal genera" else "D  Shared bacterial genera")
}
write_table(bind_rows(selections),"plotted_genus_category_means.csv")
write_table(bind_rows(membership_rows),"named_genus_membership.csv");write_table(bind_rows(venn_counts),"venn_counts.csv")
caption<-paste("A–B: union of the ten most abundant categories in each experiment; each sample has equal weight.",
  "Uncl.: genus unresolved within the indicated taxon. C–D: named genus-level labels only; unresolved assignments excluded.",
  "Venn areas are not proportional to counts. Presence and absence refer to detection after filtering.",sep="\n")
fig<-wrap_plots(plots$ITS,plots$`16S`,wrap_elements(full=venns$ITS),wrap_elements(full=venns$`16S`),ncol=2)+
  plot_layout(heights=c(2.4,1.05),guides="collect")+
  plot_annotation(title="Field and microcosm litter microbiomes",caption=caption,
    theme=theme(plot.title=element_text(size=18,face="bold",hjust=0),plot.caption=element_text(size=10,hjust=0),
                plot.margin=margin(12,15,12,12)))
fig<-fig & theme(legend.position="top")
ggsave(file.path(out,paste0(stem,".pdf")),fig,width=14.2,height=11.8,device=cairo_pdf)
ggsave(file.path(out,paste0(stem,".png")),fig,width=14.2,height=11.8,dpi=200)
writeLines(capture.output(sessionInfo()),file.path(out,"provenance","R_session_info.txt"))
file.copy(script_file,file.path(out,"provenance",basename(script_file)),overwrite=TRUE)
file.copy(file.path(root,"scripts","run_amplicon_audit.R"),file.path(out,"provenance","run_amplicon_audit.R"),overwrite=TRUE)
write_json(list(samples_and_filter_counts_checked=TRUE,relative_abundance_denominator="all retained ASV reads, before rarefaction; no additional genus-count filter",
  mean="equal-weight mean across retained samples, not equal-weight treatment cells",top_categories="union of top ten per experiment, ranked for display by larger mean",
  venn_definition="detected named genus-level labels after marker/experiment-specific ASV filtering; complexes count as one database-assigned label",
  venn_exclusions="empty genus labels, uncultured, unidentified, unknown, unclassified, norank, ambiguous, metagenome, unassigned, incertae_sedis",
  absence_interpretation="not detected in retained samples under these filters; not evidence of ecological exclusivity",
  plot_inference="descriptive comparison; no new statistical tests"),
  file.path(out,"figure_configuration.json"),pretty=TRUE,auto_unbox=TRUE)
cat("Saved",file.path(out,paste0(stem,".pdf")),"\n")
print(bind_rows(venn_counts))
