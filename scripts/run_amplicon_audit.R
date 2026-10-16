suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(permute)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

set.seed(42)
script_arg <- commandArgs(trailingOnly=FALSE)[grep("^--file=",commandArgs(trailingOnly=FALSE))]
repo_root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."), mustWork=TRUE)
data_root <- file.path(repo_root, "data", "amplicons")
out_root <- file.path(repo_root, "results", "reproduced", "amplicons")
dir.create(out_root, recursive=TRUE, showWarnings=FALSE)

colours <- c("S+"="#E3AC18", "UV-"="#77BCE8", "S-"="#555555",
             "SusS+"="#E3AC18", "SusUV-"="#77BCE8", "SusS-"="#555555")
labels <- c("S+"="Full sunlight", "UV-"="Reduced UV", "S-"="Reduced UV and blue-green",
            "SusS+"="Full sunlight", "SusUV-"="Reduced UV", "SusS-"="Reduced UV and blue-green")

split_taxonomy <- function(strings, ranks) {
  vals <- strsplit(as.character(strings), ";\\s*")
  vals <- lapply(vals, function(z) {length(z) <- length(ranks); z})
  m <- do.call(rbind, vals)
  colnames(m) <- ranks
  m <- apply(m, 2, function(z) sub("^[A-Za-z0-9_]+__", "", z))
  m[is.na(m)] <- ""
  m
}

load_16s <- function() {
  ft <- read.csv(file.path(data_root,"16S_feature_table.csv"), check.names=FALSE,
                 stringsAsFactors=FALSE)
  ids <- trimws(as.character(ft[[1]]))
  mat <- as.matrix(ft[,-1,drop=FALSE]); storage.mode(mat) <- "numeric"; rownames(mat) <- ids
  md <- read.csv(file.path(data_root,"16S_sample_metadata.csv"), check.names=FALSE,
                 stringsAsFactors=FALSE)
  rownames(md) <- md$sample_id
  md <- md[colnames(mat),,drop=FALSE]
  tx <- read.delim(file.path(data_root,"16S_taxonomy.tsv"), check.names=FALSE,
                   stringsAsFactors=FALSE)
  tx <- tx[match(rownames(mat),tx[[1]]),,drop=FALSE]
  tm <- split_taxonomy(tx$Taxon,c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))
  rownames(tm) <- rownames(mat)
  ps <- phyloseq(otu_table(mat,taxa_are_rows=TRUE),sample_data(md),tax_table(tm))
  ps <- subset_taxa(ps, grepl("Bacteria",Kingdom,ignore.case=TRUE))
  bad <- apply(as(tax_table(ps),"matrix"),1,function(z) any(grepl("chloroplast|mitochond",z,ignore.case=TRUE)))
  ps <- prune_taxa(!bad,ps)
  prune_taxa(taxa_sums(ps)>0,prune_samples(sample_sums(ps)>0,ps))
}

load_its <- function() {
  ft <- read.delim(file.path(data_root,"ITS_feature_table.tsv"),header=TRUE,
                   check.names=FALSE,comment.char="",stringsAsFactors=FALSE)
  ids <- as.character(ft[[1]])
  mat <- as.matrix(ft[,-1,drop=FALSE]); storage.mode(mat) <- "numeric"; rownames(mat) <- ids
  md <- read.delim(file.path(data_root,"ITS_sample_metadata.tsv"),header=TRUE,
                   check.names=FALSE,quote="",comment.char="",stringsAsFactors=FALSE)
  rownames(md) <- md[["sample-id"]]
  md <- md[colnames(mat),,drop=FALSE]
  tx <- read.delim(file.path(data_root,"ITS_taxonomy.tsv"),check.names=FALSE,stringsAsFactors=FALSE)
  tx <- tx[match(rownames(mat),tx[[1]]),,drop=FALSE]
  tm <- split_taxonomy(tx$Taxon,c("Kingdom","Phylum","Class","Order","Family","Genus","Species","SH"))
  rownames(tm) <- rownames(mat)
  ps <- phyloseq(otu_table(mat,taxa_are_rows=TRUE),sample_data(md),tax_table(tm))
  ps <- subset_taxa(ps, grepl("Fungi",Kingdom,ignore.case=TRUE))
  prune_taxa(taxa_sums(ps)>0,prune_samples(sample_sums(ps)>0,ps))
}

standardize_subset <- function(ps, marker, experiment) {
  md <- data.frame(sample_data(ps),check.names=FALSE)
  if(marker=="16S") {
    expv <- as.character(md$Experiment); ta <- as.character(md$TreatmentA)
    tb <- as.character(md$TreatmentB); repv <- as.character(md$Replicate)
  } else {
    expv <- as.character(md$experiment); ta <- as.character(md$treatA)
    tb <- as.character(md$treatB); repv <- as.character(md$replicate)
  }
  if(experiment=="field") {
    keep <- expv=="Inoculos" & tb=="Maize"
    treatment <- ta
    inoculum <- NA_character_
  } else {
    keep <- expv=="RecyprocalT" & !grepl("Initial",ta,ignore.case=TRUE)
    treatment <- ta
    inoculum <- tb
  }
  ps <- prune_samples(keep,ps)
  ps <- prune_taxa(taxa_sums(ps)>0,ps)
  md <- data.frame(sample_data(ps),check.names=FALSE)
  if(marker=="16S") {
    md$Treatment <- as.character(md$TreatmentA); md$Block <- as.character(md$Replicate)
    md$Inoculum2 <- if(experiment=="microcosm") as.character(md$TreatmentB) else NA_character_
  } else {
    md$Treatment <- as.character(md$treatA); md$Block <- as.character(md$replicate)
    md$Inoculum2 <- if(experiment=="microcosm") as.character(md$treatB) else NA_character_
  }
  lev <- if(experiment=="field") c("S+","UV-","S-") else c("SusS+","SusUV-","SusS-")
  md$Treatment <- factor(md$Treatment,levels=lev)
  md$Block <- factor(md$Block)
  md$Inoculum2 <- factor(md$Inoculum2)
  md$PermutationBlock <- if(experiment=="field") md$Block else interaction(md$Inoculum2,md$Block,drop=TRUE)
  rownames(md) <- sample_names(ps); sample_data(ps) <- sample_data(md)
  # A marker-independent, predeclared technical threshold avoids retaining failed
  # libraries while preventing ordination-defined sample removal.
  depth_before_prevalence <- sample_sums(ps)
  # P030 is the field shotgun identifier for the same biological sample represented
  # as 16s_109 (16S) and S109 (ITS). The investigator confirmed an extraction
  # failure, so this predeclared technical exclusion is applied before prevalence
  # filtering, rarefaction, ordination, PERMANOVA, and differential abundance.
  technical_outlier <- experiment=="field" & names(depth_before_prevalence) %in%
    if(marker=="16S") "16s_109" else "S109"
  included <- depth_before_prevalence >= 15000 & !technical_outlier
  exclusion_reason <- ifelse(technical_outlier,
    "confirmed technical extraction failure (same biological sample as shotgun P030)",
    ifelse(depth_before_prevalence < 15000,
      "<15,000 target-marker reads after taxonomy filtering",""))
  ledger <- data.frame(sample_id=names(depth_before_prevalence),experiment=experiment,marker=marker,
                       treatment=as.character(md$Treatment),inoculum=as.character(md$Inoculum2),
                       block=as.character(md$Block),reads_after_taxonomy=as.numeric(depth_before_prevalence),
                       included_primary=as.logical(included),
                       exclusion_reason=exclusion_reason)
  ps <- prune_samples(included,ps); ps <- prune_taxa(taxa_sums(ps)>0,ps)
  # The declared prevalence filter is calculated after experiment/marker subsetting.
  m <- as(otu_table(ps),"matrix"); if(!taxa_are_rows(ps)) m <- t(m)
  keep_tax <- rowSums(m>0) >= ceiling(0.05*ncol(m))
  ps <- prune_taxa(keep_tax,ps); ps <- prune_taxa(taxa_sums(ps)>0,ps)
  list(ps=ps,ledger=ledger)
}

pairwise_permanova <- function(mat,md,nperm=9999) {
  lev <- levels(md$Treatment); ans <- list(); z <- 1
  for(i in seq_len(length(lev)-1)) for(j in (i+1):length(lev)) {
    keep <- md$Treatment %in% lev[c(i,j)]
    d <- vegdist(mat[keep,,drop=FALSE],method="bray")
    mm <- droplevels(md[keep,,drop=FALSE])
    set.seed(2000+i*10+j)
    fit <- adonis2(d~Treatment,data=mm,permutations=nperm,strata=mm$PermutationBlock)
    ans[[z]] <- data.frame(group1=lev[i],group2=lev[j],F=fit$F[1],R2=fit$R2[1],p=fit$`Pr(>F)`[1]);z<-z+1
  }
  out <- do.call(rbind,ans); out$p_BH <- p.adjust(out$p,"BH"); out
}

clean_genus <- function(ps) {
  tx <- data.frame(tax_table(ps),check.names=FALSE,stringsAsFactors=FALSE)
  g <- trimws(as.character(tx$Genus)); f <- trimws(as.character(tx$Family))
  g[is.na(g)|g==""] <- paste0("Unclassified_",ifelse(is.na(f)|f=="","family",f))[is.na(g)|g==""]
  g
}

export_genus <- function(ps,marker,experiment) {
  genus <- clean_genus(ps)
  m <- as(otu_table(ps),"matrix"); if(!taxa_are_rows(ps))m<-t(m)
  agg <- rowsum(m,group=genus,reorder=FALSE)
  agg <- agg[rowSums(agg)>=10,,drop=FALSE]
  counts <- as.data.frame(t(agg),check.names=FALSE); counts$sample_id <- rownames(counts)
  counts <- counts[,c("sample_id",setdiff(names(counts),"sample_id"))]
  md <- data.frame(sample_data(ps),check.names=FALSE)
  md$sample_id <- rownames(md)
  write.csv(counts,file.path(out_root,paste0(experiment,"_",marker,"_all_genus_counts.csv")),row.names=FALSE)
  write.csv(md,file.path(out_root,paste0(experiment,"_",marker,"_metadata.csv")),row.names=FALSE)

  rel <- sweep(agg,2,colSums(agg),"/")
  top <- names(sort(rowMeans(rel),decreasing=TRUE))[seq_len(min(12,nrow(rel)))]
  long <- as.data.frame(t(rel),check.names=FALSE) %>% tibble::rownames_to_column("sample_id") %>%
    pivot_longer(-sample_id,names_to="Genus",values_to="Relative_abundance") %>%
    mutate(Display=ifelse(Genus %in% top,Genus,"Other")) %>%
    group_by(sample_id,Display) %>% summarise(Relative_abundance=sum(Relative_abundance),.groups="drop") %>%
    left_join(md[,c("sample_id","Treatment")],by="sample_id")
  write.csv(long,file.path(out_root,paste0(experiment,"_",marker,"_top12_composition_data.csv")),row.names=FALSE)
  p <- ggplot(long,aes(x=sample_id,y=Relative_abundance,fill=Display))+
    geom_col(width=.92)+facet_grid(~Treatment,scales="free_x",space="free_x")+
    scale_y_continuous(labels=scales::percent)+theme_bw(base_size=10)+
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.grid=element_blank(),
          legend.position="right")+labs(x=NULL,y="Relative abundance",fill="Genus",title=marker)
  p
}

analyse_one <- function(ps,marker,experiment) {
  prepared <- standardize_subset(ps,marker,experiment)
  ps <- prepared$ps
  md <- data.frame(sample_data(ps),check.names=FALSE)
  depth <- sample_sums(ps)
  ledger <- prepared$ledger
  ledger$reads_after_prevalence <- NA_real_
  ledger$reads_after_prevalence[match(names(depth),ledger$sample_id)] <- as.numeric(depth)
  write.csv(ledger,file.path(out_root,paste0(experiment,"_",marker,"_sample_ledger.csv")),row.names=FALSE)
  rare_depth <- min(depth)
  rare <- rarefy_even_depth(ps,sample.size=rare_depth,rngseed=42,replace=FALSE,verbose=FALSE)
  mat <- as(otu_table(rare),"matrix"); if(taxa_are_rows(rare))mat<-t(mat)
  md <- data.frame(sample_data(rare),check.names=FALSE)
  set.seed(42); ord <- metaMDS(mat,distance="bray",k=2,trymax=200,autotransform=FALSE,trace=0)
  scores <- as.data.frame(scores(ord,"sites")); scores$sample_id<-rownames(scores)
  scores <- cbind(scores,md[rownames(scores),c("Treatment","Inoculum2","Block","PermutationBlock"),drop=FALSE])
  d <- vegdist(mat,"bray")
  seed_base<-if(experiment=="field")if(marker=="16S")5101 else 5102 else if(marker=="16S")5201 else 5202
  set.seed(seed_base)
  fit <- adonis2(d~Treatment,data=md,permutations=9999,strata=md$PermutationBlock)
  ctrl <- how(nperm=9999); setBlocks(ctrl)<-md$PermutationBlock
  set.seed(seed_base+100)
  disp <- betadisper(d,md$Treatment); dp <- permutest(disp,permutations=ctrl)
  pair <- pairwise_permanova(mat,md)
  global <- data.frame(experiment=experiment,marker=marker,N=nrow(md),taxa=ncol(mat),
                       rarefaction_depth=rare_depth,NMDS_stress=ord$stress,F=fit$F[1],R2=fit$R2[1],
                       p=fit$`Pr(>F)`[1],dispersion_F=dp$tab[1,"F"],dispersion_p=dp$tab[1,"Pr(>F)"])
  write.csv(global,file.path(out_root,paste0(experiment,"_",marker,"_global_statistics.csv")),row.names=FALSE)
  write.csv(pair,file.path(out_root,paste0(experiment,"_",marker,"_pairwise_permanova.csv")),row.names=FALSE)
  write.csv(scores,file.path(out_root,paste0(experiment,"_",marker,"_nmds_scores.csv")),row.names=FALSE)
  annotation <- sprintf("Blocked PERMANOVA: F = %.2f, R² = %.3f, p = %s",
                        fit$F[1],fit$R2[1],format.pval(fit$`Pr(>F)`[1],digits=2,eps=.0001))
  xr<-range(scores$NMDS1);yr<-range(scores$NMDS2)
  p <- ggplot(scores,aes(NMDS1,NMDS2,colour=Treatment,fill=Treatment))+
    stat_ellipse(geom="polygon",type="norm",level=.75,alpha=.10,colour=NA)+
    stat_ellipse(type="norm",level=.75,linetype="dotted",linewidth=.55)+
    geom_point(shape=21,size=3.2,stroke=.45,colour="black")+
    scale_fill_manual(values=colours,labels=labels)+
    scale_colour_manual(values=colours,labels=labels,guide="none")+
    annotate("label",x=xr[2]-.02*diff(xr),y=yr[2]-.02*diff(yr),label=annotation,
             hjust=1,vjust=1,size=3.05,linewidth=.18)+
    coord_cartesian(clip="off")+theme_bw(base_size=12)+
    theme(panel.grid=element_blank(),legend.position="bottom",plot.margin=margin(6,5,4,5))+
    labs(title=marker,subtitle=sprintf("N = %d; stress = %.3f",nrow(md),ord$stress),
         colour=if(experiment=="field")"Sunlight exposure" else "Litter photodegradation",fill=NULL)
  g <- export_genus(ps,marker,experiment)
  list(nmds=p,composition=g,global=global)
}

ps16 <- load_16s(); psits <- load_its()
results <- list()
for(exp in c("field","microcosm")) {
  results[[paste(exp,"16S")]] <- analyse_one(ps16,"16S",exp)
  results[[paste(exp,"ITS")]] <- analyse_one(psits,"ITS",exp)
  pn <- results[[paste(exp,"ITS")]]$nmds + results[[paste(exp,"16S")]]$nmds +
    plot_layout(guides="collect") & theme(legend.position="bottom")
  pc <- results[[paste(exp,"ITS")]]$composition + results[[paste(exp,"16S")]]$composition +
    plot_annotation(title=if(exp=="field")"Dominant genera under field sunlight exposure" else
      "Dominant genera after litter photodegradation",
      tag_levels="A",
      theme=theme(plot.title=element_text(size=16,face="bold",hjust=.5),
                  plot.tag=element_text(size=16,face="bold"),
                  plot.margin=margin(14,5,5,5)))
  ggsave(file.path(out_root,paste0(exp,"_amplicon_NMDS.pdf")),pn,width=12,height=5.7)
  ggsave(file.path(out_root,paste0(exp,"_amplicon_NMDS.png")),pn,width=12,height=5.7,dpi=300)
  ggsave(file.path(out_root,paste0(exp,"_top12_genera.pdf")),pc,width=14,height=6)
  ggsave(file.path(out_root,paste0(exp,"_top12_genera.png")),pc,width=14,height=6,dpi=300)
}
all_global <- bind_rows(lapply(results,`[[`,"global"))
write.csv(all_global,file.path(out_root,"amplicon_global_statistics.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out_root,"R_session_info.txt"))
