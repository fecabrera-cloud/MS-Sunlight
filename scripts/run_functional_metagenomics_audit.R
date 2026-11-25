suppressPackageStartupMessages({
  library(vegan); library(permute); library(emmeans);library(lmerTest)
  library(ggplot2); library(ggrepel); library(dplyr); library(tidyr); library(patchwork)
})
set.seed(42)
script_arg <- commandArgs(trailingOnly=FALSE)[grep("^--file=",commandArgs(trailingOnly=FALSE))]
repo_root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."), mustWork=TRUE)
data_root <- file.path(repo_root, "data", "functional_metagenomics")
out_root <- file.path(repo_root, "results", "reproduced", "functional_metagenomics")
dir.create(out_root,recursive=TRUE,showWarnings=FALSE)

field_matrix_path <- file.path(data_root,"field_KO_relative_abundance.tsv")
field_meta_path <- file.path(data_root,"field_metadata.csv")
micro_matrix_path <- file.path(data_root,"microcosm_KO_relative_abundance.tsv")
micro_meta_path <- file.path(data_root,"microcosm_metadata.csv")
mapping_path <- file.path(data_root,"frozen_332_KO_functional_group_mapping.csv")
definitions_path <- file.path(data_root,"functional_group_definitions.csv")

mapping <- read.csv(mapping_path,check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(nrow(mapping)==332,length(unique(mapping$functional_group))==13,!anyDuplicated(mapping$ko_id))
definitions <- read.csv(definitions_path,check.names=FALSE,stringsAsFactors=FALSE)
frozen_universe<-read.csv(file.path(data_root,"full_7931_KO_universe_review.csv"),
                          stringsAsFactors=FALSE)
stopifnot(nrow(frozen_universe)==7931,!anyDuplicated(frozen_universe$ko_id))

colours <- c("S+"="#E3AC18","UV-"="#77BCE8","S-"="#555555",
             "SusS+"="#E3AC18","SusUV-"="#77BCE8","SusS-"="#555555")
treat_labels <- c("S+"="Full sunlight","UV-"="Reduced UV","S-"="Reduced UV and blue-green",
                  "SusS+"="Full sunlight","SusUV-"="Reduced UV","SusS-"="Reduced UV and blue-green")

pairwise_permanova <- function(std,group,block) {
  lev<-levels(group); ans<-list();z<-1
  for(i in seq_len(length(lev)-1)) for(j in (i+1):length(lev)) {
    keep<-group%in%lev[c(i,j)]; mm<-data.frame(Treatment=droplevels(group[keep]),Block=droplevels(block[keep]))
    d<-dist(std[keep,,drop=FALSE],method="euclidean")
    set.seed(1000+i*10+j)
    fit<-adonis2(d~Treatment,data=mm,permutations=9999,strata=mm$Block)
    ans[[z]]<-data.frame(group1=lev[i],group2=lev[j],F=fit$F[1],R2=fit$R2[1],p=fit$`Pr(>F)`[1]);z<-z+1
  }
  x<-bind_rows(ans);x$p_BH<-p.adjust(x$p,"BH");x
}

aggregate_groups <- function(ko) {
  kos<-intersect(mapping$ko_id,colnames(ko)); submap<-mapping[match(kos,mapping$ko_id),]
  out<-sapply(unique(mapping$functional_group),function(g){
    z<-kos[submap$functional_group==g]
    if(length(z)==0) rep(0,nrow(ko)) else rowSums(ko[,z,drop=FALSE])
  })
  out<-as.matrix(out);rownames(out)<-rownames(ko);out
}

source(file.path(repo_root,"scripts","factorial_permanova.R"))
source(file.path(repo_root,"scripts","fit_functional_groups.R"))

field_group_models <- function(raw,md,prefix) {
  tests<-pairs_out<-diag<-list()
  for(g in colnames(raw)) {
    dat<-md;dat$y<-raw[,g];fit<-lm(y~replicate+Treatment,dat);a<-drop1(fit,test="F")
    tests[[g]]<-data.frame(functional_group=g,F=a["Treatment","F value"],p=a["Treatment","Pr(>F)"])
    pp<-as.data.frame(pairs(emmeans(fit,~Treatment),adjust="tukey"));pp$functional_group<-g
    names(pp)[names(pp)=="p.value"]<-"p_Tukey";pairs_out[[g]]<-pp
    r<-rstandard(fit);fv<-fitted(fit);aux<-lm(I(r^2)~fv)
    diag[[g]]<-data.frame(functional_group=g,shapiro_p=shapiro.test(r)$p.value,bp_fitted_p=pchisq(length(r)*summary(aux)$r.squared,1,lower.tail=FALSE))
  }
  a<-bind_rows(tests);a$p_BH<-p.adjust(a$p,"BH")
  write.csv(a,file.path(prefix,"functional_group_omnibus.csv"),row.names=FALSE)
  write.csv(bind_rows(pairs_out),file.path(prefix,"functional_group_pairwise.csv"),row.names=FALSE)
  write.csv(bind_rows(diag),file.path(prefix,"model_diagnostics.csv"),row.names=FALSE)
  long<-data.frame(sample_id=rownames(raw),raw,check.names=FALSE)%>%pivot_longer(-sample_id,names_to="functional_group",values_to="relative_abundance")%>%
    left_join(data.frame(sample_id=rownames(md),Treatment=md$Treatment),by="sample_id")%>%left_join(a,by="functional_group")%>%
    mutate(label=paste0(stringr::str_wrap(functional_group,30),"\nq = ",format.pval(p_BH,digits=2,eps=.001)))
  set.seed(42)
  p<-ggplot(long,aes(Treatment,relative_abundance,fill=Treatment))+geom_boxplot(width=.58,outlier.shape=NA,alpha=.65)+
    geom_jitter(width=.10,shape=21,size=1.75,colour="black")+facet_wrap(~label,scales="free_y",ncol=3)+
    scale_fill_manual(values=colours,labels=treat_labels)+theme_bw(base_size=9)+
    theme(legend.position="top",axis.text.x=element_text(size=9,angle=35,hjust=1),strip.text=element_text(size=10,face="bold"),panel.grid=element_blank())+
    labs(x=NULL,y="Relative abundance",fill="Sunlight exposure")
  for(ext in c("pdf","png"))ggsave(file.path(prefix,paste0("supplementary_functional_groups.",ext)),p,width=11,height=12,dpi=300)
}
run_one <- function(experiment,exclude_samples=character()) {
  if(experiment=="field") {
    ko0<-read.delim(field_matrix_path,check.names=FALSE);rownames(ko0)<-ko0$id;ko0$id<-NULL
    md0<-read.csv(field_meta_path,stringsAsFactors=FALSE)
    md0<-md0[md0$TreatmentB=="Maize",];rownames(md0)<-md0$id
    ids_all<-intersect(rownames(ko0),rownames(md0));ids<-setdiff(ids_all,exclude_samples)
    ko<-as.matrix(ko0[ids,,drop=FALSE]);md<-md0[ids,,drop=FALSE]
    md$Treatment<-factor(md$TreatmentA,levels=c("S+","UV-","S-"));md$replicate<-factor(md$replicate)
    md$PermutationBlock<-md$replicate
  } else {
    ko0<-read.delim(micro_matrix_path,check.names=FALSE);rownames(ko0)<-ko0$sample_id;ko0$sample_id<-NULL
    md0<-read.csv(micro_meta_path,stringsAsFactors=FALSE);rownames(md0)<-md0$sample_id
    md0<-md0[md0$Substrate!="Initial",]
    ids_all<-intersect(rownames(ko0),rownames(md0));ids<-setdiff(ids_all,exclude_samples)
    ko<-as.matrix(ko0[ids,,drop=FALSE]);md<-md0[ids,,drop=FALSE]
    md$Treatment<-factor(md$Substrate,levels=c("SusS+","SusUV-","SusS-"));md$Inoculum<-factor(md$Inoculum);md$replicate<-factor(md$replicate)
    md$PermutationBlock<-interaction(md$Inoculum,md$replicate,drop=TRUE)
  }

  storage.mode(ko)<-"numeric"
  ko<-ko[,intersect(frozen_universe$ko_id,colnames(ko)),drop=FALSE]
  raw_all<-aggregate_groups(ko);eligible<-apply(raw_all,2,sd)>0;raw<-raw_all[,eligible,drop=FALSE]
  std<-scale(raw);d<-dist(std,"euclidean")
  prefix<-file.path(out_root,experiment);dir.create(prefix,recursive=TRUE,showWarnings=FALSE)
  write_table<-function(x,n)write.csv(x,file.path(prefix,n),row.names=FALSE,na="")
  write_table(data.frame(sample_id=rownames(raw_all),raw_all,check.names=FALSE),"functional_group_relative_abundance.csv")
  write_table(data.frame(sample_id=rownames(std),std,check.names=FALSE),"functional_group_zscores.csv")
  eligibility<-data.frame(functional_group=colnames(raw_all),mapped_KOs=vapply(colnames(raw_all),function(g)sum(mapping$functional_group==g),integer(1)),
    sd=apply(raw_all,2,sd),included=eligible)
  write_table(eligibility,"group_eligibility.csv")
  if(experiment=="field") {
    set.seed(4201);fit<-adonis2(d~Treatment,data=md,permutations=9999,strata=md$PermutationBlock)
    ctrl<-how(nperm=9999);setBlocks(ctrl)<-md$PermutationBlock;set.seed(4301)
    bp<-permutest(betadisper(d,md$Treatment),permutations=ctrl)
    F<-fit$F[1];R2<-fit$R2[1];pv<-fit$`Pr(>F)`[1];dispF<-bp$tab$F[1];dispP<-bp$tab$`Pr(>F)`[1]
    write_table(pairwise_permanova(std,md$Treatment,md$PermutationBlock),"pairwise_permanova.csv")
    field_group_models(raw,md,prefix)
  } else {
    md2<-data.frame(sample_id=rownames(md),inoculum=full_treatment(md$Inoculum),substrate=full_treatment(md$Treatment),block=factor(md$replicate))
    P<-readRDS(file.path(repo_root,"data","ordination_reference","microcosm_functional_block_permutations.rds"))
    ff<-factorial_permanova(d,md2,P,prefix,4202)
    a<-ff$anova[ff$anova$term=="substrate",];b<-ff$dispersion[ff$dispersion$grouping=="substrate",]
    F<-a$F;R2<-a$R2;pv<-a$p;dispF<-b$F;dispP<-b$p
    fit_functional_groups(raw,ff$metadata,prefix)
    stopifnot(system2("Rscript",c(shQuote(file.path(repo_root,"scripts","render_functional_groups.R")),shQuote(prefix),shQuote(prefix)))==0)
  }
  pc<-cmdscale(d,k=2,eig=TRUE,add=TRUE);var<-pc$eig[pc$eig>0]/sum(pc$eig[pc$eig>0])
  scores<-data.frame(sample_id=rownames(pc$points),PCoA1=pc$points[,1],PCoA2=pc$points[,2],Treatment=md$Treatment,
    Inoculum=if(experiment=="microcosm")md$Inoculum else NA,replicate=md$replicate)
  loads<-data.frame(functional_group=colnames(std),PCoA1=apply(std,2,function(x)cor(x,scores$PCoA1)),PCoA2=apply(std,2,function(x)cor(x,scores$PCoA2)))
  write_table(scores,"pcoa_scores.csv");write_table(loads,"pcoa_loadings.csv")
  write_table(data.frame(axis=seq_along(var),variance_fraction=var),"ordination_variance.csv")
  mapped_columns<-intersect(mapping$ko_id,colnames(ko));detected<-mapped_columns[colSums(ko[,mapped_columns,drop=FALSE])>0]
  global<-data.frame(experiment,N=nrow(md),mapped_KOs=332,available_mapped_KO_columns=length(mapped_columns),detected_mapped_KOs=length(detected),
    mapped_groups=13,analyzed_groups=ncol(std),F,R2,p=pv,dispersion_F=dispF,dispersion_p=dispP,PCoA1=var[1],PCoA2=var[2],
    mean_targeted_coverage=mean(rowSums(ko[,mapped_columns,drop=FALSE])))
  write_table(global,"global_statistics.csv")
  if(experiment=="field") {
    write_table(md,"metadata.csv")
    renderer<-"render_field_ordination.R"
  } else {
    xy<-ff$metadata;xy$Axis1<-scores$PCoA1;xy$Axis2<-scores$PCoA2;write_table(xy,"ordination_scores.csv")
    renderer<-"render_functional_ordination.R"
  }
  stopifnot(system2("Rscript",c(shQuote(file.path(repo_root,"scripts",renderer)),shQuote(prefix)))==0)
  list(global=global,ko=ko,detected=detected)
}
field<-run_one("field",exclude_samples="P030");micro<-run_one("microcosm")
write.csv(bind_rows(field$global,micro$global),file.path(out_root,"global_statistics_all_experiments.csv"),row.names=FALSE)
updated<-mapping
for(experiment in c("field","microcosm")) {
  result<-if(experiment=="field")field else micro;ko<-result$ko
  updated[[paste0("detected_",experiment)]]<-mapping$ko_id%in%result$detected
  updated[[paste0("prevalence_",experiment)]]<-vapply(mapping$ko_id,function(k)if(k%in%colnames(ko))mean(ko[,k]>0) else 0,numeric(1))
  updated[[paste0("mean_relative_abundance_",experiment)]]<-vapply(mapping$ko_id,function(k)if(k%in%colnames(ko))mean(ko[,k]) else 0,numeric(1))
  definitions[[paste0("n_kos_detected_",experiment)]]<-vapply(definitions$functional_group,function(g)sum(mapping$functional_group==g & mapping$ko_id%in%result$detected),integer(1))
}
write.csv(updated,file.path(out_root,"frozen_332_KO_functional_group_mapping.csv"),row.names=FALSE)
write.csv(definitions,file.path(out_root,"functional_group_definitions.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out_root,"R_session_info.txt"))
