suppressPackageStartupMessages({
  library(vegan); library(permute); library(nlme); library(emmeans)
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

fit_group_model <- function(dat,experiment) {
  if(experiment=="field") {
    fit<-lm(y~replicate+Treatment,data=dat)
    dr<-drop1(fit,test="F")
    return(list(fit=fit,F=dr["Treatment","F value"],p=dr["Treatment","Pr(>F)"],
                model_class="linear_model",random_variance=NA_real_,residual_variance=var(residuals(fit)),
                random_variance_ratio=NA_real_,near_zero_random_variance=NA))
  }
  fit<-try(lme(y~factor(replicate)+Treatment,random=~1|Inoculum,data=dat,
               method="REML",control=lmeControl(returnObject=TRUE)),silent=TRUE)
  if(inherits(fit,"try-error")) {
    fit<-lm(y~factor(Inoculum)+factor(replicate)+Treatment,data=dat)
    dr<-drop1(fit,test="F")
    return(list(fit=fit,F=dr["Treatment","F value"],p=dr["Treatment","Pr(>F)"],
                model_class="fixed_inoculum_fallback",random_variance=NA_real_,
                residual_variance=var(residuals(fit)),random_variance_ratio=NA_real_,
                near_zero_random_variance=NA))
  }
  a<-anova(fit);vc<-VarCorr(fit)
  rv<-as.numeric(vc[1,"Variance"]);ev<-as.numeric(vc[nrow(vc),"Variance"])
  ratio<-rv/ev
  list(fit=fit,F=a["Treatment","F-value"],p=a["Treatment","p-value"],
       model_class="random_inoculum_intercept",random_variance=rv,residual_variance=ev,
       random_variance_ratio=ratio,near_zero_random_variance=is.finite(ratio)&&ratio<1e-6)
}

bp_fitted_test <- function(resid,fitted) {
  if(length(unique(fitted))<2 || var(resid^2)==0) return(c(statistic=NA_real_,p=NA_real_))
  aux<-lm(I(resid^2)~fitted);stat<-length(resid)*summary(aux)$r.squared
  c(statistic=stat,p=pchisq(stat,df=1,lower.tail=FALSE))
}

run_group_models <- function(raw,md,experiment) {
  omni<-list();pairs<-list();diag<-list();letters_out<-list();loo_all<-list();z<-1
  for(g in colnames(raw)) {
    dat<-md;dat$y<-raw[,g]
    model<-fit_group_model(dat,experiment);fit<-model$fit;Fv<-model$F;pv<-model$p
    em<-emmeans(fit,~Treatment)
    pp<-as.data.frame(pairs(em,adjust="tukey"));names(pp)[names(pp)=="p.value"]<-"p_Tukey"
    pp$functional_group<-g
    tr<-levels(dat$Treatment);code<-setNames(paste0("T",seq_along(tr)),tr)
    clean_contrast<-gsub("[()]","",as.character(pp$contrast))
    cn<-vapply(strsplit(clean_contrast," - ",fixed=TRUE),function(x)paste(code[x],collapse="-"),character(1))
    ll<-multcompView::multcompLetters(setNames(pp$p_Tukey,cn),threshold=.05)$Letters
    letters_out[[z]]<-data.frame(functional_group=g,Treatment=tr,
                                 letter=unname(ll[code[tr]]),stringsAsFactors=FALSE)
    pp$contrast<-gsub(" - ","_vs_",pp$contrast,fixed=TRUE)
    pairs[[z]]<-pp
    omni[[z]]<-data.frame(functional_group=g,F=Fv,p=pv)
    r<-if(inherits(fit,"lme")) residuals(fit,type="normalized") else rstandard(fit)
    fv<-fitted(fit);bp<-bp_fitted_test(r,fv)
    shap<-if(length(r)>=3 && length(r)<=5000) shapiro.test(r)$p.value else NA_real_
    flig<-tryCatch(fligner.test(r,dat$Treatment)$p.value,error=function(e) NA_real_)
    cooks<-if(inherits(fit,"lm")) cooks.distance(fit) else rep(NA_real_,nrow(dat))
    cook_cut<-4/nrow(dat)
    loo<-lapply(seq_len(nrow(dat)),function(i){
      refit<-try(fit_group_model(dat[-i,,drop=FALSE],experiment),silent=TRUE)
      if(inherits(refit,"try-error")) return(data.frame(functional_group=g,removed_sample=rownames(dat)[i],F=NA,p=NA))
      data.frame(functional_group=g,removed_sample=rownames(dat)[i],F=refit$F,p=refit$p)
    })%>%bind_rows()
    loo_all[[z]]<-loo
    delta<-abs(loo$F-Fv);ii<-if(all(is.na(delta))) NA_integer_ else which.max(delta)
    most_influential<-if(inherits(fit,"lm") && any(is.finite(cooks))) names(which.max(cooks)) else
      if(is.na(ii)) NA_character_ else loo$removed_sample[ii]
    diag[[z]]<-data.frame(functional_group=g,N=nrow(dat),model_class=model$model_class,
      zero_fraction=mean(dat$y==0),min=min(dat$y),max=max(dat$y),
      shapiro_residual_p=shap,normality_flag=is.finite(shap)&&shap<.05,
      bp_fitted_statistic=unname(bp["statistic"]),bp_fitted_p=unname(bp["p"]),
      heteroscedasticity_flag=is.finite(bp["p"])&&bp["p"]<.05,
      fligner_treatment_p=flig,max_abs_standardized_residual=max(abs(r),na.rm=TRUE),
      cooks_cutoff=if(inherits(fit,"lm")) cook_cut else NA_real_,
      max_cooks_distance=if(inherits(fit,"lm")) max(cooks,na.rm=TRUE) else NA_real_,
      n_cooks_over_cutoff=if(inherits(fit,"lm")) sum(cooks>cook_cut,na.rm=TRUE) else NA_integer_,
      most_influential_sample=most_influential,
      loo_F_min=min(loo$F,na.rm=TRUE),loo_F_max=max(loo$F,na.rm=TRUE),
      loo_p_min=min(loo$p,na.rm=TRUE),loo_p_max=max(loo$p,na.rm=TRUE),
      random_intercept_variance=model$random_variance,residual_variance=model$residual_variance,
      random_to_residual_variance_ratio=model$random_variance_ratio,
      near_zero_random_variance=model$near_zero_random_variance)
    z<-z+1
  }
  o<-bind_rows(omni);o$p_BH<-p.adjust(o$p,"BH")
  p<-bind_rows(pairs)
  loo<-bind_rows(loo_all)%>%group_by(removed_sample)%>%mutate(p_BH=p.adjust(p,"BH"))%>%ungroup()
  d<-bind_rows(diag)
  stability<-loo%>%left_join(o%>%transmute(functional_group,primary_significant=p_BH<.05),by="functional_group")%>%
    group_by(functional_group)%>%summarise(
    loo_q_min=min(p_BH,na.rm=TRUE),loo_q_max=max(p_BH,na.rm=TRUE),
    loo_BH_significance_stable=all((p_BH<.05)==primary_significant),.groups="drop")
  d<-d%>%left_join(stability,by="functional_group")
  list(omnibus=o,pairwise=p,diagnostics=d,letters=bind_rows(letters_out),leave_one_out=loo)
}

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
  arrow_scale<-if(experiment=="microcosm") 1.20 else .84
  mul<-min(diff(range(scores$PCoA1))/max(diff(range(loadings$PCoA1)),1e-12),
           diff(range(scores$PCoA2))/max(diff(range(loadings$PCoA2)),1e-12))*arrow_scale
  loadings$xend<-loadings$PCoA1*mul;loadings$yend<-loadings$PCoA2*mul
  ann_prefix<-"PERMANOVA"
  ann<-sprintf("%s: F = %.2f, R² = %.3f, p = %s",ann_prefix,global$F,global$R2,
               format.pval(global$p,digits=2,eps=.0001))
  cent<-scores%>%group_by(Treatment)%>%summarise(PCoA1=mean(PCoA1),PCoA2=mean(PCoA2),.groups="drop")
  xr<-range(c(scores$PCoA1,loadings$xend));yr<-range(c(scores$PCoA2,loadings$yend))
  xpad<-max(diff(xr)*.18,.01);ypad<-max(diff(yr),.01)
  xlim<-c(xr[1]-xpad,xr[2]+xpad);ylim<-c(yr[1]-.10*ypad,yr[2]+.25*ypad)
  ann_x<-Inf;ann_y<-Inf;ann_hjust<-1.06;ann_vjust<-1.18
  point_size<-4.25;point_stroke<-.9
  centroid_size<-5.6;centroid_stroke<-1
  label_force_pull<-if(experiment=="microcosm") .55 else .32
  loading_labels<-
    geom_text_repel(data=loadings,aes(x=xend,y=yend,label=display_label),inherit.aes=FALSE,
                    size=3.15,fontface="bold",colour="#3F3F3F",alpha=.80,segment.color=NA,
                    box.padding=.42,point.padding=.35,max.overlaps=Inf,
                    force=2.8,force_pull=label_force_pull,max.time=8,seed=42)
  main<-ggplot(scores,aes(PCoA1,PCoA2,colour=Treatment,fill=Treatment))+
    stat_ellipse(geom="polygon",type="t",level=.75,alpha=.11,colour=NA)+
    stat_ellipse(type="t",level=.75,linetype="dotted",linewidth=.55)+
    geom_point(shape=21,size=point_size,stroke=point_stroke,colour="black")+
    geom_point(data=cent,shape=23,size=centroid_size,stroke=centroid_stroke,colour="black")+
    geom_segment(data=loadings,aes(x=0,y=0,xend=xend,yend=yend),inherit.aes=FALSE,
                 arrow=arrow(length=unit(.17,"cm")),colour="#3F3F3F",
                 alpha=.80,linewidth=.42)+
    loading_labels+
    scale_colour_manual(values=colours,breaks=levels(scores$Treatment),
                        labels=legend_labels,guide="none")+
    scale_fill_manual(values=colours,breaks=levels(scores$Treatment),labels=legend_labels)+
    annotate("label",x=ann_x,y=ann_y,label=ann,
             hjust=ann_hjust,vjust=ann_vjust,size=3.15,linewidth=.18)+
    coord_cartesian(xlim=xlim,ylim=ylim,clip="off")+
    theme_bw(base_size=13)+theme(
      panel.grid.major=element_line(colour="#EEEEEE",linewidth=.25),
      panel.grid.minor=element_blank(),legend.position="top",
      legend.box.margin=margin(0,0,0,0),plot.margin=margin(6,8,4,6),
      plot.title=element_text(face="bold",size=15,hjust=.5))+
    labs(title=ptitle,x=xlab,y=ylab,colour=tvar,fill=tvar)
  main
}

run_one <- function(experiment,output_name=experiment,exclude_samples=character()) {
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
  # Apply the prespecified, treatment-blind KO universe consistently to both
  # experiments. Relative-abundance values retain their original whole-KO-pool
  # TSS denominator; excluded columns cannot enter any downstream KO table/model.
  ko<-ko[,intersect(frozen_universe$ko_id,colnames(ko)),drop=FALSE]
  all_detected_kos<-colnames(ko)[colSums(ko,na.rm=TRUE)>0]
  mapped_columns<-intersect(mapping$ko_id,colnames(ko))
  detected_kos<-mapped_columns[colSums(ko[,mapped_columns,drop=FALSE],na.rm=TRUE)>0]
  detection_by_group<-mapping%>%
    transmute(ko_id,functional_group,detected=ko_id%in%detected_kos)%>%
    group_by(functional_group)%>%summarise(n_detected=sum(detected),.groups="drop")
  ko_detection<-data.frame(
    ko_id=mapping$ko_id,
    detected=mapping$ko_id%in%detected_kos,
    prevalence=vapply(mapping$ko_id,function(k) if(k%in%colnames(ko)) mean(ko[,k]>0) else 0,numeric(1)),
    mean_relative_abundance=vapply(mapping$ko_id,function(k) if(k%in%colnames(ko)) mean(ko[,k]) else 0,numeric(1)))
  raw_all<-aggregate_groups(ko)
  eligible<-apply(raw_all,2,sd)>0;raw<-raw_all[,eligible,drop=FALSE]
  std<-scale(raw);rownames(std)<-rownames(raw)
  d<-dist(std,"euclidean")
  set.seed(if(experiment=="field")4201 else 4202)
  fit<-adonis2(d~Treatment,data=md,permutations=9999,strata=md$PermutationBlock)
  ctrl<-how(nperm=9999);setBlocks(ctrl)<-md$PermutationBlock
  set.seed(if(experiment=="field")4301 else 4302)
  bd<-betadisper(d,md$Treatment);bp<-permutest(bd,permutations=ctrl)
  pair<-pairwise_permanova(std,md$Treatment,md$PermutationBlock)
  pc<-cmdscale(d,k=2,eig=TRUE,add=TRUE)
  var<-pc$eig[pc$eig>0]/sum(pc$eig[pc$eig>0])
  scores<-data.frame(sample_id=rownames(pc$points),PCoA1=pc$points[,1],PCoA2=pc$points[,2],
                     Treatment=md$Treatment,Inoculum=if(experiment=="microcosm")md$Inoculum else NA,
                     replicate=md$replicate)
  loads<-data.frame(functional_group=colnames(std),
                    PCoA1=apply(std,2,function(x)cor(x,scores$PCoA1)),
                    PCoA2=apply(std,2,function(x)cor(x,scores$PCoA2)))
  global<-data.frame(experiment=experiment,N=nrow(md),mapped_KOs=nrow(mapping),
                     available_mapped_KO_columns=length(mapped_columns),detected_mapped_KOs=length(detected_kos),
                     mapped_groups=ncol(raw_all),analyzed_groups=ncol(raw),F=fit$F[1],R2=fit$R2[1],p=fit$`Pr(>F)`[1],
                     dispersion_F=bp$tab[1,"F"],dispersion_p=bp$tab[1,"Pr(>F)"],
                     PCoA1=var[1],PCoA2=var[2],mean_targeted_coverage=mean(rowSums(ko[,intersect(mapping$ko_id,colnames(ko)),drop=FALSE])))
  models<-run_group_models(raw,md,experiment)
  log_sensitivity<-lapply(colnames(raw),function(g){
    y<-raw[,g];pc<-if(any(y>0)) min(y[y>0])/2 else 1e-12
    dat<-md;dat$y<-log(y+pc);m<-fit_group_model(dat,experiment)
    data.frame(functional_group=g,pseudocount=pc,F=m$F,p=m$p)
  })%>%bind_rows()%>%mutate(p_BH=p.adjust(p,"BH"))
  prefix<-file.path(out_root,output_name)
  dir.create(prefix,recursive=TRUE,showWarnings=FALSE)
  write.csv(data.frame(sample_id=rownames(raw),raw,check.names=FALSE),file.path(prefix,"functional_group_relative_abundance.csv"),row.names=FALSE)
  write.csv(data.frame(sample_id=rownames(std),std,check.names=FALSE),file.path(prefix,"functional_group_zscores.csv"),row.names=FALSE)
  write.csv(md,file.path(prefix,"metadata.csv"),row.names=FALSE)
  write.csv(global,file.path(prefix,"global_statistics.csv"),row.names=FALSE)
  write.csv(pair,file.path(prefix,"pairwise_permanova.csv"),row.names=FALSE)
  write.csv(scores,file.path(prefix,"pcoa_scores.csv"),row.names=FALSE)
  write.csv(loads,file.path(prefix,"pcoa_loadings.csv"),row.names=FALSE)
  write.csv(models$omnibus,file.path(prefix,"functional_group_omnibus.csv"),row.names=FALSE)
  write.csv(models$pairwise,file.path(prefix,"functional_group_pairwise.csv"),row.names=FALSE)
  write.csv(models$diagnostics,file.path(prefix,"model_diagnostics.csv"),row.names=FALSE)
  write.csv(models$leave_one_out,file.path(prefix,"model_influence_leave_one_out.csv"),row.names=FALSE)
  write.csv(log_sensitivity,file.path(prefix,"log_abundance_sensitivity.csv"),row.names=FALSE)
  write.csv(models$letters,file.path(prefix,"functional_group_tukey_letters.csv"),row.names=FALSE)
  excluded_here<-intersect(ids_all,exclude_samples)
  write.csv(data.frame(sample_id=excluded_here,
                       exclusion_reason=rep("Confirmed technical extraction failure; excluded before all downstream analyses",length(excluded_here))),
            file.path(prefix,"excluded_samples.csv"),row.names=FALSE)
  fig<-make_pcoa_figure(scores,loads,var,global,experiment)
  ggsave(file.path(prefix,"pcoa_publication.pdf"),fig,width=10.5,height=8)
  ggsave(file.path(prefix,"pcoa_publication.png"),fig,width=10.5,height=8,dpi=300)
  ggsave(file.path(prefix,"pcoa_publication.jpg"),fig,width=10.5,height=8,dpi=300,bg="white")
  writeLines(c(
    if(experiment=="field")
      "PERMANOVA permutations were restricted within field replicate/block."
    else
      "PERMANOVA permutations were restricted within inoculum-by-replicate blocks.",
    "A restricted permutation block controls the exchangeability of samples; it is not a random-effect term.",
    if(experiment=="microcosm")
      "Inoculum was included as a random intercept only in the univariate mixed-effects functional-group models."
    else
      "Replicate was included as a fixed blocking term in the univariate field functional-group models."
  ),file.path(prefix,"STATISTICAL_DESIGN_NOTE.txt"))
  long<-data.frame(sample_id=rownames(raw),raw,check.names=FALSE)%>%pivot_longer(-sample_id,names_to="functional_group",values_to="relative_abundance")%>%
    left_join(data.frame(sample_id=rownames(md),Treatment=md$Treatment),by="sample_id")%>%
    left_join(models$omnibus[,c("functional_group","p_BH")],by="functional_group")%>%
    mutate(label=paste0(stringr::str_wrap(functional_group,width=30),
                        "\nq = ",format.pval(p_BH,digits=2,eps=.001)))
  # Gate pairwise letters on the across-group BH-adjusted omnibus result. This
  # prevents nominal within-group Tukey contrasts from being mistaken for
  # significant endpoints when the corresponding omnibus q-value is >= 0.05.
  letter_pos<-long%>%
    group_by(functional_group,Treatment)%>%summarise(y=max(relative_abundance),.groups="drop")%>%
    group_by(functional_group)%>%mutate(pad=max(diff(range(y)),max(y)*.08,1e-8)*.22,y=y+pad)%>%ungroup()%>%
    left_join(models$letters,by=c("functional_group","Treatment"))%>%
    left_join(distinct(long,functional_group,label,p_BH),by="functional_group")%>%
    filter(p_BH<.05)
  supplementary_labels<-if(experiment=="field")
    setNames(levels(md$Treatment),levels(md$Treatment)) else
    treat_labels[levels(md$Treatment)]
  sf<-ggplot(long,aes(Treatment,relative_abundance,fill=Treatment))+geom_boxplot(width=.58,outlier.shape=NA,alpha=.65)+
    geom_jitter(width=.10,shape=21,size=1.75,colour="black")+
    geom_text(data=letter_pos,aes(Treatment,y,label=letter),inherit.aes=FALSE,size=4.4,fontface="bold")+
    facet_wrap(~label,scales="free_y",ncol=3)+
    scale_fill_manual(values=colours,breaks=levels(md$Treatment),labels=supplementary_labels)+theme_bw(base_size=9)+
    scale_y_continuous(expand=expansion(mult=c(.04,.18)))+
    theme(legend.position="top",axis.text.x=element_text(size=9,angle=35,hjust=1),
          strip.text=element_text(size=10,face="bold"),axis.title=element_text(size=11),
          panel.grid=element_blank(),plot.title=element_text(size=15,face="bold",hjust=.5))+
    labs(title=if(experiment=="field")"Functional-group responses to field sunlight exposure" else
           "Functional-group responses to litter photodegradation",
         x=NULL,y="Relative abundance",fill=if(experiment=="field")"Sunlight exposure" else "Litter photodegradation")
  ggsave(file.path(prefix,"supplementary_functional_groups.pdf"),sf,width=11,height=12)
  ggsave(file.path(prefix,"supplementary_functional_groups.png"),sf,width=11,height=12,dpi=300)
  writeLines(c(
    if(experiment=="field")
      "Field functional-group models used sunlight treatment as the focal effect and replicate as a fixed blocking term."
    else
      "Microcosm functional-group models used substrate treatment and replicate as fixed effects and inoculum as a random intercept.",
    "Facet q-values are Benjamini-Hochberg-adjusted omnibus treatment p-values across all analyzable functional groups.",
    "Pairwise letters are displayed only for functional groups whose across-group BH-adjusted omnibus q-value is below 0.05.",
    "For those significant groups, letters are calculated directly from Tukey-adjusted pairwise p-values within the group; no second BH adjustment is applied to the already Tukey-adjusted contrasts."
  ),file.path(prefix,"SUPPLEMENTARY_FUNCTIONAL_GROUPS_CAPTION.txt"))

  load_long<-loads%>%pivot_longer(c(PCoA1,PCoA2),names_to="axis",values_to="correlation")
  load_long$functional_group<-factor(load_long$functional_group,
    levels=loads$functional_group[order(loads$PCoA1)])
  lp<-ggplot(load_long,aes(functional_group,correlation,fill=axis))+geom_col(position="dodge",width=.72)+
    geom_hline(yintercept=0,linewidth=.35)+coord_flip()+
    scale_fill_manual(values=c(PCoA1="#276FBF",PCoA2="#D1495B"),labels=c(PCoA1="PCoA1",PCoA2="PCoA2"))+
    theme_bw(base_size=11)+theme(panel.grid.major.y=element_blank(),legend.position="top",
      plot.title=element_text(face="bold",hjust=.5))+labs(title="Functional-group correlations with ordination axes",
      x=NULL,y="Pearson correlation with PCoA score",fill=NULL)
  ggsave(file.path(prefix,"pcoa_loadings_barplot.pdf"),lp,width=9,height=7)
  ggsave(file.path(prefix,"pcoa_loadings_barplot.png"),lp,width=9,height=7,dpi=300)

  groups<-colnames(raw);upper<-list();zz<-1
  for(i in seq_len(length(groups)-1))for(j in (i+1):length(groups)){
    ct<-suppressWarnings(cor.test(raw[,groups[i]],raw[,groups[j]],method="spearman",exact=FALSE))
    upper[[zz]]<-data.frame(group1=groups[i],group2=groups[j],rho=unname(ct$estimate),p=ct$p.value);zz<-zz+1
  }
  cu<-bind_rows(upper);cu$p_BH<-p.adjust(cu$p,"BH")
  corr<-bind_rows(cu,transmute(cu,group1=.data$group2,group2=.data$group1,rho=rho,p=p,p_BH=p_BH),
                  data.frame(group1=groups,group2=groups,rho=1,p=0,p_BH=0))
  corr$mark<-ifelse(corr$p_BH<.001,"***",ifelse(corr$p_BH<.01,"**",ifelse(corr$p_BH<.05,"*","")))
  write.csv(corr,file.path(prefix,"functional_group_spearman_correlations.csv"),row.names=FALSE)
  corr$group1<-factor(corr$group1,levels=rev(groups));corr$group2<-factor(corr$group2,levels=groups)
  cp<-ggplot(corr,aes(group2,group1,fill=rho))+geom_tile(colour="white",linewidth=.25)+
    geom_text(aes(label=ifelse(group1==group2,"1",paste0(sprintf("%.2f",rho),mark))),size=2.6)+
    scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",midpoint=0,limits=c(-1,1))+
    coord_equal()+theme_minimal(base_size=9)+theme(axis.title=element_blank(),panel.grid=element_blank(),
      axis.text.x=element_text(angle=45,hjust=1),plot.title=element_text(face="bold",hjust=.5))+
    labs(title="Spearman correlations among targeted functional groups",fill="rho")
  ggsave(file.path(prefix,"functional_group_correlation_matrix.pdf"),cp,width=11,height=10)
  ggsave(file.path(prefix,"functional_group_correlation_matrix.png"),cp,width=11,height=10,dpi=300)
  writeLines(c("Correlations are descriptive associations among targeted relative-abundance sums.",
    "Asterisks denote Benjamini-Hochberg-adjusted p-values across unique group pairs (* q<0.05, ** q<0.01, *** q<0.001).",
    "Negative correlations do not by themselves demonstrate physiological trade-offs, and compositional closure can induce correlations."),
    file.path(prefix,"CORRELATION_INTERPRETATION.txt"))

  # Experiment-specific ungrouped KO audit; 75% prevalence, CLR, design-adjusted linear model.
  prev<-colMeans(ko>0);keep<-prev>=.75;kk<-ko[,keep,drop=FALSE]
  pcnt<-min(kk[kk>0],na.rm=TRUE)/2;clr<-log(kk+pcnt);clr<-clr-rowMeans(clr)
  test_one<-function(y){
    dat<-md;dat$y<-y
    if(experiment=="field")f<-lm(y~replicate+Treatment,data=dat) else f<-lm(y~Inoculum+replicate+Treatment,data=dat)
    a<-drop1(f,test="F");c(F=a["Treatment","F value"],p=a["Treatment","Pr(>F)"])
  }
  kt<-t(apply(clr,2,test_one));kores<-data.frame(ko_id=rownames(kt),prevalence=prev[rownames(kt)],F=kt[,"F"],p=kt[,"p"])
  kores$p_BH<-p.adjust(kores$p,"BH");kores<-kores[order(kores$p_BH,kores$p),]
  write.csv(kores,file.path(prefix,"ungrouped_KO_75pct_omnibus.csv"),row.names=FALSE)
  ko_screen_summary<-data.frame(experiment=experiment,N=nrow(md),prevalence_threshold=.75,
    eligible_KOs=nrow(kores),significant_KOs_BH=sum(kores$p_BH<.05),
    minimum_raw_p=min(kores$p),minimum_BH_q=min(kores$p_BH),
    status="exploratory secondary screen")
  writeLines(c(
    "The 75%-prevalence KO-level analysis is an exploratory secondary screen, not the primary functional-group endpoint.",
    "KO relative abundances were zero-replaced with half the smallest positive value, CLR transformed, and tested with design-adjusted linear models.",
    if(experiment=="field")
      "Field models included replicate and treatment as fixed effects."
    else
      "Microcosm models included inoculum, replicate, and treatment as fixed effects; this differs from the random-inoculum primary group models and should be interpreted as a sensitivity/discovery screen.",
    "Benjamini-Hochberg correction was applied across all prevalence-eligible KOs within each experiment.",
    sprintf("Result: %d of %d eligible KOs had q < 0.05 (minimum q = %s).",
            ko_screen_summary$significant_KOs_BH,ko_screen_summary$eligible_KOs,
            format.pval(ko_screen_summary$minimum_BH_q,digits=3,eps=1e-4)),
    "Individual KO hits require separate functional interpretation and are not evidence that an entire functional group changed."
  ),file.path(prefix,"KO_SCREEN_INTERPRETATION.txt"))
  list(global=global,omnibus=models$omnibus,pairwise=pair,diagnostics=models$diagnostics,
       log_sensitivity=log_sensitivity,
       detection_by_group=detection_by_group,ko_detection=ko_detection,ko_screen=ko_screen_summary,
       all_detected_kos=all_detected_kos)
}

field<-run_one("field","field",exclude_samples="P030")
micro<-run_one("microcosm")
write.csv(bind_rows(field$global,micro$global),file.path(out_root,"global_statistics_all_experiments.csv"),row.names=FALSE)
write.csv(bind_rows(field$ko_screen,micro$ko_screen),file.path(out_root,"ko_screen_summary.csv"),row.names=FALSE)
stopifnot(length(union(field$all_detected_kos,micro$all_detected_kos))==7931)

# Recompute every experiment-specific detection field from the final analyzed sample sets.
mapping_out<-mapping
mapping_out$detected_field<-field$ko_detection$detected
mapping_out$prevalence_field<-field$ko_detection$prevalence
mapping_out$mean_relative_abundance_field<-field$ko_detection$mean_relative_abundance
mapping_out$detected_microcosm<-micro$ko_detection$detected
mapping_out$prevalence_microcosm<-micro$ko_detection$prevalence
mapping_out$mean_relative_abundance_microcosm<-micro$ko_detection$mean_relative_abundance
mapping_out$detected_in_experiment<-NULL
mapping_out$included_in_group_sum<-NULL
write.csv(mapping_out,file.path(out_root,"frozen_332_KO_functional_group_mapping.csv"),row.names=FALSE)
write.csv(mapping%>%count(functional_group,name="n_KOs")%>%arrange(functional_group),
          file.path(out_root,"functional_group_membership_counts.csv"),row.names=FALSE)
definitions$n_kos_detected_field<-field$detection_by_group$n_detected[match(definitions$functional_group,field$detection_by_group$functional_group)]
definitions$n_kos_detected_microcosm<-micro$detection_by_group$n_detected[match(definitions$functional_group,micro$detection_by_group$functional_group)]
write.csv(definitions,file.path(out_root,"functional_group_definitions.csv"),row.names=FALSE)

diagnostic_summary<-function(result,experiment) {
  x<-result$diagnostics
  primary<-result$omnibus%>%transmute(functional_group,primary_significant=p_BH<.05)
  logs<-result$log_sensitivity%>%transmute(functional_group,log_significant=p_BH<.05)
  discordant<-primary%>%left_join(logs,by="functional_group")%>%filter(primary_significant!=log_significant)%>%pull(functional_group)
  data.frame(experiment=experiment,groups=nrow(x),
    residual_normality_flags=sum(x$normality_flag,na.rm=TRUE),
    heteroscedasticity_flags=sum(x$heteroscedasticity_flag,na.rm=TRUE),
    influence_flags=sum(!x$loo_BH_significance_stable,na.rm=TRUE),
    near_zero_random_variance_flags=sum(x$near_zero_random_variance,na.rm=TRUE),
    fallback_models=sum(x$model_class=="fixed_inoculum_fallback",na.rm=TRUE),
    log_sensitivity_discordant_groups=length(discordant),
    discordant_group_names=paste(discordant,collapse="; "))
}
diag_summary<-bind_rows(diagnostic_summary(field,"field"),diagnostic_summary(micro,"microcosm"))
write.csv(diag_summary,file.path(out_root,"model_diagnostics_summary.csv"),row.names=FALSE)
writeLines(c(
  "Functional-group model diagnostics",
  "Residual normality was screened with Shapiro-Wilk tests on standardized residuals.",
  "Heteroscedasticity was screened with a Breusch-Pagan-style auxiliary regression of squared standardized residuals on fitted values; Fligner-Killeen treatment tests are also reported.",
  "Field influence was summarized by Cook's distance and all models were checked by leave-one-sample-out refitting.",
  "Microcosm influence was assessed by leave-one-sample-out refitting because standard Cook's distances are not defined for these nlme fits.",
  "Leave-one-out omnibus p-values were BH-adjusted across groups for each omitted sample; stability indicates whether the q<0.05 conclusion was unchanged in every refit.",
  "Near-zero random-intercept variance was defined as an inoculum-to-residual variance ratio below 1e-6.",
  "A log(relative abundance + half-minimum-positive) model with the same design was used as a transformation sensitivity analysis; it is not a replacement primary model.",
  apply(diag_summary,1,function(r)sprintf("%s: %s normality flag(s), %s heteroscedasticity flag(s), %s BH-conclusion influence flag(s), %s near-zero random-variance flag(s), %s fallback model(s), and %s group(s) discordant under log sensitivity [%s].",
    r[["experiment"]],r[["residual_normality_flags"]],r[["heteroscedasticity_flags"]],
    r[["influence_flags"]],r[["near_zero_random_variance_flags"]],r[["fallback_models"]],
    r[["log_sensitivity_discordant_groups"]],ifelse(nchar(r[["discordant_group_names"]])==0,"none",r[["discordant_group_names"]])))
),file.path(out_root,"MODEL_DIAGNOSTICS_INTERPRETATION.txt"))

flag_names<-function(x,condition) {
  z<-x$functional_group[condition & !is.na(condition)]
  if(length(z)==0) "none" else paste(z,collapse="; ")
}
writeLines(c(
  "# Functional metagenomics: methods and diagnostic report","",
  "## Analysis scope","",
  "The primary endpoint is the relative representation of 332 curated KOs assigned once each to 13 hypothesis-driven functional groups. The profile is targeted (approximately 5% of total KO relative abundance), not an exhaustive metagenomic functional profile.","",
  sprintf("Both experiments were restricted to the same prespecified, treatment-blind universe of %d KOs before candidate screening or statistical analysis.",nrow(frozen_universe)),"",
  "A mapped KO is defined as detected only when it has nonzero abundance in at least one retained sample of the corresponding final analysis. `available_mapped_KO_columns` separately records mapped KOs represented by a matrix column, including all-zero columns. Detection, prevalence, mean abundance, and per-group detected counts are recomputed after final sample filtering.","",
  sprintf("Final detected mapped KOs: field %d of %d available columns; microcosm %d of %d available columns.",field$global$detected_mapped_KOs,field$global$available_mapped_KO_columns,micro$global$detected_mapped_KOs,micro$global$available_mapped_KO_columns),"",
  "## Primary functional-group inference","",
  "Group relative abundances are sums of total-sum-scaled KO relative abundances. Univariate field models include replicate and treatment as fixed effects; microcosm models include replicate and treatment as fixed effects and inoculum as a random intercept. Omnibus p-values are BH-adjusted across all analyzable groups.","",
  "Pairwise contrasts use Tukey-adjusted p-values directly. Figure letters are calculated from those Tukey-adjusted values. A second BH adjustment is not applied to already Tukey-adjusted pairwise p-values.","",
  "Z-score standardization is used only to give groups equal weight in the multivariate ordination/PERMANOVA. Standardizing a single response would change its units but not its F statistic or p-value, and would not resolve compositionality or residual-model assumptions; it is therefore not used for the univariate models.","",
  "## Secondary KO-level screen","",
  "The 75%-prevalence KO analysis is explicitly exploratory. Eligible KOs are zero-replaced with half the smallest positive abundance, CLR transformed, and tested with design-adjusted linear models; BH correction is applied across eligible KOs within each experiment. It is a discovery/sensitivity screen and does not supersede the curated functional-group endpoint.","",
  sprintf("Field KO screen: %d of %d eligible KOs at q < 0.05. Microcosm KO screen: %d of %d eligible KOs at q < 0.05.",field$ko_screen$significant_KOs_BH,field$ko_screen$eligible_KOs,micro$ko_screen$significant_KOs_BH,micro$ko_screen$eligible_KOs),"",
  "## Diagnostics and robustness","",
  sprintf("Field residual-normality flags: %s.",flag_names(field$diagnostics,field$diagnostics$normality_flag)),
  sprintf("Field fitted-value heteroscedasticity flags: %s.",flag_names(field$diagnostics,field$diagnostics$heteroscedasticity_flag)),
  sprintf("Field leave-one-out BH-instability flags: %s.",flag_names(field$diagnostics,!field$diagnostics$loo_BH_significance_stable)),
  sprintf("Microcosm residual-normality flags: %s.",flag_names(micro$diagnostics,micro$diagnostics$normality_flag)),
  sprintf("Microcosm fitted-value heteroscedasticity flags: %s.",flag_names(micro$diagnostics,micro$diagnostics$heteroscedasticity_flag)),
  sprintf("Microcosm leave-one-out BH-instability flags: %s.",flag_names(micro$diagnostics,!micro$diagnostics$loo_BH_significance_stable)),
  sprintf("Microcosm near-zero inoculum random-variance flags: %s.",flag_names(micro$diagnostics,micro$diagnostics$near_zero_random_variance)),"",
  "The log-abundance sensitivity preserves every field nonsignificant conclusion. In the microcosm, five of the six primary significant groups remain significant after log transformation; pectin degradation changes from q < 0.05 to q = 0.069. Leave-one-out analysis also shows that the photoreactivation/NER conclusion is sample-sensitive. The most robust microcosm univariate results are therefore cellulose degradation, hemicellulose degradation, alpha-glucan/fructan processing, and oligosaccharide/disaccharide utilization. Pectin and photoreactivation/NER should be reported as primary-model findings with explicit sensitivity caveats.","",
  "The significant field multivariate result remains subject to heterogeneous dispersion, and no field group is BH-significant in either the raw-abundance primary models or the log-abundance sensitivity."
),file.path(out_root,"FUNCTIONAL_METHODS_AND_DIAGNOSTIC_REPORT.md"))
writeLines(capture.output(sessionInfo()),file.path(out_root,"R_session_info.txt"))
