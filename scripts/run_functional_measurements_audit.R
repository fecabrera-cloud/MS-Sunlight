suppressPackageStartupMessages({library(nlme);library(emmeans);library(ggplot2);library(dplyr);library(patchwork)})
set.seed(42)
script_arg<-commandArgs(trailingOnly=FALSE)[grep("^--file=",commandArgs(trailingOnly=FALSE))]
repo_root<-normalizePath(file.path(dirname(sub("^--file=", "", script_arg)),".."),mustWork=TRUE)
out_root<-file.path(repo_root,"results","reproduced","functional_measurements");dir.create(out_root,recursive=TRUE,showWarnings=FALSE)
dat<-read.csv(file.path(repo_root,"data","functional_measurements","cleaned_functional_measurements.csv"),stringsAsFactors=FALSE)
dat<-dat[!dat$is_control,]
dat$Treatment<-factor(dat$substrate_code,levels=c("S+","UV-","S-"));dat$inoculum<-factor(dat$inoculum);dat$block<-factor(dat$block)
outs<-data.frame(key=c("decomposition","biomass","beta_glucosidase","cue_proxy"),
                 column=c("decomposition_pct","microbial_biomass_c_ug_g","beta_glucosidase","cue_proxy"),
                 label=c("Dry-mass loss (%)","SIR-derived microbial biomass C (ug C g^-1)",
                         "Beta-glucosidase activity (umol g^-1 h^-1)",
                         "CUE proxy (biomass C / cumulative respiration)"),
                 log=c(FALSE,TRUE,TRUE,TRUE),stringsAsFactors=FALSE)
omni<-list();pairs_out<-list();emm_out<-list();plots<-list();fits<-list()
pal<-c("S+"="#E3AC18","UV-"="#77BCE8","S-"="#555555")
labs<-c("S+"="Full sunlight","UV-"="Reduced UV","S-"="Reduced UV and blue-green")
for(i in seq_len(nrow(outs))){
 o<-outs[i,];d<-dat[!is.na(dat[[o$column]]),];d$y<-d[[o$column]];if(o$log)d$y<-log(d$y)
 fit<-lme(y~block+Treatment,random=~1|inoculum,data=d,method="REML",control=lmeControl(returnObject=TRUE));fits[[o$key]]<-fit
 a<-anova(fit);pv<-a["Treatment","p-value"];fv<-a["Treatment","F-value"]
 omni[[i]]<-data.frame(outcome=o$key,N=nrow(d),F=fv,p=pv,random_inoculum_variance=as.numeric(VarCorr(fit)[1,"Variance"]))
 em<-emmeans(fit,~Treatment);ee<-as.data.frame(confint(em));pp<-as.data.frame(pairs(em,adjust="tukey"));pp$outcome<-o$key
 if(o$log){ee$estimate<-exp(ee$emmean);ee$lower<-exp(ee$lower.CL);ee$upper<-exp(ee$upper.CL)}else{ee$estimate<-ee$emmean;ee$lower<-ee$lower.CL;ee$upper<-ee$upper.CL}
 # Compact-letter display from the Tukey-adjusted pairwise p-values. Treatment
 # codes are temporarily sanitized because the original codes contain hyphens.
 tr<-as.character(ee$Treatment);code<-setNames(paste0("T",seq_along(tr)),tr)
 clean_contrast<-gsub("[()]","",as.character(pp$contrast))
 cn<-vapply(strsplit(clean_contrast," - ",fixed=TRUE),function(x)paste(code[x],collapse="-"),character(1))
 lp<-setNames(pp$p.value,cn);ll<-multcompView::multcompLetters(lp,threshold=.05)$Letters
 ee$letter<-unname(ll[code[tr]]);ee$letter_y<-ee$upper*1.07
 ee$outcome<-o$key;emm_out[[i]]<-ee;pairs_out[[i]]<-pp
 plots[[i]]<-ggplot(d,aes(Treatment,.data[[o$column]],fill=Treatment))+
   geom_boxplot(width=.52,outlier.shape=NA,alpha=.25)+geom_jitter(width=.09,shape=21,size=2.1,colour="black")+
   geom_point(data=ee,aes(Treatment,estimate),inherit.aes=FALSE,shape=23,size=3.8,fill="white")+
   geom_errorbar(data=ee,aes(Treatment,ymin=lower,ymax=upper),inherit.aes=FALSE,width=.08)+
   geom_text(data=ee,aes(Treatment,letter_y,label=letter),inherit.aes=FALSE,fontface="bold",size=4)+
   scale_y_continuous(expand=expansion(mult=c(.04,.18)))+
   scale_fill_manual(values=pal,labels=labs)+theme_classic(base_size=11)+theme(legend.position="none",axis.title.x=element_blank())+
   labs(y=o$label,title=LETTERS[i])
}
om<-bind_rows(omni);om$p_BH<-p.adjust(om$p,"BH");write.csv(om,file.path(out_root,"omnibus_tests.csv"),row.names=FALSE)
write.csv(bind_rows(pairs_out),file.path(out_root,"pairwise_contrasts.csv"),row.names=FALSE)
write.csv(bind_rows(emm_out),file.path(out_root,"estimated_marginal_means.csv"),row.names=FALSE)
saveRDS(fits,file.path(out_root,"models.rds"))
for(i in seq_along(plots)){
  plots[[i]]<-plots[[i]]+annotate("text",x=2,y=Inf,
    label=paste0("Treatment q = ",format.pval(om$p_BH[i],digits=2,eps=.0001)),
    vjust=1.15,size=3.2)
}
fig<-wrap_plots(plots,ncol=2)+plot_annotation(title="Microbial responses to litter photodegradation")
ggsave(file.path(out_root,"primary_functional_responses.pdf"),fig,width=10,height=8)
ggsave(file.path(out_root,"primary_functional_responses.png"),fig,width=10,height=8,dpi=300)
write.csv(data.frame(substrate=names(table(dat$Treatment)),N=as.integer(table(dat$Treatment))),file.path(out_root,"sample_counts.csv"),row.names=FALSE)
writeLines(c(
  "Figure. Microbial responses during dark incubation of photochemically conditioned maize litter.",
  "Points are individual inoculated microcosms; boxplots summarize observed values. White diamonds and error bars are model-estimated marginal means and 95% confidence intervals from block-adjusted linear mixed-effects models with inoculum as a random intercept. Letters denote Tukey-adjusted pairwise differences within each response, and q values are Benjamini-Hochberg-adjusted omnibus treatment p-values across the four responses. The CUE proxy is SIR-derived microbial biomass C divided by cumulative respiration and should not be interpreted as a direct physiological measurement of carbon-use efficiency. N = 46 for dry-mass loss, microbial biomass C, and beta-glucosidase activity; N = 43 for the CUE proxy."
),file.path(out_root,"FIGURE_CAPTION.txt"))
writeLines(capture.output(sessionInfo()),file.path(out_root,"R_session_info.txt"))
