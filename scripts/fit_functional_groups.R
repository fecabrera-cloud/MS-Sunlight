# Primary microcosm functional-group mixed models and substrate contrasts.
suppressPackageStartupMessages({library(lme4);library(lmerTest);library(emmeans);library(dplyr)})
fit_functional_groups <- function(raw,md,out) {
 options(contrasts=c("contr.sum","contr.poly"))
 dir.create(out,recursive=TRUE,showWarnings=FALSE)
 write_table<-function(x,n)write.csv(x,file.path(out,n),row.names=FALSE,na="")
 origins<-c("Full sunlight","Reduced UV","Reduced UV and blue-green")
group_order <- c("Cellulose degradation","Hemicellulose degradation","Pectin degradation",
  "Alpha-glucan and fructan processing","Oligosaccharide and disaccharide utilization","Cutin degradation",
  "Photoreactivation and nucleotide-excision repair","General DNA-damage repair and tolerance",
  "Enzymatic ROS defense","Carotenoid biosynthesis and modification","Sensory photoreception","Aromatic-compound catabolism")
conditions <- list()
fit_model <- function(d,label,retry=FALSE) {
  warns <- msgs <- character()
  fit <- withCallingHandlers(lmer(y~inoculum*substrate+(1|block),data=d,REML=TRUE,
    control=if(retry)lmerControl(optimizer="nloptwrap",optCtrl=list(maxeval=200000,xtol_abs=1e-10,ftol_abs=1e-10)) else
      lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=200000))),
    warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart("muffleWarning")},
    message=function(m){msgs<<-c(msgs,conditionMessage(m));invokeRestart("muffleMessage")})
  if(length(c(warns,msgs)))conditions[[length(conditions)+1]] <<- data.frame(model=label,
    warnings=paste(unique(warns),collapse="; "),messages=paste(unique(msgs),collapse="; "))
  if(any(fit@optinfo$conv$opt!=0)) {
    if(retry)stop("Optimizer retry failed: ",label)
    return(fit_model(d,paste(label,"optimizer retry"),retry=TRUE))
  }
  stopifnot(ncol(getME(fit,"X"))==9,!any(grepl("block",colnames(getME(fit,"X")))))
  fit
}
anova_table <- function(fit,g,scale) {
  a <- as.data.frame(anova(fit,type=3,ddf="Satterthwaite"))
  stopifnot(setequal(rownames(a),c("inoculum","substrate","inoculum:substrate")),all(is.finite(a$`Pr(>F)`)))
  data.frame(functional_group=g,term=rownames(a),N=nobs(fit),numerator_df=a$NumDF,
    denominator_df=a$DenDF,F=a$`F value`,p=a$`Pr(>F)`,scale=scale,
    formula="relative_abundance ~ inoculum * substrate + (1 | block)",
    test="Type III, sum contrasts, Satterthwaite",stringsAsFactors=FALSE)
}
diagnostic_table <- function(fit,d,g,scale,pseudo) {
  r <- residuals(fit);fv <- fitted(fit);aux <- lm(I(r^2)~fv)
  vc <- as.data.frame(VarCorr(fit))
  data.frame(functional_group=g,scale,N=nrow(d),blocks=nlevels(d$block),zero_fraction=mean(d$value==0),
    singular_block=isSingular(fit,tol=1e-4),block_variance_internal=vc$vcov[vc$grp=="block"],
    residual_variance_internal=sigma(fit)^2,shapiro_p=shapiro.test(r)$p.value,
    bp_fitted_p=pchisq(nrow(d)*summary(aux)$r.squared,1,lower.tail=FALSE),
    fligner_cell_p=fligner.test(r,interaction(d$inoculum,d$substrate))$p.value,
    max_abs_scaled_residual=max(abs(r/sigma(fit))),fixed_rank=ncol(getME(fit,"X")),fixed_block_columns=0,
    optimizer_code=paste(fit@optinfo$conv$opt,collapse=";"),
    convergence_messages=paste(fit@optinfo$conv$lme4$messages,collapse="; "))
}
quiet_emm <- function(...)suppressMessages(emmeans(...,lmer.df="satterthwaite"))
means_table <- function(emm,g,kind,center,spread) {
  e <- as.data.frame(confint(emm));e$functional_group <- g;e$kind <- kind
  e$estimate_percent <- e$emmean*spread+center;e$SE_percent <- e$SE*spread
  e$lower_percent <- e$lower.CL*spread+center;e$upper_percent <- e$upper.CL*spread+center
  e$interval <- "Pointwise 95% confidence interval for expected abundance"
  e
}
pair_order <- combn(origins,2,simplify=FALSE)
pair_keys <- vapply(pair_order,paste,collapse=" - ",FUN.VALUE=character(1))
contrasts_table <- function(emm,g,kind,spread) {
  cp <- pairs(emm,adjust="tukey")
  p <- as.data.frame(summary(cp,infer=c(TRUE,TRUE)))
  rawp <- as.data.frame(summary(cp,adjust="none"))
  ii <- match(gsub("[()]","",p$contrast),pair_keys);stopifnot(!anyNA(ii))
  p$group1 <- vapply(pair_order[ii],`[`,character(1),1)
  p$group2 <- vapply(pair_order[ii],`[`,character(1),2)
  p$functional_group <- g;p$kind <- kind;p$p_unadjusted <- rawp$p.value
  names(p)[names(p)=="p.value"] <- "p_Tukey"
  p$difference_percent <- p$estimate*spread;p$SE_percent <- p$SE*spread
  p$lower_percent <- p$lower.CL*spread;p$upper_percent <- p$upper.CL*spread
  p$family <- if("inoculum"%in%names(p))paste(g,"substrate within",p$inoculum) else paste(g,"substrate averaged equally over inoculum")
  p$adjustment <- "Tukey within each family of three substrate comparisons"
  p
}

  models<-list();tests<-diagnostics<-residual_rows<-input_rows<-mean_rows<-pair_rows<-list()
  for(g in group_order) {
    d<-md;d$value<-100*raw[,g];center<-mean(d$value);spread<-sd(d$value);d$y<-(d$value-center)/spread
    fit<-fit_model(d,g);models[[g]]<-fit
    tests[[g]]<-anova_table(fit,g,"original percent")
    diagnostics[[g]]<-diagnostic_table(fit,d,g,"original percent",0)
    residual_rows[[g]]<-data.frame(functional_group=g,sample_id=d$sample_id,fitted=fitted(fit),residual=residuals(fit),scaled_residual=residuals(fit)/sigma(fit))
    input_rows[[g]]<-data.frame(functional_group=g,d)
    cell<-quiet_emm(fit,~inoculum*substrate);margin<-quiet_emm(fit,~substrate,weights="equal")
    conditional<-quiet_emm(fit,~substrate|inoculum)
    mean_rows[[length(mean_rows)+1]]<-means_table(cell,g,"inoculum-substrate cells",center,spread)
    mean_rows[[length(mean_rows)+1]]<-means_table(margin,g,"substrate marginal; equal inoculum weights",center,spread)
    pair_rows[[length(pair_rows)+1]]<-contrasts_table(margin,g,"substrate marginal",spread)
    pair_rows[[length(pair_rows)+1]]<-contrasts_table(conditional,g,"substrate within inoculum",spread)
  }
  tests<-bind_rows(tests)%>%group_by(term)%>%mutate(q_BH=p.adjust(p,"BH"))%>%ungroup();tests$omitted_sample<-NULL
  write_table(tests,"factorial_anova.csv");write_table(bind_rows(diagnostics),"model_diagnostics.csv")
  write_table(bind_rows(residual_rows),"residuals.csv");write_table(bind_rows(input_rows),"model_inputs_long.csv")
  write_table(bind_rows(mean_rows),"estimated_marginal_means.csv");write_table(bind_rows(pair_rows),"substrate_pairwise_Tukey.csv")
  saveRDS(models,file.path(out,"functional_group_models.rds"))
}
