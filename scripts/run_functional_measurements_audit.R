# Microcosm functional measurements: full factorial models with random block.
suppressPackageStartupMessages({library(lme4);library(lmerTest);library(emmeans);library(ggplot2);library(patchwork);library(dplyr);library(jsonlite)})
options(contrasts=c("contr.sum","contr.poly"))
script<-normalizePath(sub("^--file=","",commandArgs(FALSE)[grep("^--file=",commandArgs(FALSE))]))
repo_root<-dirname(dirname(script));source(file.path(repo_root,"scripts","factorial_permanova.R"))
out<-file.path(repo_root,"results","reproduced","functional_measurements");dir.create(out,recursive=TRUE,showWarnings=FALSE)
data_out<-out
write_table<-function(x,n)write.csv(x,file.path(out,n),row.names=FALSE,na="")
input_data<-read.csv(file.path(repo_root,"data","functional_measurements","cleaned_functional_measurements.csv"))
columns<-c(decomposition="decomposition_pct",biomass="microbial_biomass_c_ug_g",beta_glucosidase="beta_glucosidase",cue_proxy="cue_proxy")
origins <- c("Full sunlight","Reduced UV","Reduced UV and blue-green")
palette <- setNames(c("#E3AC18","#77BCE8","#555555"),origins)
shapes <- setNames(c(21,24,22),origins);offsets <- setNames(c(-.14,0,.14),origins)
font_sizes <- list(panel_title=17,axis_title=15,axis_text=13.5,ANOVA=13.5,
  bracket=12.5,legend_title=14,legend_text=13.5)
spec <- data.frame(outcome=c("decomposition","biomass","beta_glucosidase","cue_proxy"),
  title=c("A  Dry-mass loss","B  Microbial biomass C","C  β-Glucosidase activity","D  Biomass-to-respiration\n     ratio"),
  ylabel=c("Dry-mass loss (%)","SIR-derived biomass C (µg C g⁻¹)",
    "Potential activity (µmol g⁻¹ h⁻¹)","Biomass C / cumulative respiration"),
  log=c(FALSE,TRUE,TRUE,TRUE),N=c(46,46,46,43))
pair_levels <- combn(origins,2,simplify=FALSE)
pair_keys <- vapply(pair_levels,paste,collapse=" - ",FUN.VALUE=character(1))
models <- list();raw_rows <- stat_rows <- mean_rows <- pair_rows <- conditional_rows <- check_rows <- list()
for(j in seq_len(nrow(spec))) {
 s <- spec[j,];d0<-input_data[!input_data$is_control & is.finite(input_data[[columns[[s$outcome]]]]),]
 d<-data.frame(sample_id=as.character(d0$id),inoculum=full_treatment(d0$inoculum),substrate=full_treatment(d0$substrate_code),block=factor(d0$block),value=d0[[columns[[s$outcome]]]])
 d$y<-if(s$log)log(d$value) else d$value
 for(v in c("inoculum","substrate"))d[[v]] <- factor(d[[v]],levels=origins)
 d$block <- factor(d$block)
 stopifnot(nrow(d)==s$N,nlevels(d$block)==6,!anyNA(d),
   max(abs(d$y-if(s$log)log(d$value) else d$value))<1e-10)
 fit <- lmer(y~inoculum*substrate+(1|block),data=d,REML=TRUE)
 stopifnot(ncol(getME(fit,"X"))==9,!any(grepl("block",colnames(getME(fit,"X")))),all(fit@optinfo$conv$opt==0))
 models[[s$outcome]] <- fit
 a <- as.data.frame(anova(fit,type=3,ddf="Satterthwaite"))
 stats <- data.frame(outcome=s$outcome,term=rownames(a),N=nrow(d),numerator_df=a$NumDF,
   denominator_df=a$DenDF,F=a$`F value`,p=a$`Pr(>F)`,transformation=if(s$log)"natural log" else "original percentage",
   formula="response ~ inoculum * substrate + (1 | block)")
 stat_rows[[j]] <- stats
 em <- suppressMessages(emmeans(fit,~substrate,weights="equal",lmer.df="satterthwaite"))
 cell <- as.data.frame(emmeans(fit,~inoculum*substrate,lmer.df="satterthwaite"))
 ee <- as.data.frame(confint(em));ee$outcome <- s$outcome
 expected <- vapply(as.character(ee$substrate),function(v)mean(cell$emmean[cell$substrate==v]),numeric(1))
 stopifnot(max(abs(expected-ee$emmean))<1e-10)
 ee$estimate <- if(s$log)exp(ee$emmean) else ee$emmean
 ee$lower <- if(s$log)exp(ee$lower.CL) else ee$lower.CL
 ee$upper <- if(s$log)exp(ee$upper.CL) else ee$upper.CL
 ee$scale <- if(s$log)"back-transformed geometric mean; equal inoculum weights on log scale" else "original percentage"
 ee$inoculum_weights <- "1/3, 1/3, 1/3"
 mean_rows[[j]] <- ee
 get_pairs <- function(emm,kind) {
   pp <- as.data.frame(summary(pairs(emm,adjust="tukey"),infer=c(TRUE,TRUE)))
   rawp <- as.data.frame(summary(pairs(emm,adjust="none")))
   ii <- match(gsub("[()]","",pp$contrast),pair_keys);stopifnot(!anyNA(ii))
   pp$group1 <- vapply(pair_levels[ii],`[`,character(1),1)
   pp$group2 <- vapply(pair_levels[ii],`[`,character(1),2)
   pp$outcome <- s$outcome;pp$kind <- kind;pp$p_unadjusted <- rawp$p.value
   names(pp)[names(pp)=="p.value"] <- "p_Tukey"
   pp$effect_scale <- if(s$log)"log ratio" else "percentage-point difference"
   pp$effect <- if(s$log)exp(pp$estimate) else pp$estimate
   pp$effect_lower <- if(s$log)exp(pp$lower.CL) else pp$lower.CL
   pp$effect_upper <- if(s$log)exp(pp$upper.CL) else pp$upper.CL
   pp
 }
 pair_rows[[j]] <- get_pairs(em,"substrate averaged equally over inoculum")
 conditional <- get_pairs(emmeans(fit,~substrate|inoculum,lmer.df="satterthwaite"),"substrate within inoculum")
 conditional_rows[[j]] <- conditional
 d$outcome <- s$outcome;raw_rows[[j]] <- d
 check_rows[[j]] <- data.frame(outcome=s$outcome,N=nrow(d),blocks=nlevels(d$block),
   fixed_columns=ncol(getME(fit,"X")),fixed_block_columns=0)
}
statistics <- bind_rows(stat_rows)%>%group_by(term)%>%mutate(p_BH=p.adjust(p,"BH"))%>%ungroup()
raw <- bind_rows(raw_rows);means <- bind_rows(mean_rows);contrasts <- bind_rows(pair_rows)
for(v in c("substrate","inoculum"))raw[[v]] <- factor(raw[[v]],levels=origins)
means$substrate <- factor(means$substrate,levels=origins)
statistics <- statistics%>%mutate(significance=ifelse(p_BH<.05,"S","NS"))
stopifnot(all(statistics$significance[statistics$term=="substrate"]=="S"),
  all(statistics$significance[statistics$term!="substrate"]=="NS"))
set.seed(42)
raw$x <- as.integer(raw$substrate)+unname(offsets[as.character(raw$inoculum)])+runif(nrow(raw),-.035,.035)
means$x <- as.integer(means$substrate)
panels <- list();bracket_rows <- annotation_rows <- plot_checks <- list()
for(j in seq_len(nrow(spec))) {
 s <- spec[j,];d <- filter(raw,outcome==s$outcome);e <- filter(means,outcome==s$outcome)
 a <- filter(statistics,outcome==s$outcome)
 a <- a[match(c("substrate","inoculum","inoculum:substrate"),a$term),]
 a$label <- paste0(c("Substrate","Inoculum","Interaction"),": ",a$significance)
 cp <- filter(contrasts,outcome==s$outcome,p_Tukey<.05)
 stopifnot(nrow(cp)==if(s$outcome=="decomposition")0 else 2)
 yr <- range(c(d$value,e$lower,e$upper));span <- diff(yr)
 ymin <- max(0,yr[1]-.07*span)
 b <- cp
 b$x1 <- match(b$group1,origins);b$x2 <- match(b$group2,origins)
 b <- b[order(b$x2-b$x1),]
 b$y <- yr[2]+(.13+.20*seq_len(nrow(b))-.20)*span
 b$label <- as.character(ifelse(b$p_Tukey<.0001,"p < 0.0001",sprintf("p = %.4f",b$p_Tukey)))
 b$label_y <- b$y+.035*span
 b$end_y <- b$y-.025*span
 header_bottom <- if(nrow(b))yr[2]+.53*span else yr[2]+.13*span
 a$x <- .56;a$y <- header_bottom+c(.22,.11,0)*span
 a$fontface <- ifelse(a$term=="substrate","bold","plain")
 ymax <- header_bottom+.32*span
 bracket_rows[[j]] <- b;annotation_rows[[j]] <- a
 p <- ggplot()+
   geom_point(data=d,aes(x=x,y=value,fill=substrate,shape=inoculum),size=2.8,stroke=.65,alpha=.65)+
   geom_errorbar(data=e,aes(x=x,ymin=lower,ymax=upper),width=.10,linewidth=.8,colour="black",show.legend=FALSE)+
   geom_point(data=e,aes(x=x,y=estimate),shape=23,fill="white",colour="black",size=4.8,stroke=1,show.legend=FALSE)+
   geom_text(data=a,aes(x=x,y=y,label=label,fontface=fontface),inherit.aes=FALSE,hjust=0,
     size=font_sizes$ANOVA/ggplot2::.pt,colour="#333333",show.legend=FALSE)+
   scale_fill_manual(values=palette,breaks=origins,name="Microcosm substrate",
     guide=guide_legend(order=1,nrow=1,title.position="left",override.aes=list(shape=21,alpha=1,size=3.5)))+
   scale_shape_manual(values=shapes,breaks=origins,name="Inoculum origin",
     guide=guide_legend(order=2,nrow=1,title.position="left",override.aes=list(fill="white",alpha=1,size=3.5)))+
   scale_x_continuous(breaks=1:3,labels=c("Full\nsunlight","Reduced\nUV","Reduced UV and\nblue-green"),limits=c(.5,3.5))+
   scale_y_continuous(breaks=scales::breaks_pretty(n=5),expand=expansion(mult=c(0,0)))+
   coord_cartesian(ylim=c(ymin,ymax),clip="on")+
   labs(title=s$title,x=NULL,y=s$ylabel)+theme_bw(base_size=15)+
   theme(panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),
     panel.grid.major.y=element_line(colour="#EEEEEE",linewidth=.25),
     axis.text=element_text(size=font_sizes$axis_text,colour="#333333"),axis.title.y=element_text(size=font_sizes$axis_title),
     plot.title=element_text(size=font_sizes$panel_title,face="bold",margin=margin(b=7)),
     legend.position="bottom",legend.box="vertical",legend.title=element_text(size=font_sizes$legend_title),
     legend.text=element_text(size=font_sizes$legend_text),legend.margin=margin(2,0,2,0),plot.margin=margin(8,9,7,8))
 if(nrow(b))p <- p+
   geom_segment(data=b,aes(x=x1,xend=x2,y=y,yend=y),linewidth=.5,inherit.aes=FALSE)+
   geom_segment(data=b,aes(x=x1,xend=x1,y=y,yend=end_y),linewidth=.5,inherit.aes=FALSE)+
   geom_segment(data=b,aes(x=x2,xend=x2,y=y,yend=end_y),linewidth=.5,inherit.aes=FALSE)+
   geom_text(data=b,aes(x=(x1+x2)/2,y=label_y,label=label),size=font_sizes$bracket/ggplot2::.pt,vjust=0,inherit.aes=FALSE)
 built <- ggplot_build(p)
 stopifnot(nrow(built$data[[1]])==s$N,nrow(built$data[[2]])==3,nrow(built$data[[3]])==3,
   identical(sort(unique(built$data[[1]]$shape)),c(21,22,24)),
   identical(built$data[[4]]$label,a$label),all(d$x>.5 & d$x<3.5))
 if(nrow(b))stopifnot(nrow(built$data[[5]])==2,
   all(built$data[[6]]$yend<built$data[[6]]$y),all(built$data[[7]]$yend<built$data[[7]]$y))
 plot_checks[[j]] <- data.frame(outcome=s$outcome,plotted_samples=nrow(built$data[[1]]),
   plotted_marginal_means=nrow(built$data[[3]]),ANOVA_labels=nrow(built$data[[4]]),
   significant_substrate_brackets=nrow(b),visible_bracket_caps=if(nrow(b))all(b$end_y<b$y) else TRUE)
 panels[[s$outcome]] <- p
}
figure <- wrap_plots(panels,ncol=2,guides="collect")+
 plot_annotation(theme=theme(plot.margin=margin(8,8,8,8)))
figure <- figure & theme(legend.position="bottom",legend.box="vertical")
for(ext in c("pdf","png")) {
 if(ext=="pdf")ggsave(file.path(out,"Figure6.pdf"),figure,width=11.8,height=10.2,device=cairo_pdf,bg="white") else
   ggsave(file.path(out,"Figure6.png"),figure,width=11.8,height=10.2,dpi=300,bg="white")
}
ggsave(file.path(out,"Figure6-publication-preview.png"),figure,width=11.8,height=10.2,dpi=120,bg="white")
write_table(statistics,"factorial_anova.csv");write_table(means,"substrate_marginal_means.csv")
write_table(contrasts,"substrate_pairwise_Tukey.csv");write_table(bind_rows(conditional_rows),"conditional_substrate_pairwise_Tukey.csv")
write_table(raw,"plotted_samples.csv");write_table(bind_rows(bracket_rows),"displayed_brackets.csv")
write_table(bind_rows(annotation_rows),"displayed_ANOVA_labels.csv");write_table(bind_rows(check_rows),"design_checks.csv")
write_table(bind_rows(plot_checks),"plot_validation.csv")
saveRDS(models,file.path(data_out,"figure6_models.rds"));saveRDS(figure,file.path(data_out,"figure6_plot.rds"))
write_json(list(outcomes=4,N=spec$N,blocks=6,random_block_retained=TRUE,fixed_block_columns=0,
  primary_model="response ~ inoculum * substrate + (1 | block)",
  corner_labels="S: BH q < 0.05; NS: BH q >= 0.05; adjustment across 4 outcomes separately for each term",
  displayed_substrate_brackets=6,bracket_p_values="Tukey within 3 marginal substrate comparisons per outcome",
  marginal_means="Equal inoculum weights; log responses back-transformed to geometric means",
  fourth_panel="Biomass-to-respiration ratio; CUE proxy, not a direct physiological CUE measure",
  treatment_labels=origins,inoculum_shapes=as.list(shapes),substrate_colours=as.list(palette),font_sizes_pt=font_sizes),
  file.path(out,"validation.json"),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,"R_session_info.txt"))
cat("Figure 6 saved with 12 ANOVA labels, inoculum shapes, marginal means, and six Tukey brackets.\n")
