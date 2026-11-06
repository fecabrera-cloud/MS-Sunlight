# Shared factorial inference for microcosm community distances.
sunlight_levels <- c("Full sunlight","Reduced UV","Reduced UV and blue-green")
full_treatment <- function(x) {
  codes <- c("S+"="Full sunlight","UV-"="Reduced UV","S-"="Reduced UV and blue-green",
    "SusS+"="Full sunlight","SusUV-"="Reduced UV","SusS-"="Reduced UV and blue-green")
  z <- sub("^(Sus|Inoc)","",as.character(x));z[z=="R+"]<-"S+";z[z=="R-"]<-"S-"
  ii <- z%in%names(codes);z[ii] <- codes[z[ii]]
  stopifnot(all(z%in%sunlight_levels));factor(z,levels=sunlight_levels)
}
factorial_permanova <- function(D,md,P,out,seed) {
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  write_table <- function(x,n)write.csv(x,file.path(out,n),row.names=FALSE,na="")
  md$inoculum <- full_treatment(md$inoculum);md$substrate <- full_treatment(md$substrate);md$block <- factor(md$block)
  stopifnot(identical(attr(D,"Labels"),md$sample_id),nrow(P)==9999,ncol(P)==nrow(md),
    qr(model.matrix(~inoculum*substrate,md))$rank==9,
    all(vapply(seq_len(nrow(P)),function(i)all(md$block[P[i,]]==md$block),logical(1))))
  fit <- vegan::adonis2(D~inoculum*substrate,data=md,permutations=P,by="terms")
  a <- data.frame(term=rownames(fit),df=fit$Df,SumOfSqs=fit$SumOfSqs,R2=fit$R2,F=fit$F,p=fit$`Pr(>F)`)
  write_table(a,"sequential_permanova.csv")
  write_table(data.frame(component=c("Inoculum","Substrate","Interaction","Residual"),
    df=fit$Df[1:4],inertia=fit$SumOfSqs[1:4],proportion=fit$R2[1:4],percent=100*fit$R2[1:4]),"sequential_variance_partition.csv")
  dispersions <- objects <- list()
  for(key in c("inoculum","substrate","cell")) {
    group <- if(key=="cell")interaction(md$inoculum,md$substrate,drop=TRUE) else md[[key]]
    bd <- vegan::betadisper(D,group);test <- vegan::permutest(bd,permutations=P)
    dispersions[[key]] <- data.frame(grouping=key,F=test$tab$F[1],p=test$tab$`Pr(>F)`[1],n_permutations=nrow(P),minimum_group_N=min(table(group)))
    objects[[key]] <- list(model=bd,test=test)
  }
  disp <- do.call(rbind,dispersions);write_table(disp,"dispersion_tests.csv")
  pp <- list();i <- 0L
  for(direction in c("substrate_within_inoculum","inoculum_within_substrate")) {
    byvar <- if(direction=="substrate_within_inoculum")"inoculum" else "substrate"
    focal <- if(byvar=="inoculum")"substrate" else "inoculum"
    for(level in sunlight_levels)for(pair in combn(sunlight_levels,2,simplify=FALSE)) {
      i <- i+1L;idx <- md[[byvar]]==level & md[[focal]]%in%pair
      mm <- droplevels(md[idx,,drop=FALSE]);dd <- as.dist(as.matrix(D)[idx,idx])
      ctrl <- permute::how(nperm=9999);permute::setBlocks(ctrl) <- mm$block
      set.seed(seed+i);perm <- suppressMessages(permute::shuffleSet(nrow(mm),nset=9999,control=ctrl))
      mm$comparison <- mm[[focal]];aa <- vegan::adonis2(dd~comparison,data=mm,permutations=perm)
      pp[[i]] <- data.frame(direction,conditioning_level=level,group1=pair[1],group2=pair[2],N=nrow(mm),
        F=aa$F[1],R2=aa$R2[1],p=aa$`Pr(>F)`[1],n_permutations=nrow(perm),minimum_attainable_p=1/(nrow(perm)+1))
    }
  }
  pp <- do.call(rbind,pp);pp$p_BH <- p.adjust(pp$p,"BH");write_table(pp,"conditional_pairwise_permanova.csv")
  write_table(md,"metadata.csv")
  saveRDS(list(fit=fit,distance_model=vegan::dbrda(D~inoculum*substrate,data=md),dispersion=objects,metadata=md),file.path(out,"models.rds"))
  list(anova=a,dispersion=disp,metadata=md)
}
