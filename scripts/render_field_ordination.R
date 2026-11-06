suppressPackageStartupMessages({library(ggplot2);library(ggrepel);library(dplyr)})
script<-normalizePath(sub("^--file=","",commandArgs(FALSE)[grep("^--file=",commandArgs(FALSE))]))
base<-dirname(dirname(script));field_dir<-normalizePath(commandArgs(TRUE)[1]);out_dir<-field_dir
plot_env<-new.env(parent=globalenv());sys.source(file.path(base,"scripts","field_ordination_style.R"),envir=plot_env)
plot_env$colours<-c("S+"="#E3AC18","UV-"="#77BCE8","S-"="#555555")
plot_env$treat_labels<-c("S+"="Full sunlight","UV-"="Reduced UV","S-"="Reduced UV and blue-green")
scores <- read.csv(file.path(field_dir,"pcoa_scores.csv"),stringsAsFactors=FALSE)
scores$Treatment <- factor(scores$Treatment,levels=c("S+","UV-","S-"))
loadings <- read.csv(file.path(field_dir,"pcoa_loadings.csv"),stringsAsFactors=FALSE)
global <- read.csv(file.path(field_dir,"global_statistics.csv"))
variance <- as.numeric(global[1,c("PCoA1","PCoA2")])
fig <- plot_env$make_pcoa_figure(scores,loadings,variance,global,"field")
colored_build <- ggplot_build(fig)

# Fixed grey aesthetics remove process-category encoding from arrows and labels.
grey <- "#3F3F3F"
fig$layers[[5]]$mapping <- aes(x=0,y=0,xend=xend,yend=yend)
fig$layers[[5]]$aes_params$colour <- grey
fig$layers[[6]]$mapping <- aes(x=xend,y=yend,label=display_label)
fig$layers[[6]]$aes_params$colour <- grey
# Match the labels in the author-controlled publication Figure 2 PDF.
publication_labels <- c("S+"="Full sunlight","UV-"="Reduced UV",
                        "S-"="Reduced UV and BG")
fill_scale <- fig$scales$get_scales("fill")
fill_scale$labels <- publication_labels
fill_scale$name <- NULL
fig <- fig+guides(colour="none",
  fill=guide_legend(position="bottom",nrow=1,order=1,
                    override.aes=list(size=4.5)))+
  theme(legend.position="bottom")

# Only these three labels retain leaders. A zero point radius and negligible
# padding attach their leaders to the actual loading tips, without visible gaps.
leader_groups <- c("Alpha-glucan and fructan processing",
                   "Hemicellulose degradation","Enzymatic ROS defense")
label_data <- fig$layers[[6]]$data
has_leader <- label_data$functional_group %in% leader_groups
photo_label <- which(label_data$functional_group=="Photoreactivation and nucleotide-excision repair")
ros_label <- which(label_data$functional_group=="Enzymatic ROS defense")
dna_label <- which(label_data$functional_group=="General DNA-damage repair and tolerance")
label_data$dx[photo_label] <- label_data$dx[photo_label]-.9
label_data$dy[photo_label] <- label_data$dy[photo_label]+.2
label_data$dx[ros_label] <- label_data$dx[ros_label]-.5
label_data$dx[dna_label] <- label_data$dx[dna_label]+.3
fig$layers[[6]]$data <- label_data
fig$layers[[6]]$position <- ggrepel::position_nudge_repel(x=label_data$dx,y=label_data$dy)
fig$layers[[6]]$aes_params$segment.colour <- ifelse(has_leader,grey,NA_character_)
point_sizes <- ifelse(nzchar(label_data$display_label),4.25,6.5)
point_sizes[has_leader] <- 0
fig$layers[[6]]$aes_params$point.size <- point_sizes
fig$layers[[6]]$geom_params$point.padding <- grid::unit(1e-6,"lines")
fig$layers[[6]]$geom_params$min.segment.length <- grid::unit(0,"lines")
grey_build <- ggplot_build(fig)

# Geometry, treatment styling, loading lengths, and variance labels are retained.
stopifnot(identical(colored_build$data[1:4],grey_build$data[1:4]),
          identical(colored_build$data[[5]][,c("x","y","xend","yend")],
                    grey_build$data[[5]][,c("x","y","xend","yend")]),
          all(grey_build$data[[5]]$colour==grey),
          all(grey_build$data[[6]]$colour==grey),
          all(grey_build$data[[6]]$size==13/ggplot2::.pt),
          sum(!is.na(grey_build$data[[6]]$segment.colour))==3,
          setequal(label_data$functional_group[has_leader],leader_groups),
          identical(as.character(grey_build$plot$scales$get_scales("fill")$get_labels()),
                    unname(publication_labels)),
          is.null(fig$labels$title),length(fig$layers)==6)

for(ext in c("pdf","png"))ggsave(file.path(out_dir,paste0("pcoa_publication.",ext)),fig,width=10.5,height=8,dpi=300,bg="white")
