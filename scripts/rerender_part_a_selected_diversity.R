#!/usr/bin/env Rscript
# Figure-only renderer for the 3,434-selected-cohort diversity panel.
# Requires existing validated tables; does not rerun diversity statistics.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(patchwork)
})
root <- Sys.getenv("SEN_ROOT")
snap <- Sys.getenv("SEN_PARTA_SNAPSHOTS", unset=root)
source(file.path(root,"scripts/part_a_figure_label_helpers.R"))
src <- file.path(snap,"results/part_a/diversity")
out <- file.path(root,"results/part_a/mdpi_label_corrections")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
div_file <- file.path(src,"QR_diversity_q0_q1_q2_bySource_with95CI.csv")
comp_file <- file.path(src,"QR_composition_bySource.csv")
if (!file.exists(div_file) || !file.exists(comp_file))
  stop("Selected 3,434-genome diversity summaries are unavailable in ",src,
       ". Run workflow/A2_Genotype/05_diversity_composition.R first.")
div_sum <- readr::read_csv(div_file,show_col_types=FALSE)
comp <- readr::read_csv(comp_file,show_col_types=FALSE)
source_levels <- c("Human","Chicken-US","Other-US","International")
stopifnot(nrow(div_sum)==4L,
          setequal(as.character(div_sum$Source),source_levels),
          setequal(as.character(comp$Source),source_levels),
          all(c("q1","q1_lo","q1_hi","q2","q2_lo","q2_hi") %in% names(div_sum)),
          all(c("marker_group","p") %in% names(comp)))
div_sum$Source <- factor(as.character(div_sum$Source),levels=source_levels)
comp$Source <- factor(as.character(comp$Source),levels=source_levels)
comp$marker_group <- factor(comp$marker_group,levels=unique(as.character(comp$marker_group)))
source_display <- function(x) gsub("-", "\u2014", x, fixed=TRUE)

plot_div <- div_sum %>%
  tidyr::pivot_longer(c(q1,q2),names_to="metric",values_to="value") %>%
  dplyr::mutate(
    lo=dplyr::if_else(metric=="q1",q1_lo,q2_lo),
    hi=dplyr::if_else(metric=="q1",q1_hi,q2_hi),
    metric=dplyr::recode(metric,q1="Exponential Shannon (q1)",
                        q2="Inverse Simpson (q2)")
  )

p_div <- ggplot2::ggplot(plot_div,ggplot2::aes(Source,value)) +
  ggplot2::geom_point(size=2.8) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin=lo,ymax=hi),width=.15) +
  ggplot2::facet_wrap(~metric,nrow=1) +
  ggplot2::labs(x="Source",y="Effective number of QR/FQR categories") +
  ggplot2::scale_x_discrete(labels=source_display) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1))

p_comp <- ggplot2::ggplot(comp,ggplot2::aes(Source,p,fill=marker_group)) +
  ggplot2::geom_col() +
  ggplot2::scale_y_continuous(labels=scales::percent) +
  ggplot2::scale_x_discrete(labels=source_display) +
  ggplot2::scale_fill_discrete(
    labels=function(x) parse(text=sen_marker_plotmath(x))) +
  ggplot2::labs(x="Source",y="Normalized QR/FQR composition",fill="QR category") +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1))

combo <- p_div | p_comp
target <- file.path(out,"selected_3434_QR_diversity_and_composition.tif")
ggplot2::ggsave(target,combo,width=15,height=6,dpi=600,compression="lzw")
stopifnot(file.exists(target),file.info(target)$size>0)
cat("[PASS] Corrected selected 3,434-genome diversity TIFF: ",target,"\n",sep="")
