#!/usr/bin/env Rscript
# ADDITIVE outputs only: same seven-category temporal stacks for master/selected,
# complete enriched metadata for all 33,831 and fixed 3,434 selected isolates.
# The original Rmd's statistics, plots, selection and source logic are unchanged.
source(file.path(Sys.getenv("SEN_ROOT", unset=getwd()), "analysis/part_a/lib/common.R"))

root <- Sys.getenv("SEN_ROOT", unset=getwd())
raw_file <- sen_path("results/part_a/master_33831/00_parsed_raw/SEN_33831_PARSED_51COL.csv")
rmd_file <- sen_path("results/part_a/full_master_original/SEN_33831_original_Rmd_derived.csv")
selected_file <- sen_path("results/part_a/cohort/Finalupdated_SEN_Genomes1_QC_merged.csv")
out <- ensure_dir(sen_path("results/part_a/reconciled_metadata"))
plot_out <- ensure_dir(sen_path("results/part_a/temporal_top5"))

raw <- readr::read_csv(raw_file, col_types=readr::cols(.default=readr::col_character()),
                       show_col_types=FALSE, progress=FALSE)
rmd <- readr::read_csv(rmd_file, col_types=readr::cols(.default=readr::col_character()),
                       show_col_types=FALSE, progress=FALSE)
sel <- readr::read_csv(selected_file, col_types=readr::cols(.default=readr::col_character()),
                       show_col_types=FALSE, progress=FALSE)
require_cols <- function(df, cols, label) {
  stopifnot(all(cols %in% names(df)))
  if (anyDuplicated(df[,cols,drop=FALSE])) stop(label, ": duplicate Run+Isolate")
}
require_cols(raw,c("Run","Isolate"),"Raw master")
require_cols(rmd,c("Run","Isolate"),"Original Rmd derived")
require_cols(sel,c("Run","Isolate"),"QC merged selection")
stopifnot(nrow(raw)==33831L, ncol(raw)==51L,
          nrow(rmd)==33831L, nrow(sel)==3434L,
          setequal(raw$Run,rmd$Run),
          all(sel$Run %in% raw$Run),
          setequal(paste(raw$Run,raw$Isolate,sep="\t"),
                   paste(rmd$Run,rmd$Isolate,sep="\t")))
if (any(names(rmd)[-c(1,2)] %in% names(raw)))
  stop("Original Rmd derived export unexpectedly includes raw master columns")
names(rmd)[!names(rmd) %in% c("Run","Isolate")] <-
  paste0("Rmd_",names(rmd)[!names(rmd) %in% c("Run","Isolate")])
full <- dplyr::left_join(raw,rmd,by=c("Run","Isolate"))
stopifnot(nrow(full)==33831L, all(!is.na(full$Rmd_Final_source)),
          all(full$Rmd_Final_source %in%
              c("Human","Chicken-US","Other-US","International")))
src <- full %>% dplyr::count(Rmd_Final_source,name="n") %>%
  dplyr::arrange(Rmd_Final_source)
observed <- setNames(src$n,src$Rmd_Final_source)
expect <- c("Human"=26541L,"Chicken-US"=6637L,"Other-US"=582L,
            "International"=71L)
stopifnot(identical(as.integer(observed[names(expect)]),as.integer(expect)))

# Same curated marker-token extraction for the master and selected population.
# IMPORTANT: each genome occupies exactly ONE bucket, including multi-marker
# genomes (priority in major list order); multi-marker status is exported.
major <- c("gyrA_D87Y","gyrA_D87N","gyrA_S83Y","gyrA_D87G","qnrB19")
levels7 <- c("QS/FQS",major,"Others")
tokens <- extract_qr_core_markers(full[["AMR genotypes core"]])
bucket <- function(t) {
  if (!length(t)) return("QS/FQS")
  hits <- which(tolower(major) %in% tolower(t))
  if (length(hits)) major[hits[1L]] else "Others"
}
year <- as.integer(derive_collection_year(full[["Collection date"]]))
period <- dplyr::case_when(year < 2018L ~ "Pre-2018",
  year >= 2018L & year <= 2020L ~ "2018-2020",
  year >= 2021L & year <= 2025L ~ "2021-2025",
  TRUE ~ NA_character_)
full <- full %>% dplyr::mutate(
  Audit_Collection_year=year,
  Audit_Temporal_block=period,
  Audit_Source=Rmd_Final_source,
  Audit_Animal_source=dplyr::if_else(Rmd_Final_source %in%
    c("International","Other-US"),Rmd_Other_source_group,NA_character_),
  Audit_QR_marker_count=lengths(tokens),
  Audit_QR_marker_tokens=vapply(tokens,function(x) paste(sort(unique(x)),collapse=";"),
                               character(1)),
  Audit_QR_status=dplyr::if_else(Audit_QR_marker_count>0L,
                               "QR/FQR marker-positive","QS/FQS (no marker)"),
  Audit_Top5_bucket=vapply(tokens,bucket,character(1)),
  Audit_Multi_marker=Audit_QR_marker_count>1L,
  Audit_Selected_pre_QC=Run %in% sel$Run)
stopifnot(sum(full$Audit_Selected_pre_QC)==3434L,
          sum(!is.na(full$Audit_Temporal_block))==33508L,
          identical(sum(full$Audit_QR_marker_count>0L),sum(full$Audit_Top5_bucket!="QS/FQS")))

# Prefix every selected QC value when appending to the full master. The
# selected original 78 columns remain unchanged in the selected export.
qc_cols <- c("QR_genotype_archived","FQR_genotype_archived","QR_genotype",
             "FQR_genotype","QR_status","QR_marker_count","QR_markers",
             "QR_corrected","Clade","Final_phylogeny","AUDIT")
stopifnot(all(qc_cols %in% names(sel)))
qc <- sel %>% dplyr::select(Run,Isolate,dplyr::all_of(qc_cols))
names(qc)[-c(1,2)] <- paste0("QC_",names(qc)[-c(1,2)])
full <- dplyr::left_join(full,qc,by=c("Run","Isolate"))
full$Audit_Analyzed_post_QC <- !is.na(full$QC_Final_phylogeny) &
                               full$QC_Final_phylogeny=="Included"
stopifnot(nrow(full)==33831L,sum(full$Audit_Analyzed_post_QC)==3306L,
          sum(!is.na(full$QC_Final_phylogeny))==3434L)

# Shared metadata columns are taken FROM the full master so corresponding
# selected and full rows reconcile exactly, without altering archived selected
# columns or correcting anything besides the three already curated QR edits.
extra <- names(full)[grepl("^(Rmd_|Audit_)", names(full))]
sel_extra <- full %>% dplyr::select(Run,Isolate,dplyr::all_of(extra))
stopifnot(!any(extra %in% names(sel)))
sel_aug <- dplyr::left_join(sel,sel_extra,by=c("Run","Isolate"))
stopifnot(nrow(sel_aug)==3434L,ncol(sel)==78L,
          all(sel_aug$Audit_Selected_pre_QC),
          sum(sel_aug$Audit_Analyzed_post_QC)==3306L,
          sum(sel_aug$QR_corrected=="YES")==3L,
          all(sel_aug$Run==sel$Run))
if (!all(vapply(names(sel),function(n) identical(as.character(sel[[n]]),
                                                 as.character(sel_aug[[n]])),
                logical(1))))
  stop("Selected original/QC fields changed in the enriched export")

# Mismatches are reported at isolate level: never silently overwrite prior
# historical selection metadata when regenerated master classifications differ.
reconciliation <- sel_aug %>%
  dplyr::transmute(
    Run,Isolate,
    Selected_Source=Source,Master_Source=Audit_Source,
    Source_mismatch=!is.na(Source) & Source != Audit_Source,
    Selected_year=Collection_year,Master_year=as.character(Audit_Collection_year),
    Year_mismatch=!is.na(Collection_year) &
      Collection_year!=as.character(Audit_Collection_year),
    Selected_QR_genotype=QR_genotype,
    Master_QR_genotype=as.character(as.integer(Audit_QR_marker_count>0L)),
    QR_mismatch=QR_genotype!=as.character(as.integer(Audit_QR_marker_count>0L)),
    QC_Final_phylogeny=Final_phylogeny,QC_Clade=Clade,
    QC_QR_corrected=QR_corrected) %>%
  dplyr::filter(Source_mismatch | Year_mismatch | QR_mismatch)
readr::write_csv(reconciliation,file.path(out,"SEN_3434_master_selected_field_discrepancies.csv"),na="")
if (nrow(reconciliation))
  message("[REVIEW] ",nrow(reconciliation),
          " selected-versus-master field discrepancies exported; none overwritten")

full_path <- file.path(out,"SEN_33831_full_enriched_metadata.csv")
sel_path <- file.path(out,"SEN_3434_selected_enriched_QC_metadata.csv")
readr::write_csv(full,full_path,na="")
readr::write_csv(sel_aug,sel_path,na="")
summary <- tibble::tibble(
  Metric=c("Full genome rows","Original full metadata columns","Master Rmd derived columns",
    "Selected QC original columns","Selected genome rows","Final phylogeny included",
    "Final phylogeny excluded","QR genotype corrections",
    "Full master dated genomes","Full master missing collection year",
    "Selected-master source/year/QR discrepancies"),
  Value=c(nrow(full),ncol(raw),ncol(rmd)-2L,ncol(sel),nrow(sel_aug),
    sum(full$Audit_Analyzed_post_QC),
    sum(full$Audit_Selected_pre_QC)-sum(full$Audit_Analyzed_post_QC),
    sum(sel$QR_corrected=="YES"),sum(!is.na(year)),sum(is.na(year)),
    nrow(reconciliation)))
readr::write_csv(summary,file.path(out,"SEN_metadata_reconciliation_summary.csv"))

# Calculate mutually exclusive seven-category stacks. The same priority
# categories, dates, years, denominator logic and palette apply to both sets.
counts_for <- function(df, cohort) {
  df %>%
    dplyr::filter(!is.na(Audit_Temporal_block)) %>%
    dplyr::mutate(Period=factor(Audit_Temporal_block,
      levels=c("Pre-2018","2018-2020","2021-2025")),
      Category=factor(Audit_Top5_bucket,levels=levels7)) %>%
    dplyr::count(Period,Category,name="Count") %>%
    tidyr::complete(
      Period=factor(c("Pre-2018","2018-2020","2021-2025"),
                    levels=c("Pre-2018","2018-2020","2021-2025")),
      Category=factor(levels7,levels=levels7),
      fill=list(Count=0L)) %>%
    dplyr::group_by(Period) %>%
    dplyr::mutate(Total_N=sum(Count),
                  Percent=ifelse(Total_N>0,100*Count/Total_N,NA_real_)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Cohort=cohort)
}
total_counts <- counts_for(full,"Total available 33,831-genome master")
selected_counts <- counts_for(sel_aug,"Selected 3,434-genome cohort")
stopifnot(sum(total_counts$Count)==33508L,
          sum(selected_counts$Count)==3434L,
          all(abs(total_counts %>% dplyr::group_by(Period) %>%
            dplyr::summarise(p=sum(Percent),.groups="drop") %>%
            dplyr::pull(p)-100)<1e-9),
          all(abs(selected_counts %>% dplyr::group_by(Period) %>%
            dplyr::summarise(p=sum(Percent),.groups="drop") %>%
            dplyr::pull(p)-100)<1e-9))
readr::write_csv(dplyr::bind_rows(total_counts,selected_counts),
                 file.path(plot_out,"SEN_Total_Selected_TemporalBlocks_Top5_QSFQS_data.csv"))
multi <- full %>% dplyr::filter(Audit_Multi_marker) %>%
  dplyr::select(Run,Isolate,Audit_Collection_year,Audit_Source,
                Audit_QR_marker_count,Audit_QR_marker_tokens,Audit_Top5_bucket,
                Audit_Selected_pre_QC,Audit_Analyzed_post_QC)
readr::write_csv(multi,file.path(plot_out,"SEN_master_multimarker_bucket_priority_audit.csv"))

palette <- c("QS/FQS"="#D9D9D9","gyrA_D87Y"="#0072B2",
             "gyrA_D87N"="#D55E00","gyrA_S83Y"="#009E73",
             "gyrA_D87G"="#CC79A7","qnrB19"="#E69F00",
             "Others"="#666666")
make_plot <- function(df, cohort_label, y_key) {
  pp <- ggplot2::ggplot(df,ggplot2::aes(x=Period,y=.data[[y_key]],fill=Category)) +
    ggplot2::geom_col(width=0.72) +
    ggplot2::scale_fill_manual(values=palette,drop=FALSE) +
    ggplot2::labs(
      title=cohort_label,
      subtitle="One mutually exclusive bucket per genome; QR/FQR marker priority: top five, then Others",
      x="Collection period",
      y=if (y_key=="Percent") "Percentage of genomes (%)" else "Number of genomes",
      fill="QR/FQR genotype") +
    ggplot2::theme_classic(base_size=13) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(face="bold"),
                   legend.position="right")
  if (y_key=="Percent") pp <- pp +
    ggplot2::scale_y_continuous(labels=function(y) paste0(y,"%"),limits=c(0,105))
  else pp <- pp + ggplot2::scale_y_continuous(labels=scales::label_comma())
  pp
}
save_plots <- function(df,prefix,label) {
  for (metric in c("Count","Percent")) {
    p <- make_plot(df,label,metric)
    ggplot2::ggsave(file.path(plot_out,paste0(prefix,"_",metric,".png")),p,
                    width=10,height=6.2,dpi=600,bg="white")
    ggplot2::ggsave(file.path(plot_out,paste0(prefix,"_",metric,".tif")),p,
                    width=10,height=6.2,dpi=600,compression="lzw",bg="white")
    ggplot2::ggsave(file.path(plot_out,paste0(prefix,"_",metric,".pdf")),p,
                    width=10,height=6.2,bg="white")
  }
}
save_plots(total_counts,"SEN_Total_33831_TemporalBlocks_Top5_QSFQS",
           "Total available SEN genomes")
save_plots(selected_counts,"SEN_Selected_3434_TemporalBlocks_Top5_QSFQS",
           "Selected SEN genomes")

cat("[PASS] Full enriched metadata: ",nrow(full)," rows, ",
    ncol(full)," columns, all original 51 retained\n",sep="")
cat("[PASS] Selected enriched metadata: ",nrow(sel_aug)," rows, ",
    ncol(sel_aug)," columns, all archived 69 + 9 QC retained\n",sep="")
cat("[PASS] Seven-category full master temporal stacks: 33,508 dated genomes; ",
    "selected stacks: 3,434 genomes; marker multi-hits assigned once\n",sep="")
cat("[OUTPUT] ",full_path,"\n[OUTPUT] ",sel_path,"\n",sep="")
