#!/usr/bin/env Rscript
# Original Part A2 master-population workflow: 33,831 records, NOT the 3,434 sample.
# Stage 0 must run first and writes the only raw parsed 51-column snapshot.
source(file.path(Sys.getenv("SEN_ROOT", unset=getwd()), "analysis/part_a/lib/common.R"))
root <- Sys.getenv("SEN_ROOT", unset=getwd())
raw_path <- sen_path("results/part_a/master_33831/00_parsed_raw/SEN_33831_PARSED_51COL.csv")
if (!file.exists(raw_path)) stop("Run analysis/part_a/00_parse_master_33831.R first.")
out <- ensure_dir(sen_path("results/part_a/master_33831"))
raw <- readr::read_csv(raw_path,col_types=readr::cols(.default=readr::col_character()),show_col_types=FALSE,progress=FALSE)
stopifnot(nrow(raw)==33831L,ncol(raw)==51L,dplyr::n_distinct(raw$Run)==33831L)
needed <- c("Run","Isolate","Strain","Collection date","Source type","Isolation source",
 "Host","IFSAC category","Collected by","Isolate identifiers","Food origin","AMR genotypes core")
assert_columns(raw,needed,"parsed master")

# Canonical functions from Part A2, shared with the selected-cohort reconstruction.
primary <- classify_primary_source(raw[["Source type"]],raw[["Isolation source"]],
 raw[["Host"]],raw[["IFSAC category"]],raw[["Collected by"]])
origin <- derive_origin_and_source(primary,raw[["Strain"]],
 raw[["Isolate identifiers"]],raw[["Food origin"]])
master <- raw %>% dplyr::mutate(
 Collection_year=as.integer(derive_collection_year(raw[["Collection date"]])),
 Source=origin$Source,Origin_group=origin$Origin_group,Source_group=origin$Source_group,
 QR_genotype=as.integer(derive_qr_genotype_core(raw[["AMR genotypes core"]])),
 QR_markers=collapse_qr_core_markers(raw[["AMR genotypes core"]]),
 QR_status=ifelse(QR_genotype==1L,"QR/FQR marker-positive","QR/FQR marker-negative"),
 Temporal_eligible_1990_2025=!is.na(Collection_year) &
       Collection_year>=1990 & Collection_year<=2025
)
# Selected and post-QC markers are attached ONLY for audits and coverage analyses.
keys <- readr::read_csv(sen_path("data/part_a/selection/Selected_updated_Isolate.csv"),
 show_col_types=FALSE,col_types=readr::cols(.default=readr::col_character()))
keys <- keys %>% dplyr::mutate(dplyr::across(c(Strain,Isolate),~stringr::str_squish(as.character(.x))))
master <- master %>% dplyr::mutate(dplyr::across(c(Strain,Isolate),~stringr::str_squish(as.character(.x))))
selected_rows <- master %>% dplyr::semi_join(keys,by=c("Strain","Isolate"))
stopifnot(nrow(selected_rows)==3434L,dplyr::n_distinct(selected_rows$Run)==3434L)
clades <- readr::read_tsv(sen_path("metadata/final_clade_metadata.tsv"),show_col_types=FALSE)
final_runs <- unique(as.character(clades$Sample[clades$Sample!="Reference"]))
stopifnot(length(final_runs)==3306L,all(final_runs %in% selected_rows$Run))
master$Selection_status <- ifelse(master$Run %in% selected_rows$Run,"Selected","Not selected")
master$Final_phylogeny_status <- ifelse(master$Run %in% final_runs,"Analyzed",
 ifelse(master$Run %in% selected_rows$Run,"Selected, excluded","Not selected"))
stopifnot(sum(master$Selection_status=="Selected")==3434L,
 sum(master$Final_phylogeny_status=="Analyzed")==3306L)

# Optional cross-check against the regenerated selected cohort if Part A ran earlier.
rebuilt_path <- sen_path("results/part_a/cohort/Finalupdated_SEN_Genomes_rebuilt.csv")
if (file.exists(rebuilt_path)) {
 sel <- readr::read_csv(rebuilt_path,show_col_types=FALSE,progress=FALSE)
 compare <- master %>% dplyr::filter(Run %in% sel$Run) %>%
   dplyr::select(Run,Source,Collection_year,QR_genotype) %>%
   dplyr::inner_join(sel %>% dplyr::transmute(Run,sel_source=as.character(Source),
     sel_year=as.integer(Collection_year),sel_qr=as.integer(QR_genotype)),by="Run")
 bad <- compare %>% dplyr::filter(Source!=sel_source |
   xor(is.na(Collection_year),is.na(sel_year)) |
   dplyr::coalesce(Collection_year!=sel_year,FALSE) | QR_genotype!=sel_qr)
 readr::write_csv(bad,file.path(out,"master_vs_rebuilt_discrepancies.csv"))
 if(nrow(compare)!=3434L || nrow(bad)>0) stop("Master/selected audit failed: ",
 nrow(compare)," matched, ",nrow(bad)," disagreeing; inspect master_vs_rebuilt_discrepancies.csv")
}
audit <- master %>% dplyr::group_by(Source) %>% dplyr::summarise(
 N_available=dplyr::n(),N_missing_year=sum(is.na(Collection_year)),
 N_outside_1990_2025=sum(!is.na(Collection_year)&!Temporal_eligible_1990_2025),
 N_temporal_1990_2025=sum(Temporal_eligible_1990_2025),
 N_selected=sum(Selection_status=="Selected"),
 N_analyzed=sum(Final_phylogeny_status=="Analyzed"),.groups="drop")
stopifnot(sum(audit$N_available)==33831L,sum(audit$N_selected)==3434L,
 sum(audit$N_analyzed)==3306L)
readr::write_csv(audit,file.path(out,"master_source_temporal_qc_audit.csv"))
readr::write_csv(master,file.path(out,"SEN_33831_DERIVED_ANALYSIS_READY.csv"),na="")
if (!requireNamespace("writexl",quietly=TRUE)) stop("Install r-writexl into sen_stats")
writexl::write_xlsx(list(Derived_master_33831=master,QC_summary=audit),
 file.path(out,"SEN_33831_DERIVED_ANALYSIS_READY.xlsx"))

# Status prevalence: denominator is the available MASTER genomes per source and year.
valid <- master %>% dplyr::filter(Temporal_eligible_1990_2025) %>%
 dplyr::mutate(year_plot=ifelse(Collection_year<=2000,"≤2000",as.character(Collection_year)),
 year5=dplyr::case_when(Collection_year<=2000~"≤2000",
 Collection_year<=2005~"2001-2005",Collection_year<=2010~"2006-2010",
 Collection_year<=2015~"2011-2015",Collection_year<=2020~"2016-2020",
 TRUE~"2021-2025"))
annual <- valid %>% dplyr::count(Source,Collection_year,QR_status,name="n") %>%
 dplyr::group_by(Source,Collection_year) %>%
 dplyr::mutate(N_year=sum(n),pct=100*n/N_year) %>% dplyr::ungroup()
readr::write_csv(annual,file.path(out,"master_QRFQR_by_source_and_year.csv"))
stack <- valid %>% dplyr::count(Source,year_plot,QR_status,name="n") %>%
 dplyr::group_by(Source,year_plot) %>%
 dplyr::mutate(N_year=sum(n),pct=100*n/N_year) %>% dplyr::ungroup()
readr::write_csv(stack,file.path(out,"master_QRFQR_by_source_and_plot_year.csv"))
strata <- valid %>% dplyr::group_by(Source,Collection_year,QR_status) %>%
 dplyr::summarise(N_available=dplyr::n(),N_selected=sum(Selection_status=="Selected"),
 N_analyzed=sum(Final_phylogeny_status=="Analyzed"),.groups="drop") %>%
 dplyr::mutate(pct_selected=100*N_selected/N_available,
 pct_analyzed=100*N_analyzed/N_available,
 weight_available_over_analyzed=ifelse(N_analyzed>0,N_available/N_analyzed,NA_real_))
readr::write_csv(strata,file.path(out,"master_selection_weights_year_source_marker.csv"))
src_n <- valid %>% dplyr::count(Source,name="N_temporal")
src_label <- stats::setNames(paste0(src_n$Source," (N=",format(src_n$N_temporal,big.mark=",",trim=TRUE),")"),src_n$Source)
valid$Source_label <- factor(src_label[valid$Source],levels=src_label[unique(valid$Source)])
# Exact production token extraction, one marker per genome (no duplicate inflation).
tokens <- master %>% dplyr::transmute(Run,Source,Collection_year,
 token=extract_qr_core_markers(.data[["AMR genotypes core"]])) %>%
 tidyr::unnest_longer(token,values_to="token") %>%
 dplyr::filter(!is.na(token),token!="") %>%
 dplyr::mutate(marker=pretty_token(tolower(token))) %>%
 dplyr::distinct(Run,Source,Collection_year,marker)
marker_frequency <- tokens %>% dplyr::count(marker,name="N_genomes",sort=TRUE) %>%
 dplyr::mutate(pct_all_available=100*N_genomes/33831)
readr::write_csv(marker_frequency,file.path(out,"master_marker_frequency.csv"))
burden <- tokens %>% dplyr::count(Run,name="n_markers")
readr::write_csv(burden,file.path(out,"master_marker_burden_per_isolate.csv"))
top5 <- head(marker_frequency$marker,5)
marker_plot <- tokens %>% dplyr::filter(!is.na(Collection_year),
 Collection_year>=1990,Collection_year<=2025) %>%
 dplyr::mutate(marker_group=ifelse(marker %in% top5,marker,"others"),
 year_plot=ifelse(Collection_year<=2000,"≤2000",as.character(Collection_year)),
 year5=dplyr::case_when(Collection_year<=2000~"≤2000",
 Collection_year<=2005~"2001-2005",Collection_year<=2010~"2006-2010",
 Collection_year<=2015~"2011-2015",Collection_year<=2020~"2016-2020",TRUE~"2021-2025")) %>%
 dplyr::distinct(Run,Source,year_plot,year5,marker_group)
plot_dots <- function(bin,filename) {
 den <- valid %>% dplyr::count(Source,.data[[bin]],name="N")
 names(den)[2] <- "period"
 dots <- marker_plot %>% dplyr::count(Source,.data[[bin]],marker_group,name="n")
 names(dots)[2] <- "period"
 dots <- dots %>% dplyr::left_join(den,by=c("Source","period")) %>%
 dplyr::mutate(pct=100*n/N,Source_label=src_label[Source])
 readr::write_csv(dots,file.path(out,paste0(filename,".csv")))
 p <- ggplot2::ggplot(dots,ggplot2::aes(x=period,y=marker_group))+
 ggplot2::geom_point(ggplot2::aes(size=pct,color=pct),alpha=.85)+
 ggplot2::facet_wrap(~Source_label,ncol=1)+
 ggplot2::labs(x="Collection year",y="QR/FQR marker (top five + others)",
 size="Year-specific %",color="Year-specific %")+
 ggplot2::theme_bw()+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=.5))
 p
}
p_dot_year <- plot_dots("year_plot","master_dotplot_QRFQR_markers_by_year")
p_dot_5 <- plot_dots("year5","master_dotplot_QRFQR_markers_by_5year")
stack$Source_label <- src_label[stack$Source]
p_pct <- ggplot2::ggplot(stack,ggplot2::aes(x=year_plot,y=pct,fill=QR_status))+
 ggplot2::geom_col()+ggplot2::facet_wrap(~Source_label,ncol=1)+
 ggplot2::labs(x="Collection year",y="Genomes (%)",fill="")+
 ggplot2::theme_bw()+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90,hjust=1))
p_count <- ggplot2::ggplot(stack,ggplot2::aes(x=year_plot,y=n,fill=QR_status))+
 ggplot2::geom_col()+ggplot2::facet_wrap(~Source_label,ncol=1,scales="free_y")+
 ggplot2::labs(x="Collection year",y="Genomes (n)",fill="")+
 ggplot2::theme_bw()+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90,hjust=1))
by_year <- valid %>% dplyr::group_by(Source,Collection_year) %>%
 dplyr::summarise(Available=dplyr::n(),Selected=sum(Selection_status=="Selected"),
 Analyzed=sum(Final_phylogeny_status=="Analyzed"),.groups="drop")
p_cohorts <- tidyr::pivot_longer(by_year,c(Available,Selected,Analyzed),
 names_to="Cohort",values_to="n") %>%
 ggplot2::ggplot(ggplot2::aes(Collection_year,n,color=Cohort))+
 ggplot2::geom_line()+ggplot2::facet_wrap(~Source,ncol=1,scales="free_y")+
 ggplot2::labs(x="Collection year",y="Genomes (n)")+ggplot2::theme_bw()
p_coverage <- by_year %>% dplyr::mutate(
 pct_selected=100*Selected/Available,pct_analyzed=100*Analyzed/Available) %>%
 tidyr::pivot_longer(c(pct_selected,pct_analyzed),names_to="Stage",values_to="Percent") %>%
 ggplot2::ggplot(ggplot2::aes(Collection_year,Percent,color=Stage))+
 ggplot2::geom_line()+ggplot2::geom_point(size=.8)+
 ggplot2::facet_wrap(~Source,ncol=1)+
 ggplot2::labs(x="Collection year",y="Share of available genomes (%)")+ggplot2::theme_bw()
figs <- list(master_dotplot_QRFQR_markers_by_year=p_dot_year,
 master_dotplot_QRFQR_markers_by_5year=p_dot_5,
 master_QRFQR_percent_source_year=p_pct,
 master_QRFQR_counts_source_year=p_count,
 master_available_selected_analyzed_year=p_cohorts,
 master_sampling_coverage_year=p_coverage)
for(nm in names(figs)) {
 ggplot2::ggsave(file.path(out,paste0(nm,".png")),figs[[nm]],width=12,height=8,dpi=600)
 ggplot2::ggsave(file.path(out,paste0(nm,".pdf")),figs[[nm]],width=12,height=8)
}
cat("[PASS] Master 33,831: derived XLSX/CSV, QC, six figures and sampling weights at ",out,"\n",sep="")
print(audit)
