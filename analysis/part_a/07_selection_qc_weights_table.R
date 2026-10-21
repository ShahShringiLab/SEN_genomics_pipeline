#!/usr/bin/env Rscript
# Selection/QC accounting, independent of original Rmd figures and statistics.
# Full 33,831 master -> 3,434 pre-QC selected -> 3,306 post-QC final phylogeny.
source(file.path(Sys.getenv("SEN_ROOT", unset=getwd()), "analysis/part_a/lib/common.R"))
infile <- sen_path("results/part_a/master_33831/00_parsed_raw/SEN_33831_PARSED_51COL.csv")
if (!file.exists(infile)) stop("First run Stage 0 master parser.")
master <- readr::read_csv(infile, col_types=readr::cols(.default=readr::col_character()), show_col_types=FALSE)
assert_columns(master,c("Run","Isolate","Strain","Collection date","Source type","Isolation source",
 "Host","IFSAC category","Collected by","Isolate identifiers","Food origin","AMR genotypes core"))
stopifnot(nrow(master)==33831L,ncol(master)==51L,dplyr::n_distinct(master$Run)==33831L)
primary <- classify_primary_source(master[["Source type"]],master[["Isolation source"]],
 master[["Host"]],master[["IFSAC category"]],master[["Collected by"]])
groups <- derive_origin_and_source(primary, master$Strain,master[["Isolate identifiers"]],master[["Food origin"]])
master <- master %>% dplyr::mutate(
 Year_of_isolation=as.integer(derive_collection_year(.data[["Collection date"]])),
 Source=as.character(groups$Source),
 QR_FQR_Marker_status=ifelse(derive_qr_genotype_core(.data[["AMR genotypes core"]])==1L,
  "QR/FQR marker-positive","QS/FQS (no marker)"),
 Strain_join=stringr::str_squish(as.character(Strain)),
 Isolate_join=stringr::str_squish(as.character(Isolate)))
keys <- readr::read_csv(sen_path("data/part_a/selection/Selected_updated_Isolate.csv"),
 col_types=readr::cols(.default=readr::col_character()),show_col_types=FALSE)
assert_columns(keys,c("Strain","Isolate"))
keys <- keys %>% dplyr::transmute(Strain_join=stringr::str_squish(Strain),
 Isolate_join=stringr::str_squish(Isolate)) %>% dplyr::distinct()
stopifnot(nrow(keys)==3434L)
master <- master %>% dplyr::mutate(Selected_pre_QC=FALSE)
selected <- master %>% dplyr::semi_join(keys,by=c("Strain_join","Isolate_join"))
stopifnot(nrow(selected)==3434L,dplyr::n_distinct(selected$Run)==3434L)
master$Selected_pre_QC <- master$Run %in% selected$Run
clade <- readr::read_tsv(sen_path("metadata/final_clade_metadata.tsv"),show_col_types=FALSE)
assert_columns(clade,c("Sample","Lineage"))
analyzed_runs <- unique(as.character(clade$Sample[!is.na(clade$Sample) & clade$Sample!="Reference" &
 !is.na(clade$Lineage) & nzchar(clade$Lineage)]))
stopifnot(length(analyzed_runs)==3306L, all(analyzed_runs %in% selected$Run))
master$Analyzed_post_QC <- master$Run %in% analyzed_runs
stopifnot(sum(master$Selected_pre_QC)==3434L,sum(master$Analyzed_post_QC)==3306L,
 all(!master$Analyzed_post_QC | master$Selected_pre_QC))
# Do not manufacture year for missing-year genomes. Keep a separate explicit
# missing-year stratum so that all 33,831 master records are reconciled.
master$Year_of_isolation_label <- ifelse(is.na(master$Year_of_isolation),
 "Missing year",as.character(master$Year_of_isolation))
# One row per year x source x QR status, with every source-status combination
# represented for observed years, including zeros (weights undefined if n=0).
years <- unique(master$Year_of_isolation_label)
years <- c(sort(as.integer(setdiff(years,"Missing year"))),if("Missing year" %in% years)"Missing year")
sources <- c("Human","Chicken-US","Other-US","International")
statuses <- c("QS/FQS (no marker)","QR/FQR marker-positive")
stopifnot(all(master$Source %in% sources))
tab <- master %>%
 dplyr::group_by(Year_of_isolation_label,Source,QR_FQR_Marker_status) %>%
 dplyr::summarise(Available_genomes_N=dplyr::n(),
 Selected_genomes_pre_QC=sum(Selected_pre_QC),
 Analyzed_genomes_post_QC_n=sum(Analyzed_post_QC),.groups="drop") %>%
 tidyr::complete(Year_of_isolation_label=as.character(years),Source=sources,
 QR_FQR_Marker_status=statuses,
 fill=list(Available_genomes_N=0L,Selected_genomes_pre_QC=0L,Analyzed_genomes_post_QC_n=0L)) %>%
 dplyr::mutate(
 Pct_of_available_selected=ifelse(Available_genomes_N>0,
  100*Selected_genomes_pre_QC/Available_genomes_N,NA_real_),
 Pct_of_available_analyzed=ifelse(Available_genomes_N>0,
  100*Analyzed_genomes_post_QC_n/Available_genomes_N,NA_real_),
 Weight_N_div_n=ifelse(Analyzed_genomes_post_QC_n>0,
  Available_genomes_N/Analyzed_genomes_post_QC_n,NA_real_),
 .order=match(Year_of_isolation_label,as.character(years)),
 .source_order=match(Source,sources),
 .status_order=match(QR_FQR_Marker_status,statuses)) %>%
 dplyr::arrange(.order,.source_order,.status_order) %>%
 dplyr::select(-.order,-.source_order,-.status_order) %>%
 dplyr::rename("Year of isolation"=Year_of_isolation_label,
 "QR/FQR Marker status"=QR_FQR_Marker_status,
 "Available genomes (N)"=Available_genomes_N,
 "Selected genomes (pre-QC)"=Selected_genomes_pre_QC,
 "Analyzed genomes (post-QC) (n)"=Analyzed_genomes_post_QC_n,
 "% of available selected"=Pct_of_available_selected,
 "% of available analyzed"=Pct_of_available_analyzed,
 "Weight (N ÷ n)"=Weight_N_div_n)
stopifnot(sum(tab[["Available genomes (N)"]])==33831L,
 sum(tab[["Selected genomes (pre-QC)"]])==3434L,
 sum(tab[["Analyzed genomes (post-QC) (n)"]])==3306L)
out <- ensure_dir(sen_path("results/part_a/master_33831/sampling_qc_table"))
readr::write_csv(tab,file.path(out,"year_source_marker_selection_QC_weights.csv"),na="")
if (!requireNamespace("writexl",quietly=TRUE)) stop("Install r-writexl for XLSX export")
writexl::write_xlsx(list(Year_source_marker=tab),file.path(out,"year_source_marker_selection_QC_weights.xlsx"))
cat("[PASS] Sampling and QC table generated: ",nrow(tab)," strata, 33,831 master / 3,434 selected / 3,306 analyzed.\n",sep="")
cat("[INFO] Undefined percentages (N=0) and weights (n=0) exported blank; missing years retained separately.\n")
cat("[OUTPUT] ",out,"\n",sep="")
