#!/usr/bin/env Rscript
# Original SEN methods SS1 table: isolation year 2018-2025 x QR/FQR marker genotype.
# Exactly one, mutually exclusive genotype category per genome; this avoids
# double-counting genomes carrying multiple QR/FQR markers.
source(file.path(Sys.getenv("SEN_ROOT", unset = getwd()), "analysis/part_a/lib/common.R"))
root <- Sys.getenv("SEN_ROOT", unset = getwd())
master_path <- sen_path("results/part_a/master_33831/00_parsed_raw/SEN_33831_PARSED_51COL.csv")
if (!file.exists(master_path)) stop("Run scripts/run_master_part_a2.sh to generate parsed master first.")
master <- readr::read_csv(master_path,
  col_types = readr::cols(.default = readr::col_character()), show_col_types = FALSE)
assert_columns(master, c("Run","Strain","Isolate","Collection date","AMR genotypes core"))
stopifnot(nrow(master) == 33831L, ncol(master) == 51L,
          dplyr::n_distinct(master$Run) == 33831L)

key <- readr::read_csv(sen_path("data/part_a/selection/Selected_updated_Isolate.csv"),
  col_types=readr::cols(.default=readr::col_character()),show_col_types=FALSE)
assert_columns(key,c("Strain","Isolate"))
norm <- function(z) stringr::str_squish(as.character(z))
key <- key %>% dplyr::transmute(Strain=norm(Strain),Isolate=norm(Isolate)) %>% dplyr::distinct()
stopifnot(nrow(key)==3434L)
mkey <- master %>% dplyr::transmute(Run,Strain=norm(Strain),Isolate=norm(Isolate))
sel <- dplyr::semi_join(mkey,key,by=c("Strain","Isolate"))
stopifnot(nrow(sel)==3434L,dplyr::n_distinct(sel$Run)==3434L)
clade <- readr::read_tsv(sen_path("metadata/final_clade_metadata.tsv"),show_col_types=FALSE)
assert_columns(clade,c("Sample","Lineage"))
analyzed <- unique(as.character(clade$Sample[
  !is.na(clade$Sample) & clade$Sample!="Reference" &
  !is.na(clade$Lineage) & nzchar(clade$Lineage)]))
stopifnot(length(analyzed)==3306L,all(analyzed %in% sel$Run))

major <- c("gyrA_D87Y","gyrA_D87N","gyrA_S83Y","gyrA_D87G","qnrB19")
levels <- c(major,"Multiple major QR/FQR markers",
            "Other QR/FQR marker(s)","QS/FQS (no marker)")
get_group <- function(markers) {
  if (length(markers)==0L) return("QS/FQS (no marker)")
  m <- unique(tolower(markers))
  hits <- major[vapply(tolower(major),function(t) t %in% m,logical(1))]
  if (length(hits)>1L) return("Multiple major QR/FQR markers")
  if (length(hits)==1L) return(hits)
  "Other QR/FQR marker(s)"
}
markers <- extract_qr_core_markers(master[["AMR genotypes core"]])
master$QR_marker_group <- vapply(markers,get_group,character(1))
master$QR_status <- ifelse(lengths(markers)==0L,
 "QS/FQS (no marker)","QR/FQR marker-positive")
stopifnot(identical(master$QR_status,
 ifelse(derive_qr_genotype_core(master[["AMR genotypes core"]])==1L,
 "QR/FQR marker-positive","QS/FQS (no marker)")))
master$Year <- as.integer(derive_collection_year(master[["Collection date"]]))
master$Selected_pre_QC <- master$Run %in% sel$Run
master$Analyzed_post_QC <- master$Run %in% analyzed
stopifnot(sum(master$Selected_pre_QC)==3434L,sum(master$Analyzed_post_QC)==3306L,
 all(!master$Analyzed_post_QC | master$Selected_pre_QC))
records <- master %>% dplyr::filter(!is.na(Year),dplyr::between(Year,2018L,2025L))
tab <- records %>% dplyr::group_by(Year,QR_marker_group) %>%
 dplyr::summarise(N_available=dplyr::n(),
  N_selected=sum(Selected_pre_QC),N_analyzed=sum(Analyzed_post_QC),.groups="drop") %>%
 tidyr::complete(Year=2018:2025,QR_marker_group=levels,
  fill=list(N_available=0L,N_selected=0L,N_analyzed=0L)) %>%
 dplyr::mutate(.sort=match(QR_marker_group,levels)) %>%
 dplyr::arrange(Year,.sort) %>% dplyr::select(-.sort) %>%
 dplyr::rename("Year of isolation"=Year,
 "QR/FQR marker genotype"=QR_marker_group,
 "Total available genomes (N)"=N_available,
 "Total genomes selected (pre-QC)"=N_selected,
 "Total genomes analyzed (post-QC)"=N_analyzed)
stopifnot(sum(tab[["Total available genomes (N)"]])==nrow(records),
 sum(tab[["Total genomes selected (pre-QC)"]])==sum(records$Selected_pre_QC),
 sum(tab[["Total genomes analyzed (post-QC)"]])==sum(records$Analyzed_post_QC))
# Reconcile the derived marker groupings against the already validated
# year/source/QR status accounting; this prevents silent genotype-mapping drift.
annual_path <- sen_path("results/part_a/master_33831/sampling_qc_table/year_source_marker_selection_QC_weights.csv")
if (file.exists(annual_path)) {
 annual <- readr::read_csv(annual_path,show_col_types=FALSE)
 check <- records %>% dplyr::group_by(Year,QR_status) %>%
  dplyr::summarise(N=dplyr::n(),selected=sum(Selected_pre_QC),
    analyzed=sum(Analyzed_post_QC),.groups="drop")
 recorded <- annual %>% dplyr::filter(`Year of isolation` %in% as.character(2018:2025)) %>%
  dplyr::mutate(Year=as.integer(`Year of isolation`)) %>%
  dplyr::group_by(Year,QR_status=`QR/FQR Marker status`) %>%
  dplyr::summarise(N=sum(`Available genomes (N)`),
    selected=sum(`Selected genomes (pre-QC)`),
    analyzed=sum(`Analyzed genomes (post-QC) (n)`),.groups="drop")
 if (!isTRUE(all.equal(dplyr::arrange(check,Year,QR_status),
                 dplyr::arrange(recorded,Year,QR_status),check.attributes=FALSE)))
   stop("Genotype table does not reconcile against existing yearly QR-status table.")
}
# Full genotype signatures (not grouped by major marker). These capture
# co-occurring determinants without pretending each genome has one mutation.
master$QR_genotype_signature <- vapply(markers,function(x) {
  if (!length(x)) "QS/FQS (no marker)" else
    paste(sort(unique(pretty_token(x))),collapse=" + ")
},character(1))
full <- master %>% dplyr::filter(!is.na(Year),dplyr::between(Year,2018L,2025L)) %>%
  dplyr::group_by(Year,QR_genotype_signature) %>%
  dplyr::summarise(N_available=dplyr::n(),
    N_selected=sum(Selected_pre_QC),N_analyzed=sum(Analyzed_post_QC),
    .groups="drop") %>%
  dplyr::arrange(Year,dplyr::desc(N_available),QR_genotype_signature) %>%
  dplyr::rename("Year of isolation"=Year,
    "Exact QR/FQR marker genotype(s)"=QR_genotype_signature,
    "Total available genomes (N)"=N_available,
    "Total genomes selected (pre-QC)"=N_selected,
    "Total genomes analyzed (post-QC)"=N_analyzed)
stopifnot(sum(full[["Total available genomes (N)"]])==nrow(records),
 sum(full[["Total genomes selected (pre-QC)"]])==sum(records$Selected_pre_QC),
 sum(full[["Total genomes analyzed (post-QC)"]])==sum(records$Analyzed_post_QC))
out <- ensure_dir(sen_path("results/part_a/master_33831/sampling_qc_table"))
readr::write_csv(tab,file.path(out,"Table_SS1_2018_2025_year_marker_genotype.csv"))
readr::write_csv(full,file.path(out,"Table_SS1_2018_2025_exact_genotype_signatures.csv"))
if (!requireNamespace("writexl",quietly=TRUE)) stop("r-writexl required for XLSX export.")
rules <- data.frame(
 Note=c("Study population","Time window","One genome per row",
        "Major QR/FQR markers","Multi-marker classification",
        "Other QR/FQR category","QS/FQS category","Count meanings",
        "Outbreak control"),
 Definition=c(
  "N is from the original 33,831-genome master; selected n=3,434; analyzed n=3,306.",
  "2018-2025 inclusive, consistent with the 2018-2020 and 2021-2025 selection blocks.",
  "The categories are mutually exclusive. Counts may be added across genotype rows within a year.",
  paste(major,collapse=", "),
  "If a genome carries two or more of the listed major markers, assign Multiple major QR/FQR markers.",
  "Any other curated QR/FQR determinant, with none of the five major markers, is grouped as Other QR/FQR marker(s).",
  "No curated QR/FQR marker in the AMR genotypes core field. Marker-negative is not an AST-confirmed susceptibility phenotype.",
  "Available=master genomes; selected=3,434-key membership; analyzed=3,306 final clade metadata.",
  "The previously verified outbreak restriction is a feature of the archived selected cohort, not reapplied here."
 ),check.names=FALSE)
writexl::write_xlsx(list("Table_SS1"=tab,"Exact_genotype_signatures"=full,"Definitions"=rules),
 file.path(out,"Table_SS1_2018_2025_year_marker_genotype.xlsx"))
cat("[PASS] Table SS1 ",nrow(tab)," rows / 2018-2025 / genotypes + QS/FQS; ",
 sum(tab[["Total available genomes (N)"]])," available, ",
 sum(tab[["Total genomes selected (pre-QC)"]])," selected, ",
 sum(tab[["Total genomes analyzed (post-QC)"]])," analyzed.\\n",sep="")
cat("[OUTPUT] ",file.path(out,"Table_SS1_2018_2025_year_marker_genotype.xlsx"),"\\n",sep="")
