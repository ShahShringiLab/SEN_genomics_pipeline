#!/usr/bin/env Rscript
# Publication-ready breakdown of the 71 International SEN isolates by
# ANIMAL/DECOMPOSED source, separate from the original QR/FQR stacked plots.
# Derived entirely from original full-master Table S1 (chunk 17); no resampling.
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
})

root <- Sys.getenv("SEN_ROOT", unset=getwd())
source_csv <- file.path(root, "results/part_a/full_master_original",
                        "SEN_TableS1_final_source_decomposition.csv")
out_dir <- file.path(root, "results/part_a/full_master_original",
                     "international_animal_source")
if (!file.exists(source_csv)) {
  stop("Original full-master Table S1 is missing: ", source_csv,
       "\nFirst complete original full-master analysis (chunks 3-22).")
}

source_table <- readr::read_csv(source_csv, show_col_types=FALSE)
required <- c("Source category", "Decomposed source information",
              "Number of genomes (n)")
if (!all(required %in% names(source_table))) {
  stop("Unexpected original Table S1 columns: ",
       paste(setdiff(required, names(source_table)), collapse=", "))
}

international <- source_table %>%
  dplyr::transmute(
    Source = trimws(as.character(.data[["Source category"]])),
    Animal_source = trimws(as.character(.data[["Decomposed source information"]])),
    n = as.integer(.data[["Number of genomes (n)"]])
  ) %>%
  dplyr::filter(Source == "International") %>%
  dplyr::arrange(match(Animal_source, c("Others", "Chicken")))

# Strict checks ensure this represents exactly the original 71 isolates
# instead of an accidental selected (3,434) or master-wide (33,831) table.
stopifnot(
  nrow(international) == 2L,
  !anyDuplicated(international$Animal_source),
  setequal(international$Animal_source, c("Others", "Chicken")),
  !anyNA(international$n),
  all(international$n >= 0L),
  sum(international$n) == 71L,
  international$n[international$Animal_source == "Others"] == 67L,
  international$n[international$Animal_source == "Chicken"] == 4L
)

international <- international %>%
  dplyr::mutate(
    Percentage = 100 * n / sum(n),
    Label = sprintf("%s (%.1f%%)", scales::comma(n), Percentage),
    Animal_source = factor(Animal_source, levels=c("Chicken", "Others"))
  )

dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
base <- file.path(out_dir, "SEN_International_71_decomposed_animal_sources")
readr::write_csv(international, paste0(base, "_data.csv"))

p <- ggplot2::ggplot(international,
                     ggplot2::aes(x=Animal_source, y=n, fill=Animal_source)) +
  ggplot2::geom_col(width=0.57, show.legend=FALSE) +
  ggplot2::geom_text(ggplot2::aes(label=Label),
                     hjust=-0.12, fontface="bold", size=4.5) +
  ggplot2::coord_flip(clip="off") +
  ggplot2::scale_fill_manual(values=c("Chicken"="#3575A6",
                                      "Others"="#8595A6")) +
  ggplot2::scale_y_continuous(
    limits=c(0,80),
    breaks=seq(0,70,10),
    labels=scales::label_comma(),
    expand=ggplot2::expansion(mult=c(0,0))
  ) +
  ggplot2::labs(
    title="International SEN genomes (n = 71)",
    subtitle="Animal/decomposed isolation sources",
    x=NULL,
    y="Number of genomes"
  ) +
  ggplot2::theme_classic(base_size=13) +
  ggplot2::theme(
    plot.title=ggplot2::element_text(hjust=0.5,face="bold",size=16),
    plot.subtitle=ggplot2::element_text(hjust=0.5,size=12),
    axis.text.y=ggplot2::element_text(face="bold",size=12),
    axis.text.x=ggplot2::element_text(size=11),
    plot.margin=ggplot2::margin(12,34,12,12)
  )

for (ext in c("png","tif","pdf")) {
  args <- list(filename=paste0(base,".",ext), plot=p,
               width=8.6, height=4.5, units="in", bg="white")
  if (ext != "pdf") args$dpi <- 600
  if (ext == "tif") args$compression <- "lzw"
  do.call(ggplot2::ggsave,args)
}

print(international)
cat("[PASS] International decomposed animal-source figure: ",
    "Others=67 (94.4%); Chicken=4 (5.6%); N=71\n",sep="")
cat("[OUTPUT] ",base,".png/.tif/.pdf; source-data CSV\n",sep="")


# ============================================================
# DETAILED decomposition of the 67 "Others", retaining the 4 chicken samples
# and the exact original isolation-source text for isolate-level traceability.
# No generated category is claimed to be an original surveillance taxonomy.
# ============================================================

selected_enriched_path <- file.path(
  root,"results/part_a/reconciled_metadata",
  "SEN_3434_selected_enriched_QC_metadata.csv"
)
selected_qc_path <- file.path(
  root,"results/part_a/cohort",
  "Finalupdated_SEN_Genomes1_QC_merged.csv"
)
selected_archived_path <- file.path(
  root,"data/part_a/selection","Finalupdated_SEN_Genomes1.csv"
)
master_enriched_path <- file.path(
  root,"results/part_a/reconciled_metadata",
  "SEN_33831_full_enriched_metadata.csv"
)
master_raw_path <- file.path(
  root,"results/part_a/master_33831/00_parsed_raw",
  "SEN_33831_PARSED_51COL.csv"
)

# Prefer all curated selected + QC + Rmd fields, then the 78-column QC
# merge, then the original archived 69-column selected metadata.
sel_input <- c(selected_enriched_path, selected_qc_path, selected_archived_path)
sel_input <- sel_input[file.exists(sel_input)]
if (!length(sel_input)) stop("No selected metadata found for International audit.")
sel_input <- sel_input[1]
read_all_text <- function(p) readr::read_csv(
  p, col_types=readr::cols(.default=readr::col_character()),
  na=character(), show_col_types=FALSE, progress=FALSE
)
sel_meta <- read_all_text(sel_input)
needed <- c("Run","Isolate","Source","Source_group",
            "Isolation source","IFSAC category","Food origin")
stopifnot(all(needed %in% names(sel_meta)))
international_meta <- sel_meta %>%
  dplyr::filter(Source=="International") %>%
  dplyr::arrange(Run)
stopifnot(nrow(international_meta)==71L,
          dplyr::n_distinct(international_meta$Run)==71L,
          dplyr::n_distinct(international_meta$Isolate)==71L,
          all(international_meta$Source_group %in%
                c("Chicken-Non-US","Other-Non-US")),
          sum(international_meta$Source_group=="Other-Non-US")==67L,
          sum(international_meta$Source_group=="Chicken-Non-US")==4L)

# A deliberately simple rule-based grouping of the ACTUAL literal source
# strings. Original source strings are never edited. Keyword priority resolves
# processed commodities such as "cheese roll cake" into baked/sweet foods,
# while "imitation crab meat" is NOT classified as crustacean.
source_group_detail <- function(isolation, original_group) {
  x <- tolower(trimws(as.character(isolation)))
  if (identical(original_group,"Chicken-Non-US"))
    return("Chicken-containing pet treats")
  has <- function(pattern) grepl(pattern,x,perl=TRUE)
  if (has("rabbit")) return("Rabbit meat/products")
  if (has("\\b(shrimp|crab|lobster)\\b") &&
      !has("imitation\\s+crab")) return("Crustaceans")
  if (has("\\b(squid|cuttlefish|snail|mussel|clam)\\b"))
    return("Mollusks")
  if (has("imitation\\s+crab")) return("Processed/imitated seafood")
  if (has("\\b(fish|tilapia|catfish|pompano|pomfret|pangasius|scad|esomus)\\b"))
    return("Finfish")
  if (has("\\b(guinea\\s+pig|quail)\\b"))
    return("Other terrestrial animals")
  if (has("\\b(cake|biscuit|cookie|cookies|marshmallow|snack\\s+bar)\\b"))
    return("Bakery/confectionery")
  if (has("\\b(cheese|morbier)\\b")) return("Dairy/cheese")
  if (has("\\bdog\\s+chew\\b")) return("Pet chew (meat unspecified)")
  if (has("\\b(pepper|chili|fenugreek|thyme|spice|spiced|almond|dandelion|cilantro|melon\\s+seed|coriander|kava)\\b"))
    return("Plants/herbs/spices/nuts/seeds")
  return("Unclassified — review literal source")
}

source_details <- international_meta %>%
  dplyr::transmute(
    Run, Isolate,
    Original_source_category=Source,
    Original_decomposed_source=dplyr::if_else(
      Source_group=="Chicken-Non-US","Chicken","Others"),
    Original_source_group=Source_group,
    Isolation_source_original=.data[["Isolation source"]],
    IFSAC_category_original=.data[["IFSAC category"]],
    Food_origin_original=.data[["Food origin"]]
  )
source_details$Detailed_source_group <- vapply(
  seq_len(nrow(source_details)),
  function(i) source_group_detail(
    source_details$Isolation_source_original[i],
    source_details$Original_source_group[i]
  ), character(1)
)
stopifnot(nrow(source_details)==71L,
          sum(source_details$Original_decomposed_source=="Others")==67L,
          sum(source_details$Original_decomposed_source=="Chicken")==4L)
unknown <- source_details %>%
  dplyr::filter(Detailed_source_group=="Unclassified — review literal source")
if(nrow(unknown)) {
  warning(nrow(unknown),
          " unknown source descriptions retained explicitly. Review the audit sheet.")
}

# No synthetic, anonymous, or inferred "animal" descriptions. We group exact
# isolation descriptions into clearly labeled commodity/source classes.
group_counts <- source_details %>%
  dplyr::count(Detailed_source_group,name="Genomes") %>%
  dplyr::mutate(Pct_of_71=round(100*Genomes/71,2)) %>%
  dplyr::arrange(dplyr::desc(Genomes),Detailed_source_group)
other_details <- source_details %>%
  dplyr::filter(Original_decomposed_source=="Others")
other_counts <- other_details %>%
  dplyr::count(Detailed_source_group,name="Genomes") %>%
  dplyr::mutate(Pct_of_67=round(100*Genomes/67,2),
                Pct_of_71=round(100*Genomes/71,2)) %>%
  dplyr::arrange(dplyr::desc(Genomes),Detailed_source_group)
exact_counts <- other_details %>%
  dplyr::mutate(Exact_isolation_description=dplyr::if_else(
    trimws(Isolation_source_original)=="",
    "[Missing Isolation source]",trimws(Isolation_source_original))) %>%
  dplyr::count(Exact_isolation_description,name="Genomes") %>%
  dplyr::mutate(Pct_of_67=round(100*Genomes/67,2)) %>%
  dplyr::arrange(dplyr::desc(Genomes),Exact_isolation_description)
ifsac_counts <- source_details %>%
  dplyr::count(IFSAC_category_original,name="Genomes") %>%
  dplyr::mutate(Pct_of_71=round(100*Genomes/71,2)) %>%
  dplyr::arrange(dplyr::desc(Genomes))
country_counts <- source_details %>%
  dplyr::count(Food_origin_original,name="Genomes") %>%
  dplyr::mutate(Pct_of_71=round(100*Genomes/71,2)) %>%
  dplyr::arrange(dplyr::desc(Genomes))
original_counts <- source_details %>%
  dplyr::count(Original_decomposed_source,name="Genomes") %>%
  dplyr::mutate(Pct_of_71=round(100*Genomes/71,2))
stopifnot(
  sum(group_counts$Genomes)==71L,
  sum(other_counts$Genomes)==67L,
  sum(exact_counts$Genomes)==67L,
  sum(ifsac_counts$Genomes)==71L,
  sum(country_counts$Genomes)==71L,
  sum(original_counts$Genomes)==71L
)

# Attach the explicit decomposition alongside EVERY original selected metadata
# column and, where available, all QC/phylogeny and derived Rmd fields.
selected_sheet <- dplyr::left_join(
  source_details, international_meta, by=c("Run","Isolate")
)
stopifnot(nrow(selected_sheet)==71L)
master_path <- if(file.exists(master_enriched_path))
  master_enriched_path else master_raw_path
if (!file.exists(master_path))
  stop("Full master original metadata unavailable: ",master_path)
master <- read_all_text(master_path)
stopifnot(all(c("Run","Isolate") %in% names(master)),
          dplyr::n_distinct(master$Run)==33831L)
master_sheet <- master %>%
  dplyr::semi_join(source_details,by=c("Run","Isolate")) %>%
  dplyr::left_join(source_details,by=c("Run","Isolate")) %>%
  dplyr::arrange(Run)
stopifnot(nrow(master_sheet)==71L,
          setequal(master_sheet$Run,source_details$Run))

summary_sheet <- tibble::tibble(
  Metric=c("International total", "Original Others", "Original Chicken",
    "Distinct Others isolation descriptions", "Distinct decomposed source groups",
    "Distinct exporting countries (Food origin)", "Selected metadata file",
    "Full master metadata file", "Interpretation"),
  Value=c("71","67","4",
    as.character(nrow(exact_counts)),
    as.character(nrow(group_counts)),
    as.character(nrow(country_counts)),
    sel_input,master_path,
    "Detailed groups use transparent isolation-source keyword rules; exact original terms preserved.")
)
definitions <- tibble::tibble(
  Field=c("Source_group","Isolation source","Detailed_source_group",
          "IFSAC category","Food origin","Original Others"),
  Interpretation=c(
    "Historical source-group field: Chicken-Non-US or Other-Non-US",
    "Verbatim NCBI Pathogen Detection isolation-source description",
    "Rule-based commodity category, NOT an original surveillance classification; see isolated-source values",
    "Original IFSAC commodity taxonomy, preserved as provided",
    "Reported food origin, not necessarily isolation/collection country",
    "All 67 non-chicken International samples in the original published decomposition")
)
workbook <- list(
  Summary=summary_sheet,
  Groups_all71=group_counts,
  Groups_Others67=other_counts,
  Exact_sources_Others67=exact_counts,
  Original_67_4=original_counts,
  Countries_71=country_counts,
  IFSAC_71=ifsac_counts,
  Isolate_classification=source_details,
  All71_selected_QC=selected_sheet,
  All71_master_metadata=master_sheet,
  Definitions=definitions
)
if(!requireNamespace("writexl",quietly=TRUE))
  stop("Install R package 'writexl' in sen_stats environment to create Excel workbook.")
xlsx_path <- file.path(out_dir,"SEN_International_71_ALL_sources_and_metadata.xlsx")
writexl::write_xlsx(workbook,xlsx_path)

plot_breakdown <- function(df, category_col, count_col, denom, title,
                           subtitle, outfile, height=NULL, small_font=11L) {
  plotdf <- df %>%
    dplyr::transmute(
      Category=as.character(.data[[category_col]]),
      N=as.integer(.data[[count_col]])
    ) %>%
    dplyr::arrange(N,Category) %>%
    dplyr::mutate(
      Category=factor(Category,levels=Category),
      Label=sprintf("%d (%.1f%%)",N,100*N/denom)
    )
  if(is.null(height)) height <- max(5,2.5+0.38*nrow(plotdf))
  xmax <- max(plotdf$N)+max(1,0.35*max(plotdf$N))
  graph <- ggplot2::ggplot(plotdf,ggplot2::aes(x=Category,y=N)) +
    ggplot2::geom_col(width=0.74,fill="#3979A8") +
    ggplot2::geom_text(ggplot2::aes(label=Label),
                       hjust=-0.08,size=3.4) +
    ggplot2::coord_flip(clip="off") +
    ggplot2::scale_y_continuous(
      limits=c(0,xmax),breaks=scales::breaks_pretty(n=5),
      expand=ggplot2::expansion(mult=c(0,0))
    ) +
    ggplot2::labs(title=title,subtitle=subtitle,
                  x=NULL,y="Number of genomes") +
    ggplot2::theme_classic(base_size=12) +
    ggplot2::theme(
      plot.title=ggplot2::element_text(hjust=0.5,face="bold",size=16),
      plot.subtitle=ggplot2::element_text(hjust=0.5,size=11),
      axis.text.y=ggplot2::element_text(size=small_font),
      plot.margin=ggplot2::margin(14,40,14,14)
    )
  for(ext in c("png","tif","pdf")) {
    args <- list(filename=file.path(out_dir,paste0(outfile,".",ext)),
                 plot=graph,width=12,height=height,units="in",bg="white")
    if(ext!="pdf") args$dpi <- 400
    if(ext=="tif") args$compression <- "lzw"
    do.call(ggplot2::ggsave,args)
  }
}

plot_breakdown(
  group_counts,"Detailed_source_group","Genomes",71,
  "International SEN genomes (n = 71): detailed source breakdown",
  "Source groups derived from original isolation descriptions; includes Chicken (4)",
  "SEN_International_71_detailed_source_groups"
)
plot_breakdown(
  other_counts,"Detailed_source_group","Genomes",67,
  "International SEN: what is inside Others (n = 67)?",
  "Original non-chicken International samples; transparent commodity groups",
  "SEN_International_Others_67_group_breakdown"
)
plot_breakdown(
  exact_counts,"Exact_isolation_description","Genomes",67,
  "International SEN Others (n = 67): exact isolation sources",
  "Every original isolation-source description, without merging categories",
  "SEN_International_Others_67_exact_isolation_sources",
  height=max(12,2+0.36*nrow(exact_counts)),small_font=8L
)

readr::write_csv(group_counts,file.path(out_dir,
  "SEN_International_71_detailed_source_groups.csv"))
readr::write_csv(other_counts,file.path(out_dir,
  "SEN_International_Others_67_groups.csv"))
readr::write_csv(exact_counts,file.path(out_dir,
  "SEN_International_Others_67_exact_isolation_sources.csv"))
readr::write_csv(source_details,file.path(out_dir,
  "SEN_International_71_isolate_source_audit.csv"))

cat("[PASS] International Others (67) decomposed by literal source terms: ",
    nrow(other_counts)," commodity groups, ",
    nrow(exact_counts)," distinct exact isolation descriptions\n",sep="")
cat("[PASS] International workbook: 71 full master rows, 71 selected",
    " rows with available QC, explicit source traceability\n")
cat("[OUTPUT] ",xlsx_path,"\n",sep="")
cat("[OUTPUT] Three detailed breakdown charts each in PNG/TIFF/PDF: ",out_dir,"\n",sep="")
