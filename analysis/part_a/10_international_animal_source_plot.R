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
