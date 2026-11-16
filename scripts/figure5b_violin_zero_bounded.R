#!/usr/bin/env Rscript
# Regenerate SEN Figure 5B from the original r_snps.R method.
# Source: https://github.com/ajulojays/SENBio-Final/blob/78abf1e0efb3ed2eb91d7de72fc928a3b90e64f5/r_snps.R
# Change to original plotting: geom_violin(trim = TRUE) instead of trim = FALSE.
# Default inputs: data/phylogeny/SSLAB_FINAL.treefile,
#   metadata/final_clade_metadata.tsv and metadata/SEN_Genomes.csv.
# Optional positional arguments override treefile, clade metadata, output directory,
#   and master CSV, in that order. Run from anywhere using Rscript.
# NOTE: This plots patristic tree distances * 67,736 alignment sites,
#       NOT directly counted pairwise SNP differences.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 4) {
  stop("Usage: Rscript scripts/figure5b_violin_zero_bounded.R [TREEFILE [CLADE_METADATA.tsv [OUTPUT_DIR [SEN_Genomes.csv]]]]", call. = FALSE)
}
script_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_args) < 1) stop("Run this script with Rscript.", call. = FALSE)
script_file <- normalizePath(sub("^--file=", "", script_args[[1]]), mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)

for (pkg in c("ape", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Required R package '%s' is missing. Install it first.", pkg), call. = FALSE)
  }
}

treefile <- if (length(args) >= 1) path.expand(args[[1]]) else file.path(repo_root, "data", "phylogeny", "SSLAB_FINAL.treefile")
metadata_file <- if (length(args) >= 2) path.expand(args[[2]]) else file.path(repo_root, "metadata", "final_clade_metadata.tsv")
outdir <- if (length(args) >= 3) path.expand(args[[3]]) else file.path(repo_root, "results", "figure5b_zero_bounded")
master_default <- file.path(repo_root, "metadata", "SEN_Genomes.csv")
master_file <- if (length(args) >= 4) path.expand(args[[4]]) else if (file.exists(master_default)) master_default else NULL
for (f in c(treefile, metadata_file, master_file)) {
  if (!is.null(f) && !file.exists(f)) stop("File not found: ", f, call. = FALSE)
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

alignment_length <- 67736L  # same as original r_snps.R
clade_levels <- c("Clade 1A", "Clade 1B", "Clade 2", "Clade 3")
clade_colors <- c(
  "Clade 1A" = "#6B0000",
  "Clade 1B" = "#FF0000",
  "Clade 2"  = "#0000FF",
  "Clade 3"  = "#008000"
)

cat("[INFO] Reading tree: ", treefile, "\n", sep = "")
tree <- ape::read.tree(treefile)
if (anyDuplicated(tree$tip.label)) stop("Duplicate tree tip names found.", call. = FALSE)
cat("[INFO] Tree tips: ", length(tree$tip.label), "\n", sep = "")

metadata <- read.delim(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("Sample", "Lineage") %in% names(metadata))) {
  stop("Clade metadata must have columns 'Sample' and 'Lineage'.", call. = FALSE)
}
metadata$Sample <- trimws(as.character(metadata$Sample))
metadata$Sample <- sub(".*?(SRR[0-9]+).*", "\\1", metadata$Sample, perl = TRUE)
metadata$Lineage <- trimws(as.character(metadata$Lineage))
metadata <- metadata[metadata$Lineage %in% clade_levels, c("Sample", "Lineage")]
if (anyDuplicated(metadata$Sample)) stop("Duplicate sample IDs in clade metadata.", call. = FALSE)

if (!is.null(master_file)) {
  master <- read.csv(master_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"Run" %in% names(master)) stop("Master CSV must contain a Run column.", call. = FALSE)
  master_ids <- sub(".*?(SRR[0-9]+).*", "\\1", as.character(master$Run), perl = TRUE)
  metadata <- metadata[metadata$Sample %in% master_ids, ]
  cat("[INFO] Restricted to samples present in master CSV, matching original r_snps.R.\n")
}

metadata <- metadata[metadata$Sample %in% tree$tip.label, ]
cat("[INFO] Matched focal-clade genomes: ", nrow(metadata), "\n", sep = "")
if (!nrow(metadata)) stop("No clade IDs matched the tree tip names.", call. = FALSE)

cat("[INFO] Calculating tree patristic distances × 67,736...\n")
patristic_as_sites <- ape::cophenetic.phylo(tree) * alignment_length

# The original script used all ordered pairs (A,B) and (B,A), excluding
# diagonal self-pairs. Build those within each clade without a huge long table
# of between-clade pairs; duplication reproduces the original KDE weighting.
parts <- list()
audit <- list()
for (group in clade_levels) {
  ids <- unique(metadata$Sample[metadata$Lineage == group])
  n <- length(ids)
  if (n < 2) stop("Fewer than 2 matched genomes in ", group, call. = FALSE)
  submat <- patristic_as_sites[ids, ids, drop = FALSE]
  unique_pairs <- as.numeric(submat[lower.tri(submat)])
  if (any(!is.finite(unique_pairs)) || any(unique_pairs < 0)) {
    stop("Invalid/nonfinite or negative tree-derived distances in ", group, call. = FALSE)
  }
  parts[[group]] <- data.frame(
    Lineage = group,
    SNPs = rep(unique_pairs, times = 2L),
    stringsAsFactors = FALSE
  )
  audit[[group]] <- data.frame(
    Lineage = group,
    N_genomes = n,
    N_unique_pairs = length(unique_pairs),
    Min = min(unique_pairs),
    Median = median(unique_pairs),
    Max = max(unique_pairs)
  )
  cat(sprintf("[INFO] %-8s genomes=%4d unique_pairs=%8d min=%8.3f median=%8.3f max=%8.3f\n",
              group, n, length(unique_pairs), min(unique_pairs),
              median(unique_pairs), max(unique_pairs)))
  rm(submat, unique_pairs)
}
rm(patristic_as_sites)

within_data <- do.call(rbind, parts)
within_data$Lineage <- factor(within_data$Lineage, levels = clade_levels)
summary <- do.call(rbind, audit)
write.csv(summary, file.path(outdir, "Figure5B_distance_summary.csv"), row.names = FALSE)

p <- ggplot2::ggplot(within_data, ggplot2::aes(x = Lineage, y = SNPs, fill = Lineage)) +
  ggplot2::geom_violin(alpha = 0.8, trim = TRUE, color = "black", linewidth = 0.3) +
  ggplot2::theme_classic() +
  ggplot2::scale_fill_manual(values = clade_colors, breaks = clade_levels, drop = FALSE) +
  ggplot2::coord_cartesian(ylim = c(0, NA)) +
  ggplot2::labs(x = "Lineage", y = "Pairwise SNP distance", fill = "Lineage") +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = ggplot2::element_text(size = 20),
    axis.title.x = ggplot2::element_text(size = 28),
    axis.title.y = ggplot2::element_text(size = 28),
    legend.title = ggplot2::element_text(size = 24),
    legend.text = ggplot2::element_text(size = 24)
  )

stem <- file.path(outdir, "Figure5B_violin_zero_bounded")
ggplot2::ggsave(paste0(stem, ".png"), p, width = 10, height = 6, dpi = 300)
ggplot2::ggsave(paste0(stem, ".pdf"), p, width = 10, height = 6)
ggplot2::ggsave(paste0(stem, ".tif"), p, width = 10, height = 6, dpi = 300,
                device = "tiff", compression = "lzw")
cat("[PASS] Saved PNG, PDF, TIFF and distance summary to: ", outdir, "\n", sep = "")