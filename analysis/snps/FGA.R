library(tidyverse)

# =========================
# CONFIRMED PATHS
# =========================
base_dir <- "/home/samuelajulo/SENBio/Final/iqtree_final/FGA"

eggnog_file <- file.path(base_dir, "FGA_SplitClade1AB_Annotations.tabular")
master_file <- file.path(base_dir, "Folder_Master_Defining_Genes.csv")

# output dirs
main_outdir <- file.path(base_dir, "COG_plots")
high90_outdir <- file.path(main_outdir, "ge90_only")

dir.create(main_outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(high90_outdir, showWarnings = FALSE, recursive = TRUE)

# =========================
# CLADE COLORS
# =========================
clade_fill_colors <- c(
  "Clade 1A"   = "#5A0000",
  "Clade 1B"   = "#FF0000",
  "Clade 2"    = "#0000FF",
  "Clade 3"    = "#008000",
  "Unassigned" = "#FFFFFF"
)

clade_line_colors <- c(
  "Clade 1A"   = "#5A0000",
  "Clade 1B"   = "#FF0000",
  "Clade 2"    = "#0000FF",
  "Clade 3"    = "#008000",
  "Unassigned" = "#000000"
)

clade_order <- c("Clade 1A", "Clade 1B", "Clade 2", "Clade 3", "Unassigned")

# =========================
# ROBUST EGGNOG READER
# =========================
read_eggnog_tabular <- function(path) {
  lines <- readLines(path, warn = FALSE)

  header_idx <- which(grepl("^#query\\b|^#query_name\\b|^#\\s*query\\b", lines))[1]
  if (is.na(header_idx)) stop("Could not find EggNOG header line starting with #query")

  header_line <- sub("^#", "", lines[header_idx])
  header <- strsplit(header_line, "\t")[[1]]

  data_lines <- lines[(header_idx + 1):length(lines)]
  data_lines <- data_lines[!grepl("^#", data_lines)]
  data_lines <- data_lines[nzchar(data_lines)]

  tf <- tempfile(fileext = ".tsv")
  writeLines(c(paste(header, collapse = "\t"), data_lines), tf)

  df <- read.delim(
    tf,
    sep = "\t",
    header = TRUE,
    quote = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  unlink(tf)
  df
}

# =========================
# LOAD DATA
# =========================
eggnog <- read_eggnog_tabular(eggnog_file)
master <- read.csv(master_file, stringsAsFactors = FALSE, check.names = FALSE)

cat("\nEggNOG columns:\n")
print(colnames(eggnog))

cat("\nMaster columns:\n")
print(colnames(master))

# =========================
# DETECT COLUMNS
# =========================
query_col <- grep("^query$|^#?query$|query", colnames(eggnog), ignore.case = TRUE, value = TRUE)[1]
cog_col   <- grep("COG.*category|COG_category|cog_category", colnames(eggnog), ignore.case = TRUE, value = TRUE)[1]

if (is.na(query_col) || is.na(cog_col)) {
  stop("Could not detect EggNOG query or COG category columns")
}

eggnog2 <- eggnog %>%
  rename(
    Gene_Key = all_of(query_col),
    COG = all_of(cog_col)
  ) %>%
  mutate(
    Gene_Key = as.character(Gene_Key),
    COG = as.character(COG),
    Gene_Key = sub("^gene_query=", "", Gene_Key)
  )

# =========================
# MASTER KEY HANDLING
# =========================
if (!"Gene_Key" %in% colnames(master)) {
  if ("GENE" %in% colnames(master)) {
    master <- master %>% rename(Gene_Key = GENE)
  } else {
    stop("Master file must contain Gene_Key or GENE")
  }
}

if (!"Clade" %in% colnames(master)) {
  stop("Master file must contain Clade")
}

if (!"Coverage_Bins" %in% colnames(master)) {
  stop("Master file must contain Coverage_Bins")
}

if (!"Protein_Sequence" %in% colnames(master)) {
  stop("Master file must contain Protein_Sequence")
}

master2 <- master %>%
  mutate(
    Gene_Key = as.character(Gene_Key),
    Clade = as.character(Clade),
    Coverage_Bins = as.character(Coverage_Bins),
    Protein_Sequence = as.character(Protein_Sequence)
  )

# =========================
# MERGE
# =========================
df <- master2 %>%
  left_join(eggnog2 %>% select(Gene_Key, COG), by = "Gene_Key")

# =========================
# HANDLE NO PROTEIN / NO COG
# =========================
df$COG[is.na(df$COG) | df$COG == "" | df$COG == "-"] <- "NO_PROTEIN"

df <- df %>%
  mutate(
    Protein_Available = ifelse(
      is.na(Protein_Sequence) | Protein_Sequence == "" | Protein_Sequence == "N/A",
      FALSE, TRUE
    )
  )

# =========================
# SAFE SPLIT MULTI-COG
# =========================
split_cog <- function(x) {
  if (is.na(x) || x == "" || x == "-" || x == "NO_PROTEIN") {
    return("NO_PROTEIN")
  }
  strsplit(x, "")[[1]]
}

df_long <- df %>%
  rowwise() %>%
  mutate(COG = list(split_cog(COG))) %>%
  unnest(COG) %>%
  ungroup()

# =========================
# LABELS
# =========================
cog_labels <- c(
  "A" = "RNA processing and modification",
  "B" = "Chromatin structure and dynamics",
  "C" = "Energy production and conversion",
  "D" = "Cell cycle control, cell division, chromosome partitioning",
  "E" = "Amino acid transport and metabolism",
  "F" = "Nucleotide transport and metabolism",
  "G" = "Carbohydrate transport and metabolism",
  "H" = "Coenzyme transport and metabolism",
  "I" = "Lipid transport and metabolism",
  "J" = "Translation, ribosomal structure and biogenesis",
  "K" = "Transcription",
  "L" = "Replication, recombination and repair",
  "M" = "Cell wall/membrane/envelope biogenesis",
  "N" = "Cell motility",
  "O" = "Posttranslational modification, protein turnover, chaperones",
  "P" = "Inorganic ion transport and metabolism",
  "Q" = "Secondary metabolites biosynthesis, transport and catabolism",
  "R" = "General function prediction only",
  "S" = "Function unknown",
  "T" = "Signal transduction mechanisms",
  "U" = "Intracellular trafficking, secretion, and vesicular transport",
  "V" = "Defense mechanisms",
  "W" = "Extracellular structures",
  "Y" = "Nuclear structure",
  "Z" = "Cytoskeleton",
  "NO_PROTEIN" = "No protein / unresolved"
)

cog_order <- c(
  "A","B","C","D","E","F","G","H","I","J","K","L","M",
  "N","O","P","Q","R","S","T","U","V","W","Y","Z","NO_PROTEIN"
)

label_order <- unname(cog_labels[cog_order])

df_long <- df_long %>%
  mutate(
    COG_label = ifelse(COG %in% names(cog_labels), cog_labels[COG], paste("Other:", COG)),
    COG = factor(COG, levels = cog_order),
    COG_label = factor(COG_label, levels = label_order)
  )

# =========================
# HELPERS
# =========================
make_counts <- function(dat) {
  dat %>%
    group_by(Clade, COG, COG_label) %>%
    summarise(Count = n(), .groups = "drop") %>%
    arrange(COG, Clade)
}

make_plot <- function(plot_df) {
  # keep only clades actually present in this iteration
  present_clades <- plot_df %>%
    filter(Count > 0) %>%
    pull(Clade) %>%
    unique()

  present_clades <- clade_order[clade_order %in% present_clades]

  plot_df2 <- plot_df %>%
    filter(Clade %in% present_clades) %>%
    mutate(Clade = factor(Clade, levels = present_clades))

  ggplot(plot_df2, aes(x = COG_label, y = Count, fill = Clade, color = Clade)) +
    geom_bar(
      stat = "identity",
      position = position_dodge(width = 0.8),
      linewidth = 0.6
    ) +
    coord_flip() +
    scale_fill_manual(
      values = clade_fill_colors[present_clades],
      breaks = present_clades,
      drop = TRUE
    ) +
    scale_color_manual(
      values = clade_line_colors[present_clades],
      breaks = present_clades,
      drop = TRUE
    ) +
    guides(
      fill = guide_legend(override.aes = list(
        color = unname(clade_line_colors[present_clades]),
        fill = unname(clade_fill_colors[present_clades])
      )),
      color = "none"
    ) +
    theme_minimal(base_size = 20) +
    labs(
      x = "COG Category",
      y = "Count",
      fill = "Clade"
    ) +
    theme(
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 14),
      legend.position = "right",
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14)
    )
}

write_title_tsv <- function(outdir, stem, title_text) {
  title_df <- tibble(
    File_Stem = stem,
    Plot_Title = title_text
  )
  write.table(
    title_df,
    file = file.path(outdir, paste0(stem, "_TITLE.tsv")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

save_plot_set <- function(plot_df, outdir, stem, title_text) {
  p <- make_plot(plot_df)

  out_png <- file.path(outdir, paste0(stem, ".png"))
  out_tif <- file.path(outdir, paste0(stem, ".tif"))
  out_pdf <- file.path(outdir, paste0(stem, ".pdf"))
  out_csv <- file.path(outdir, paste0(stem, ".csv"))
  out_tsv <- file.path(outdir, paste0(stem, ".tsv"))
  out_title <- file.path(outdir, paste0(stem, "_TITLE.tsv"))

  ggsave(out_png, p, width = 18, height = 12, dpi = 600)
  ggsave(out_tif, p, width = 18, height = 12, dpi = 600, compression = "lzw")
  ggsave(out_pdf, p, width = 18, height = 12, dpi = 600)

  write.csv(plot_df, out_csv, row.names = FALSE)
  write.table(plot_df, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  write_title_tsv(outdir, stem, title_text)

  cat("\nSaved:\n")
  cat("  PNG:   ", out_png, "\n")
  cat("  TIF:   ", out_tif, "\n")
  cat("  PDF:   ", out_pdf, "\n")
  cat("  CSV:   ", out_csv, "\n")
  cat("  TSV:   ", out_tsv, "\n")
  cat("  Title: ", out_title, "\n")
}

make_summary <- function(dat, level_label) {
  gene_df <- dat %>%
    distinct(Clade, Gene_Key, Protein_Available, COG, .keep_all = FALSE)

  gene_df <- gene_df %>%
    mutate(
      COG_clean = ifelse(is.na(COG) | COG == "" | COG == "-" | COG == "NO_PROTEIN", "NO_PROTEIN", COG),
      COG_Letter_Count = ifelse(COG_clean == "NO_PROTEIN", 0L, nchar(COG_clean))
    )

  clade_summary <- gene_df %>%
    group_by(Clade) %>%
    summarise(
      Analysis_Level = level_label,
      Number_of_Genes = n_distinct(Gene_Key),
      Number_of_Proteins_Available = sum(Protein_Available, na.rm = TRUE),
      Number_of_COG_Categories = n_distinct(unlist(strsplit(paste(unique(COG_clean[COG_clean != "NO_PROTEIN"]), collapse = ""), ""))),
      .groups = "drop"
    )

  multi_cog <- gene_df %>%
    filter(COG_Letter_Count >= 2) %>%
    mutate(
      Multi_COG_Summary = paste0(Gene_Key, " (", COG_Letter_Count, " COG)")
    ) %>%
    group_by(Clade, COG_Letter_Count) %>%
    summarise(
      N_Genes = n(),
      Genes = paste(sort(Multi_COG_Summary), collapse = "; "),
      .groups = "drop"
    ) %>%
    mutate(
      Multiple_COG_Group = paste0(N_Genes, " proteins (", COG_Letter_Count, " COG)")
    ) %>%
    arrange(Clade, COG_Letter_Count)

  multi_cog_collapsed <- multi_cog %>%
    group_by(Clade) %>%
    summarise(
      Multiple_COG_Categories = ifelse(
        n() == 0,
        "None",
        paste(paste0(Multiple_COG_Group, ": ", Genes), collapse = " | ")
      ),
      .groups = "drop"
    )

  summary_df <- clade_summary %>%
    left_join(multi_cog_collapsed, by = "Clade") %>%
    mutate(
      Multiple_COG_Categories = ifelse(is.na(Multiple_COG_Categories), "None", Multiple_COG_Categories)
    )

  list(
    summary = summary_df,
    multi_cog_detail = multi_cog
  )
}

save_summary_set <- function(dat, outdir, stem, level_label) {
  s <- make_summary(dat, level_label)

  summary_csv <- file.path(outdir, paste0(stem, "_Summary.csv"))
  summary_tsv <- file.path(outdir, paste0(stem, "_Summary.tsv"))
  multi_csv <- file.path(outdir, paste0(stem, "_Summary_MultiCOG_Detail.csv"))
  multi_tsv <- file.path(outdir, paste0(stem, "_Summary_MultiCOG_Detail.tsv"))

  write.csv(s$summary, summary_csv, row.names = FALSE)
  write.table(s$summary, summary_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

  write.csv(s$multi_cog_detail, multi_csv, row.names = FALSE)
  write.table(s$multi_cog_detail, multi_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

  cat("\nSummary saved:\n")
  cat("  CSV:   ", summary_csv, "\n")
  cat("  TSV:   ", summary_tsv, "\n")
  cat("  Multi: ", multi_csv, "\n")
  cat("  Multi: ", multi_tsv, "\n")
}

# =========================
# MAIN ANALYSIS SETS
# =========================
plot_df_with <- make_counts(df_long)

plot_df_without <- df_long %>%
  filter(Clade != "Unassigned") %>%
  make_counts()

df_long_ge90 <- df_long %>%
  filter(grepl(">=90%", Coverage_Bins, fixed = TRUE))

plot_df_ge90_with <- make_counts(df_long_ge90)

plot_df_ge90_without <- df_long_ge90 %>%
  filter(Clade != "Unassigned") %>%
  make_counts()

# =========================
# SANITY CHECKS
# =========================
cat("\nUnique COG letters detected in full data:\n")
print(sort(unique(as.character(plot_df_with$COG))))

cat("\nMatched rows with non-missing EggNOG COG:\n")
print(sum(df$COG != "NO_PROTEIN"))

cat("\nGenes in >=90% subset:\n")
print(n_distinct(df_long_ge90$Gene_Key))

# =========================
# SAVE MAIN OUTPUTS
# =========================
save_plot_set(
  plot_df_with,
  main_outdir,
  "COG_SideBySide_Flipped_With_Unassigned",
  "COG categories by clade (with Unassigned)"
)
save_summary_set(
  df,
  main_outdir,
  "COG_SideBySide_Flipped_With_Unassigned",
  "With_Unassigned"
)

save_plot_set(
  plot_df_without,
  main_outdir,
  "COG_SideBySide_Flipped_Without_Unassigned",
  "COG categories by clade (without Unassigned)"
)
save_summary_set(
  df %>% filter(Clade != "Unassigned"),
  main_outdir,
  "COG_SideBySide_Flipped_Without_Unassigned",
  "Without_Unassigned"
)

# =========================
# SAVE >=90% SUBFOLDER OUTPUTS
# =========================
save_plot_set(
  plot_df_ge90_with,
  high90_outdir,
  "COG_SideBySide_Flipped_With_Unassigned_ge90",
  "COG categories by clade (>=90% defining SNP genes; with Unassigned)"
)
save_summary_set(
  df %>% filter(grepl(">=90%", Coverage_Bins, fixed = TRUE)),
  high90_outdir,
  "COG_SideBySide_Flipped_With_Unassigned_ge90",
  "With_Unassigned_ge90"
)

save_plot_set(
  plot_df_ge90_without,
  high90_outdir,
  "COG_SideBySide_Flipped_Without_Unassigned_ge90",
  "COG categories by clade (>=90% defining SNP genes; without Unassigned)"
)
save_summary_set(
  df %>% filter(grepl(">=90%", Coverage_Bins, fixed = TRUE), Clade != "Unassigned"),
  high90_outdir,
  "COG_SideBySide_Flipped_Without_Unassigned_ge90",
  "Without_Unassigned_ge90"
)

cat("\n✅ DONE — plots, titles, counts, and summaries generated\n")
