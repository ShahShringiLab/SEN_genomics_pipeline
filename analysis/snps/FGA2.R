library(tidyverse)
library(grid)

# ============================================================
# SEN FGA / COG PIPELINE
#
# Outputs:
#   1. One master package for each of the 8 iterations
#   2. Gene-level master
#   3. SNP-level master
#   4. Gene × COG master
#   5. SNP × COG master
#   6. SNP coverage audit
#   7. QC table
#   8. COG plots
#
# IMPORTANT:
# Multi-SNP genes are expanded to one row per SNP before
# coverage-bin and impact-class analysis.
#
# Example:
# SEN_RS14875:
#   SNP 1 -> 60-<80%
#   SNP 2 -> >=90%
#
# This prevents the previous 417-vs-419 discrepancy.
# ============================================================


# ============================================================
# PATHS
# ============================================================

sen_root <- Sys.getenv("SEN_ROOT", unset = normalizePath(getwd(), mustWork = FALSE))

fga_dir <- Sys.getenv(
  "SEN_FGA_DIR",
  unset = file.path(sen_root, "iqtree_final", "FGA")
)

snps_root <- Sys.getenv(
  "SEN_SNP_INTEGRATED_ROOT",
  unset = file.path(sen_root, "iqtree_final", "Snps_clade_integrated")
)

master_root <- file.path(
  fga_dir,
  "MASTER_BY_ITERATION"
)

dir.create(
  master_root,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# EGGNOG FILES
# ============================================================

anno_map <- list(

  Combined_Clade1 =
    file.path(
      fga_dir,
      "FGA_Combined_Annotations.tabular"
    ),

  Split_Clade1A1B =
    file.path(
      fga_dir,
      "FGA_SplitClade1AB_Annotations.tabular"
    )
)


# ============================================================
# ORDERS
# ============================================================

impact_order <- c(
  "High",
  "Moderate",
  "Low"
)

coverage_order <- c(
  "40-<60%",
  "60-<80%",
  "80-<90%",
  ">=90%"
)

cog_order <- c(
  "A","B","C","D","E","F","G","H","I","J","K","L","M",
  "N","O","P","Q","R","S","T","U","V","W","Y","Z",
  "NO_PROTEIN"
)


# ============================================================
# COG LABELS
# ============================================================

cog_labels <- c(

  A = "RNA processing and modification",

  B = "Chromatin structure and dynamics",

  C = "Energy production and conversion",

  D = "Cell cycle control, cell division, chromosome partitioning",

  E = "Amino acid transport and metabolism",

  F = "Nucleotide transport and metabolism",

  G = "Carbohydrate transport and metabolism",

  H = "Coenzyme transport and metabolism",

  I = "Lipid transport and metabolism",

  J = "Translation, ribosomal structure and biogenesis",

  K = "Transcription",

  L = "Replication, recombination and repair",

  M = "Cell wall/membrane/envelope biogenesis",

  N = "Cell motility",

  O = "Posttranslational modification, protein turnover, chaperones",

  P = "Inorganic ion transport and metabolism",

  Q = "Secondary metabolites biosynthesis, transport and catabolism",

  R = "General function prediction only",

  S = "Function unknown",

  T = "Signal transduction mechanisms",

  U = "Intracellular trafficking, secretion, and vesicular transport",

  V = "Defense mechanisms",

  W = "Extracellular structures",

  Y = "Nuclear structure",

  Z = "Cytoskeleton",

  NO_PROTEIN =
    "No protein / unresolved"
)

label_order <- unname(
  cog_labels[cog_order]
)


# ============================================================
# IMPACT FUNCTION
# ============================================================

classify_impact <- function(effect_string) {

  if (
    is.na(effect_string) ||
    trimws(effect_string) == ""
  ) {
    return("Low")
  }

  effects <- unique(
    trimws(
      unlist(
        strsplit(
          effect_string,
          ";",
          fixed = TRUE
        )
      )
    )
  )

  if (
    any(
      effects %in% c(
        "frameshift_variant",
        "stop_gained",
        "stop_lost",
        "start_lost",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "exon_loss_variant",
        "feature_ablation",
        "transcript_ablation"
      )
    )
  ) {

    return("High")
  }

  if (
    any(
      effects %in% c(
        "missense_variant",
        "protein_altering_variant",
        "inframe_insertion",
        "inframe_deletion"
      )
    )
  ) {

    return("Moderate")
  }

  return("Low")
}


# ============================================================
# READ EGGNOG
# ============================================================

read_eggnog <- function(path) {

  if (!file.exists(path)) {

    stop(
      "Missing EggNOG file: ",
      path
    )
  }

  lines <- readLines(
    path,
    warn = FALSE
  )

  header_idx <- grep(
    "^#query",
    lines
  )[1]

  if (is.na(header_idx)) {

    stop(
      "Could not find #query header in: ",
      path
    )
  }

  header <- strsplit(
    sub(
      "^#",
      "",
      lines[header_idx]
    ),
    "\t"
  )[[1]]

  data_lines <- lines[
    !grepl(
      "^#",
      lines
    )
  ]

  tmp <- tempfile(
    fileext = ".tsv"
  )

  writeLines(
    c(
      paste(
        header,
        collapse = "\t"
      ),
      data_lines
    ),
    tmp
  )

  out <- read.delim(
    tmp,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  unlink(tmp)

  return(out)
}


# ============================================================
# PALETTE
# ============================================================

make_distinct_palette <- function(labels) {

  labels <- sort(
    unique(
      as.character(labels)
    )
  )

  labels <- labels[
    !is.na(labels)
  ]

  base_pal <- c(

    "#1F77B4",
    "#D62728",
    "#2CA02C",
    "#9467BD",
    "#FF7F0E",
    "#17BECF",
    "#8C564B",
    "#E377C2",
    "#BCBD22",
    "#7F7F7F",
    "#AEC7E8",
    "#FF9896",
    "#98DF8A",
    "#C5B0D5",
    "#FFBB78",
    "#9EDAE5",
    "#C49C94",
    "#F7B6D2",
    "#DBDB8D",
    "#393B79",
    "#637939",
    "#8C6D31",
    "#843C39",
    "#7B4173",
    "#3182BD",
    "#E6550D",
    "#31A354",
    "#756BB1"
  )

  if (
    length(labels) <=
    length(base_pal)
  ) {

    pal <- base_pal[
      seq_along(labels)
    ]

  } else {

    pal <- c(

      base_pal,

      grDevices::hcl.colors(
        length(labels) -
          length(base_pal),
        palette = "Dark 3"
      )
    )
  }

  names(pal) <- labels

  if (
    "No protein / unresolved" %in%
    names(pal)
  ) {

    pal[
      "No protein / unresolved"
    ] <- "#FFFFFF"
  }

  return(pal)
}


# ============================================================
# PARSE ITERATION
# ============================================================

parse_iteration_meta <- function(
    mf,
    root
) {

  rel <- sub(
    paste0(
      "^",
      root,
      "/"
    ),
    "",
    mf
  )

  parts <- strsplit(
    rel,
    "/",
    fixed = TRUE
  )[[1]]

  if (
    length(parts) < 3
  ) {

    stop(
      "Unexpected master path: ",
      mf
    )
  }

  tibble(

    Scheme =
      parts[1],

    SNP_Mode =
      parts[2],

    Assignment =
      parts[3]
  )
}


# ============================================================
# SPLIT SEMICOLON VALUES
# ============================================================

split_semicolon <- function(x) {

  if (
    is.na(x) ||
    trimws(
      as.character(x)
    ) == ""
  ) {

    return(
      character(0)
    )
  }

  trimws(
    strsplit(
      as.character(x),
      ";",
      fixed = TRUE
    )[[1]]
  )
}


# ============================================================
# NORMALIZE NUMBER OF VALUES
# ============================================================

normalize_vector <- function(
    values,
    expected_n,
    field_name,
    gene_key
) {

  if (
    length(values) ==
    expected_n
  ) {

    return(values)
  }

  if (
    length(values) == 1 &&
    expected_n > 1
  ) {

    return(
      rep(
        values,
        expected_n
      )
    )
  }

  if (
    length(values) == 0
  ) {

    return(
      rep(
        NA_character_,
        expected_n
      )
    )
  }

  stop(

    "Cannot safely expand ",
    field_name,
    " for gene ",
    gene_key,
    ". Expected ",
    expected_n,
    " values but found ",
    length(values)
  )
}


# ============================================================
# EXPAND GENE ROWS TO SNP ROWS
# ============================================================

expand_gene_rows_to_snps <- function(df) {

  out <- vector(
    "list",
    nrow(df)
  )

  for (
    i in seq_len(
      nrow(df)
    )
  ) {

    row <- df[
      i,
      ,
      drop = FALSE
    ]

    ids <- split_semicolon(
      row$Variant_IDs
    )

    effects <- split_semicolon(
      row$Variant_Effects
    )

    bins <- split_semicolon(
      row$Coverage_Bins
    )

    expected_n <- suppressWarnings(
      as.integer(
        row$Defining_SNP_Count_In_Gene
      )
    )

    if (
      is.na(expected_n) ||
      expected_n < 1
    ) {

      expected_n <- max(

        length(ids),

        length(effects),

        length(bins),

        1
      )
    }


    if (
      length(ids) > 1 &&
      length(ids) !=
      expected_n
    ) {

      stop(

        "Variant_ID mismatch for ",
        row$Gene_Key,
        ": Defining_SNP_Count_In_Gene = ",
        expected_n,
        " but Variant_IDs contains ",
        length(ids),
        " values."
      )
    }


    ids <- normalize_vector(

      ids,
      expected_n,
      "Variant_IDs",
      row$Gene_Key
    )


    effects <- normalize_vector(

      effects,
      expected_n,
      "Variant_Effects",
      row$Gene_Key
    )


    bins <- normalize_vector(

      bins,
      expected_n,
      "Coverage_Bins",
      row$Gene_Key
    )


    expanded <- row[
      rep(
        1,
        expected_n
      ),
      ,
      drop = FALSE
    ]


    expanded$SNP_Index_In_Gene <-
      seq_len(
        expected_n
      )


    expanded$Variant_ID <-
      ids


    expanded$Variant_Effect <-
      effects


    expanded$Coverage_Bin <-
      bins


    expanded$Impact_Class <-
      vapply(

        effects,

        classify_impact,

        character(1)
      )


    expanded$Defining_SNP_Count_In_Gene_Original <-
      expected_n


    out[[i]] <- expanded
  }

  bind_rows(out)
}


# ============================================================
# EXPAND MULTIPLE COG CATEGORIES
# ============================================================

expand_cog <- function(df) {

  df %>%

    mutate(

      COG_original =
        as.character(
          COG_original
        ),

      COG =
        purrr::map(

          COG_original,

          function(x) {

            if (
              is.na(x) ||
              x == "" ||
              x == "-" ||
              x == "NO_PROTEIN"
            ) {

              return(
                "NO_PROTEIN"
              )
            }

            unique(
              strsplit(
                x,
                ""
              )[[1]]
            )
          }
        )
    ) %>%

    unnest(COG) %>%

    mutate(

      COG =
        as.character(COG),

      COG_label =
        ifelse(

          COG %in%
            names(
              cog_labels
            ),

          unname(
            cog_labels[
              COG
            ]
          ),

          paste(
            "Other:",
            COG
          )
        )
    )
}


# ============================================================
# FIND THE 8 MASTER FILES
# ============================================================

master_files <- list.files(

  snps_root,

  pattern =
    "Folder_Master_Defining_Genes.csv$",

  recursive = TRUE,

  full.names = TRUE
)


master_files <- master_files[

  grepl(

    "Combined_Clade1|Split_Clade1A1B",

    master_files
  )
]


master_files <- master_files[

  !grepl(

    "/old/",

    master_files,

    ignore.case = TRUE
  )
]


if (
  length(master_files) == 0
) {

  stop(
    "No master files found."
  )
}


cat(
  "\nFound",
  length(master_files),
  "iteration master files\n"
)


# ============================================================
# VERIFY EXACTLY ONE FILE PER ITERATION
# ============================================================

meta_check <- map_dfr(

  master_files,

  function(mf) {

    parse_iteration_meta(
      mf,
      snps_root
    ) %>%

      mutate(
        Master_File = mf
      )
  }
)


duplicate_iterations <- meta_check %>%

  count(
    Scheme,
    SNP_Mode,
    Assignment
  ) %>%

  filter(
    n > 1
  )


if (
  nrow(
    duplicate_iterations
  ) > 0
) {

  print(
    duplicate_iterations
  )

  stop(
    "Duplicate master files found."
  )
}


# ============================================================
# STORAGE
# ============================================================

all_counts_for_final <- list()

all_snps_for_final <- list()


# ============================================================
# MAIN LOOP
# ============================================================

for (
  mf in master_files
) {

  cat(
    "\n============================================================\n"
  )

  cat(
    "Processing:\n",
    mf,
    "\n"
  )

  cat(
    "============================================================\n"
  )


  # ----------------------------------------------------------
  # META
  # ----------------------------------------------------------

  meta <- parse_iteration_meta(
    mf,
    snps_root
  )


  scheme <- meta$Scheme[[1]]

  snp_mode <- meta$SNP_Mode[[1]]

  assignment <- meta$Assignment[[1]]


  anno_file <- anno_map[[scheme]]


  if (
    is.null(
      anno_file
    )
  ) {

    stop(
      "No EggNOG annotation file mapped for scheme: ",
      scheme
    )
  }


  # ----------------------------------------------------------
  # READ MASTER
  # ----------------------------------------------------------

  master <- read.csv(

    mf,

    stringsAsFactors = FALSE,

    check.names = FALSE
  )


  # ----------------------------------------------------------
  # STANDARDIZE Gene_Key
  # ----------------------------------------------------------

  if (
    !"Gene_Key" %in%
    names(master)
  ) {

    if (
      "GENE" %in%
      names(master)
    ) {

      master <- master %>%

        rename(
          Gene_Key =
            GENE
        )

    } else {

      stop(
        "Gene_Key/GENE column absent in ",
        mf
      )
    }
  }


  # ----------------------------------------------------------
  # STANDARDIZE COVERAGE COLUMN
  # ----------------------------------------------------------

  if (
    !"Coverage_Bins" %in%
    names(master)
  ) {

    if (
      "Coverage_Bin" %in%
      names(master)
    ) {

      master <- master %>%

        rename(
          Coverage_Bins =
            Coverage_Bin
        )

    } else {

      stop(
        "Coverage_Bins absent in ",
        mf
      )
    }
  }


  # ----------------------------------------------------------
  # STANDARDIZE EFFECT COLUMN
  # ----------------------------------------------------------

  if (
    !"Variant_Effects" %in%
    names(master)
  ) {

    if (
      "EFFECT" %in%
      names(master)
    ) {

      master <- master %>%

        rename(
          Variant_Effects =
            EFFECT
        )

    } else {

      stop(
        "Variant_Effects absent in ",
        mf
      )
    }
  }


  # ----------------------------------------------------------
  # REQUIRED FIELDS
  # ----------------------------------------------------------

  required_cols <- c(

    "Clade",

    "Gene_Key",

    "Variant_IDs",

    "Variant_Effects",

    "Coverage_Bins",

    "Defining_SNP_Count_In_Gene"
  )


  missing_cols <- setdiff(

    required_cols,

    names(master)
  )


  if (
    length(
      missing_cols
    ) > 0
  ) {

    stop(

      "Missing required columns: ",

      paste(
        missing_cols,
        collapse = ", "
      )
    )
  }


  # ==========================================================
  # READ EGGNOG
  # ==========================================================

  egg <- read_eggnog(
    anno_file
  )


  query_col <- grep(

    "^query$|query",

    colnames(egg),

    ignore.case = TRUE,

    value = TRUE
  )[1]


  cog_col <- grep(

    "COG",

    colnames(egg),

    ignore.case = TRUE,

    value = TRUE
  )[1]


  if (
    is.na(query_col) ||
    is.na(cog_col)
  ) {

    stop(
      "Could not identify EggNOG query/COG columns."
    )
  }


  # ==========================================================
  # ONE EGGNOG ROW PER GENE
  # ==========================================================

  egg2 <- egg %>%

    rename(

      Gene_Key =
        all_of(
          query_col
        ),

      COG =
        all_of(
          cog_col
        )
    ) %>%

    mutate(

      Gene_Key =
        sub(

          "^gene_query=",

          "",

          as.character(
            Gene_Key
          )
        ),

      COG =
        as.character(COG),

      COG =
        ifelse(

          is.na(COG) |
            COG == "" |
            COG == "-",

          NA_character_,

          COG
        )
    ) %>%

    group_by(
      Gene_Key
    ) %>%

    summarise(

      COG = {

        vals <- na.omit(
          COG
        )

        if (
          length(vals) == 0
        ) {

          NA_character_

        } else {

          chars <- unique(

            unlist(

              strsplit(

                paste(
                  vals,
                  collapse = ""
                ),

                ""
              )
            )
          )

          chars <- chars[
            chars != ""
          ]

          paste(
            chars,
            collapse = ""
          )
        }
      },

      .groups =
        "drop"
    )


  # ==========================================================
  # GENE-LEVEL MASTER
  # ==========================================================

  df_gene <- master %>%

    mutate(

      Gene_Row_ID =
        row_number(),

      Gene_Key =
        as.character(
          Gene_Key
        ),

      Clade =
        as.character(
          Clade
        )
    ) %>%

    left_join(

      egg2,

      by =
        "Gene_Key"
    ) %>%

    mutate(

      COG_original =
        ifelse(

          is.na(COG) |
            COG == "" |
            COG == "-",

          "NO_PROTEIN",

          COG
        ),

      Gene_Impact_Class =
        vapply(

          Variant_Effects,

          classify_impact,

          character(1)
        ),

      Analysis_Scheme =
        scheme,

      SNP_Mode =
        snp_mode,

      Assignment =
        assignment,

      Source_Master_File =
        mf
    ) %>%

    select(
      -COG
    )


  # ==========================================================
  # SNP-LEVEL MASTER
  # ==========================================================

  df_snp <- expand_gene_rows_to_snps(
    df_gene
  )


  # ==========================================================
  # VALIDATE COVERAGE BINS
  # ==========================================================

  invalid_bins <- df_snp %>%

    filter(

      is.na(
        Coverage_Bin
      ) |

      !Coverage_Bin %in%
        coverage_order
    )


  if (
    nrow(
      invalid_bins
    ) > 0
  ) {

    cat(
      "\nUnexpected coverage-bin values:\n"
    )

    print(

      as.data.frame(

        invalid_bins %>%

          select(

            Clade,

            Gene_Key,

            Variant_ID,

            Coverage_Bin
          )
      ),

      row.names = FALSE
    )

    stop(
      "Coverage-bin expansion failed."
    )
  }


  # ==========================================================
  # SNP COUNT CHECK
  # ==========================================================

  expected_snps <- sum(

    df_gene$Defining_SNP_Count_In_Gene,

    na.rm = TRUE
  )


  observed_snps <- nrow(
    df_snp
  )


  if (
    expected_snps !=
    observed_snps
  ) {

    stop(

      "SNP expansion mismatch. Expected ",

      expected_snps,

      " but recovered ",

      observed_snps
    )
  }


  # ==========================================================
  # EXPAND COG
  # ==========================================================

  df_gene_cog <- expand_cog(
    df_gene
  )


  df_snp_cog <- expand_cog(
    df_snp
  )


  # ==========================================================
  # CLADE ORDER
  # ==========================================================

  if (
    scheme ==
    "Split_Clade1A1B"
  ) {

    clade_order <- c(

      "Clade 1A",

      "Clade 1B",

      "Clade 2",

      "Clade 3",

      "Unassigned"
    )

  } else {

    clade_order <- c(

      "Clade 1",

      "Clade 2",

      "Clade 3",

      "Unassigned"
    )
  }


  # ==========================================================
  # FACTORS
  # ==========================================================

  df_snp_cog <- df_snp_cog %>%

    mutate(

      Coverage_Bin =
        factor(

          Coverage_Bin,

          levels =
            coverage_order,

          ordered =
            TRUE
        ),

      Impact_Class =
        factor(

          Impact_Class,

          levels =
            impact_order,

          ordered =
            TRUE
        ),

      Clade =
        factor(

          Clade,

          levels =
            clade_order
        ),

      COG =
        factor(

          COG,

          levels =
            cog_order
        ),

      COG_label =
        factor(

          COG_label,

          levels =
            label_order
        )
    )


  # ==========================================================
  # PLOT COUNTS
  # ==========================================================

  counts <- df_snp_cog %>%

    count(

      Coverage_Bin,

      Impact_Class,

      Clade,

      COG_label,

      name =
        "n"
    ) %>%

    filter(
      n > 0
    )


  # ==========================================================
  # ITERATION NAME
  # ==========================================================

  iter_name <- paste(

    scheme,

    snp_mode,

    assignment,

    sep =
      "__"
  )


  # ==========================================================
  # MASTER SUBFOLDER
  # ==========================================================

  master_iter_dir <- file.path(

    master_root,

    iter_name
  )


  dir.create(

    master_iter_dir,

    showWarnings = FALSE,

    recursive = TRUE
  )


  # ==========================================================
  # SAVE GENE MASTER
  # ==========================================================

  write.csv(

    df_gene,

    file.path(

      master_iter_dir,

      "Genes_Annotated.csv"
    ),

    row.names = FALSE
  )


  write.table(

    df_gene,

    file.path(

      master_iter_dir,

      "Genes_Annotated.tsv"
    ),

    sep = "\t",

    row.names = FALSE,

    quote = FALSE
  )


  # ==========================================================
  # SAVE GENE × COG MASTER
  # ==========================================================

  write.csv(

    df_gene_cog,

    file.path(

      master_iter_dir,

      "Genes_Expanded_By_COG.csv"
    ),

    row.names = FALSE
  )


  write.table(

    df_gene_cog,

    file.path(

      master_iter_dir,

      "Genes_Expanded_By_COG.tsv"
    ),

    sep = "\t",

    row.names = FALSE,

    quote = FALSE
  )


  # ==========================================================
  # SAVE SNP MASTER
  # ==========================================================

  write.csv(

    df_snp,

    file.path(

      master_iter_dir,

      "SNPs_Annotated.csv"
    ),

    row.names = FALSE
  )


  write.table(

    df_snp,

    file.path(

      master_iter_dir,

      "SNPs_Annotated.tsv"
    ),

    sep = "\t",

    row.names = FALSE,

    quote = FALSE
  )


  # ==========================================================
  # SAVE SNP × COG MASTER
  # ==========================================================

  write.csv(

    df_snp_cog,

    file.path(

      master_iter_dir,

      "SNPs_Expanded_By_COG.csv"
    ),

    row.names = FALSE
  )


  write.table(

    df_snp_cog,

    file.path(

      master_iter_dir,

      "SNPs_Expanded_By_COG.tsv"
    ),

    sep = "\t",

    row.names = FALSE,

    quote = FALSE
  )


  # ==========================================================
  # SNP COVERAGE AUDIT
  # ==========================================================

  coverage_audit <- df_snp %>%

    count(

      Clade,

      Coverage_Bin,

      name =
        "Defining_SNPs"
    ) %>%

    arrange(

      Clade,

      factor(

        Coverage_Bin,

        levels =
          coverage_order
      )
    )


  write.table(

    coverage_audit,

    file.path(

      master_iter_dir,

      "SNP_COVERAGE_AUDIT.tsv"
    ),

    sep = "\t",

    row.names = FALSE,

    quote = FALSE
  )


  # ==========================================================
  # CLADE FUNCTIONAL SUMMARY
  #
  # Captures ALL clades present in this iteration.
  # One row per:
  #   Clade × Coverage_Bin × Impact_Class × COG
  #
  # Counts include:
  #   - unique SNPs in the stratum
  #   - unique genes in the stratum
  #   - total SNP × COG assignments
  #   - COG-specific SNP/gene counts
  # ==========================================================

  base_functional_counts <- df_snp_cog %>%
    group_by(
      Clade,
      Coverage_Bin,
      Impact_Class
    ) %>%
    summarise(
      Unique_SNPs =
        n_distinct(
          Variant_ID
        ),

      Unique_Genes =
        n_distinct(
          Gene_Key
        ),

      SNP_x_COG_Assignments =
        n(),

      .groups =
        "drop"
    )


  cog_functional_counts <- df_snp_cog %>%
    group_by(
      Clade,
      Coverage_Bin,
      Impact_Class,
      COG,
      COG_label
    ) %>%
    summarise(
      COG_Assignment_n =
        n(),

      Unique_SNPs_in_COG =
        n_distinct(
          Variant_ID
        ),

      Unique_Genes_in_COG =
        n_distinct(
          Gene_Key
        ),

      .groups =
        "drop"
    )


  clade_functional_summary <- cog_functional_counts %>%
    left_join(
      base_functional_counts,
      by = c(
        "Clade",
        "Coverage_Bin",
        "Impact_Class"
      )
    ) %>%
    mutate(
      Clade =
        as.character(
          Clade
        ),

      Coverage_Bin =
        as.character(
          Coverage_Bin
        ),

      Impact_Class =
        as.character(
          Impact_Class
        ),

      COG =
        as.character(
          COG
        ),

      COG_label =
        as.character(
          COG_label
        )
    ) %>%
    select(
      Clade,
      Coverage_Bin,
      Impact_Class,
      Unique_SNPs,
      Unique_Genes,
      SNP_x_COG_Assignments,
      COG,
      COG_label,
      COG_Assignment_n,
      Unique_SNPs_in_COG,
      Unique_Genes_in_COG
    ) %>%
    arrange(
      Clade,

      factor(
        Coverage_Bin,
        levels = coverage_order
      ),

      factor(
        Impact_Class,
        levels = impact_order
      ),

      desc(
        COG_Assignment_n
      ),

      COG_label
    )


  write.table(
    clade_functional_summary,

    file.path(
      master_iter_dir,
      "Clade_Functional_Summary.tsv"
    ),

    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )


  write.csv(
    clade_functional_summary,

    file.path(
      master_iter_dir,
      "Clade_Functional_Summary.csv"
    ),

    row.names = FALSE
  )


  # ==========================================================
  # CLADE OVERVIEW
  #
  # One row per clade in the iteration.
  # ==========================================================

  clade_overview <- df_snp_cog %>%
    mutate(
      Clade =
        as.character(
          Clade
        ),

      Impact_Class =
        as.character(
          Impact_Class
        )
    ) %>%
    group_by(
      Clade
    ) %>%
    summarise(
      Total_Unique_SNPs =
        n_distinct(
          Variant_ID
        ),

      Total_Unique_Genes =
        n_distinct(
          Gene_Key
        ),

      Total_SNP_x_COG_Assignments =
        n(),

      High_Impact_SNPs =
        n_distinct(
          Variant_ID[
            Impact_Class == "High"
          ]
        ),

      Moderate_Impact_SNPs =
        n_distinct(
          Variant_ID[
            Impact_Class == "Moderate"
          ]
        ),

      Low_Impact_SNPs =
        n_distinct(
          Variant_ID[
            Impact_Class == "Low"
          ]
        ),

      .groups =
        "drop"
    )


  write.table(
    clade_overview,

    file.path(
      master_iter_dir,
      "Clade_Overview.tsv"
    ),

    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )


  write.csv(
    clade_overview,

    file.path(
      master_iter_dir,
      "Clade_Overview.csv"
    ),

    row.names = FALSE
  )


  cat(
    "\nFUNCTIONAL SUMMARY SAVED:",
    file.path(
      master_iter_dir,
      "Clade_Functional_Summary.tsv"
    ),
    "\n"
  )



  # ==========================================================
  # QC
  # ==========================================================

  master_qc <- tibble(

    Scheme =
      scheme,

    SNP_Mode =
      snp_mode,

    Assignment =
      assignment,

    Gene_rows =
      nrow(
        df_gene
      ),

    Unique_Gene_Keys =
      n_distinct(
        df_gene$Gene_Key
      ),

    Expected_Defining_SNPs =
      expected_snps,

    Recovered_SNP_rows =
      observed_snps,

    Gene_x_COG_rows =
      nrow(
        df_gene_cog
      ),

    SNP_x_COG_rows =
      nrow(
        df_snp_cog
      ),

    Additional_Gene_COG_rows =
      nrow(
        df_gene_cog
      ) -
      nrow(
        df_gene
      ),

    Additional_SNP_COG_rows =
      nrow(
        df_snp_cog
      ) -
      nrow(
        df_snp
      )
  )


  write.table(

    master_qc,

    file.path(

      master_iter_dir,

      "MASTER_QC.tsv"
    ),

    sep = "\t",

    row.names = FALSE,

    quote = FALSE
  )


  cat(

    "\nGENES:",
    nrow(df_gene),

    "| DEFINING SNPs:",
    observed_snps,

    "| Gene x COG:",
    nrow(df_gene_cog),

    "| SNP x COG:",
    nrow(df_snp_cog),

    "\n"
  )


  # ==========================================================
  # ITERATION PLOT DIRECTORY
  # ==========================================================

  iter_outdir <- file.path(

    fga_dir,

    paste0(

      "ITER_",

      iter_name
    )
  )


  dir.create(

    iter_outdir,

    showWarnings = FALSE,

    recursive = TRUE
  )


  # ==========================================================
  # ITERATION PLOT
  # ==========================================================

  present_cogs <- sort(

    unique(

      as.character(

        counts$COG_label
      )
    )
  )


  present_cogs <- present_cogs[
    !is.na(
      present_cogs
    )
  ]


  cog_palette <- make_distinct_palette(
    present_cogs
  )


  p <- ggplot(

    counts,

    aes(

      x =
        Clade,

      y =
        n,

      fill =
        COG_label
    )
  ) +

    geom_bar(

      stat =
        "identity",

      position =
        "stack",

      width =
        0.7,

      color =
        "black",

      linewidth =
        0.25
    ) +

    coord_flip() +

    facet_grid(

      Impact_Class ~
        Coverage_Bin,

      scales =
        "free_x"
    ) +

    scale_fill_manual(

      values =
        cog_palette,

      name =
        "COG",

      drop =
        TRUE
    ) +

    labs(

      x =
        "Clade",

      y =
        "Defining SNP × COG assignments"
    ) +

    theme_minimal(
      base_size = 20
    ) +

    theme(

      axis.text.x =
        element_text(
          size = 20
        ),

      axis.text.y =
        element_text(
          size = 20
        ),

      axis.title =
        element_text(
          size = 24
        ),

      legend.text =
        element_text(
          size = 20
        ),

      legend.title =
        element_text(
          size = 20
        ),

      legend.key.size =
        unit(
          0.8,
          "cm"
        ),

      strip.text =
        element_text(
          size = 18,
          face = "bold"
        ),

      panel.border =
        element_rect(
          color = "black",
          fill = NA,
          linewidth = 1.2
        ),

      panel.spacing =
        unit(
          0.6,
          "lines"
        )
    )


  # ==========================================================
  # SAVE ITERATION PLOTS
  # ==========================================================

  ggsave(

    file.path(

      iter_outdir,

      "plot.png"
    ),

    p,

    width =
      24,

    height =
      18,

    dpi =
      600
  )


  ggsave(

    file.path(

      iter_outdir,

      "plot.tif"
    ),

    p,

    width =
      24,

    height =
      18,

    dpi =
      600,

    compression =
      "lzw"
  )


  ggsave(

    file.path(

      iter_outdir,

      "plot.pdf"
    ),

    p,

    width =
      24,

    height =
      18
  )


  write.csv(

    counts,

    file.path(

      iter_outdir,

      "counts.csv"
    ),

    row.names =
      FALSE
  )


  write.table(

    counts,

    file.path(

      iter_outdir,

      "counts.tsv"
    ),

    sep =
      "\t",

    row.names =
      FALSE,

    quote =
      FALSE
  )


  write.table(

    tibble(

      File_Stem =
        "plot",

      Plot_Title =
        paste0(

          "COG distribution of defining SNPs by clade, impact class, and coverage bin: ",

          iter_name
        )
    ),

    file.path(

      iter_outdir,

      "plot_TITLE.tsv"
    ),

    sep =
      "\t",

    row.names =
      FALSE,

    quote =
      FALSE
  )


  # ==========================================================
  # FINAL STORAGE
  #
  # No in-memory append is needed here.
  # The final section below rebuilds all final objects directly
  # from the 8 saved MASTER_BY_ITERATION folders.
  # ==========================================================


}

# ============================================================
# IMPORTANT:
# The two list assignments above must use valid R [[ ]]
# indexing. Rebuild them explicitly here to avoid accidental
# formatting corruption.
# ============================================================

# Recreate final storage safely from iteration master folders.
# This avoids relying on any malformed multiline [[ syntax.

all_counts_for_final <- list()
all_snps_for_final <- list()


master_iteration_dirs <- list.dirs(
  master_root,
  recursive = FALSE,
  full.names = TRUE
)


for (
  d in master_iteration_dirs
) {

  qc_file <- file.path(
    d,
    "MASTER_QC.tsv"
  )

  snp_file <- file.path(
    d,
    "SNPs_Annotated.tsv"
  )

  snp_cog_file <- file.path(
    d,
    "SNPs_Expanded_By_COG.tsv"
  )


  if (
    !file.exists(qc_file) ||
    !file.exists(snp_file) ||
    !file.exists(snp_cog_file)
  ) {

    next
  }


  qc <- read.delim(
    qc_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


  snp_tmp <- read.delim(
    snp_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


  snp_cog_tmp <- read.delim(
    snp_cog_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


  counts_tmp <- snp_cog_tmp %>%

    mutate(

      Coverage_Bin =
        factor(
          Coverage_Bin,
          levels =
            coverage_order,
          ordered =
            TRUE
        ),

      Impact_Class =
        factor(
          Impact_Class,
          levels =
            impact_order,
          ordered =
            TRUE
        ),

      COG_label =
        factor(
          COG_label,
          levels =
            label_order
        )
    ) %>%

    count(

      Coverage_Bin,

      Impact_Class,

      Clade,

      COG_label,

      name =
        "n"
    ) %>%

    filter(
      n > 0
    ) %>%

    mutate(

      Scheme =
        qc$Scheme[1],

      SNP_Mode =
        qc$SNP_Mode[1],

      Assignment =
        qc$Assignment[1]
    )


  all_counts_for_final <-
    append(
      all_counts_for_final,
      list(
        counts_tmp
      )
    )


  snp_tmp <- snp_tmp %>%

    mutate(

      Scheme =
        qc$Scheme[1],

      SNP_Mode =
        qc$SNP_Mode[1],

      Assignment =
        qc$Assignment[1]
    )


  all_snps_for_final <-
    append(
      all_snps_for_final,
      list(
        snp_tmp
      )
    )
}


# ============================================================
# FINAL CANONICAL ANALYSIS
#
# Clade 1:
#   Combined_Clade1
#
# Clades 1A, 1B, 2, 3, Unassigned:
#   Split_Clade1A1B
#
# vertical_only
# with_unassigned
# ============================================================

final_counts_raw <- bind_rows(
  all_counts_for_final
)


combined_part <- final_counts_raw %>%

  filter(

    Scheme ==
      "Combined_Clade1",

    SNP_Mode ==
      "vertical_only",

    Assignment ==
      "with_unassigned",

    as.character(
      Clade
    ) ==
      "Clade 1"
  )


split_part <- final_counts_raw %>%

  filter(

    Scheme ==
      "Split_Clade1A1B",

    SNP_Mode ==
      "vertical_only",

    Assignment ==
      "with_unassigned",

    as.character(
      Clade
    ) %in% c(

      "Clade 1A",

      "Clade 1B",

      "Clade 2",

      "Clade 3",

      "Unassigned"
    )
  )


final_counts <- bind_rows(

  combined_part,

  split_part
) %>%

  select(

    Coverage_Bin,

    Impact_Class,

    Clade,

    COG_label,

    n
  ) %>%

  group_by(

    Coverage_Bin,

    Impact_Class,

    Clade,

    COG_label
  ) %>%

  summarise(

    n =
      sum(n),

    .groups =
      "drop"
  )


final_clade_order <- c(

  "Clade 1",

  "Clade 1A",

  "Clade 1B",

  "Clade 2",

  "Clade 3",

  "Unassigned"
)


final_counts <- final_counts %>%

  mutate(

    Coverage_Bin =
      factor(

        Coverage_Bin,

        levels =
          coverage_order,

        ordered =
          TRUE
      ),

    Impact_Class =
      factor(

        Impact_Class,

        levels =
          impact_order,

        ordered =
          TRUE
      ),

    Clade =
      factor(

        as.character(
          Clade
        ),

        levels =
          final_clade_order
      ),

    COG_label =
      factor(

        as.character(
          COG_label
        ),

        levels =
          label_order
      )
  ) %>%

  filter(
    !is.na(
      Clade
    )
  )


# ============================================================
# FINAL SNP MASTER
# ============================================================

final_snp_raw <- bind_rows(
  all_snps_for_final
)


final_snp <- bind_rows(

  final_snp_raw %>%

    filter(

      Scheme ==
        "Combined_Clade1",

      SNP_Mode ==
        "vertical_only",

      Assignment ==
        "with_unassigned",

      Clade ==
        "Clade 1"
    ),


  final_snp_raw %>%

    filter(

      Scheme ==
        "Split_Clade1A1B",

      SNP_Mode ==
        "vertical_only",

      Assignment ==
        "with_unassigned",

      Clade %in% c(

        "Clade 1A",

        "Clade 1B",

        "Clade 2",

        "Clade 3",

        "Unassigned"
      )
    )
)


# ============================================================
# EXPECTED VERTICAL SNP COUNTS
# ============================================================

expected_final <- tribble(

  ~Clade,
  ~Coverage_Bin,
  ~Expected_SNPs,


  "Clade 1",
  "40-<60%",
  33,


  "Clade 1",
  "60-<80%",
  3,


  "Clade 1",
  "80-<90%",
  4,


  "Clade 1",
  ">=90%",
  0,


  "Clade 1A",
  "40-<60%",
  4,


  "Clade 1A",
  "60-<80%",
  2,


  "Clade 1A",
  "80-<90%",
  0,


  "Clade 1A",
  ">=90%",
  0,


  "Clade 1B",
  "40-<60%",
  3,


  "Clade 1B",
  "60-<80%",
  31,


  "Clade 1B",
  "80-<90%",
  0,


  "Clade 1B",
  ">=90%",
  5,


  "Clade 2",
  "40-<60%",
  15,


  "Clade 2",
  "60-<80%",
  65,


  "Clade 2",
  "80-<90%",
  0,


  "Clade 2",
  ">=90%",
  22,


  "Clade 3",
  "40-<60%",
  3,


  "Clade 3",
  "60-<80%",
  5,


  "Clade 3",
  "80-<90%",
  1,


  "Clade 3",
  ">=90%",
  81,


  "Unassigned",
  "40-<60%",
  142,


  "Unassigned",
  "60-<80%",
  0,


  "Unassigned",
  "80-<90%",
  0,


  "Unassigned",
  ">=90%",
  0
)


# ============================================================
# OBSERVED VERTICAL SNP COUNTS
# ============================================================

final_snp_audit <- final_snp %>%

  count(

    Clade,

    Coverage_Bin,

    name =
      "Observed_SNPs"
  )


# ============================================================
# VALIDATION
# ============================================================

final_validation <- expected_final %>%

  left_join(

    final_snp_audit,

    by = c(

      "Clade",

      "Coverage_Bin"
    )
  ) %>%

  mutate(

    Observed_SNPs =
      replace_na(

        Observed_SNPs,

        0L
      ),

    Difference =
      Observed_SNPs -
      Expected_SNPs,

    Match =
      Difference ==
      0
  )


cat(
  "\n============================================================\n"
)

cat(
  "FINAL SEN VERTICAL SNP VALIDATION\n"
)

cat(
  "============================================================\n"
)

print(
  final_validation,
  n = Inf
)


expected_total <- sum(
  expected_final$Expected_SNPs
)


observed_total <- nrow(
  final_snp
)


cat(
  "\nEXPECTED TOTAL:",
  expected_total,
  "\n"
)

cat(
  "OBSERVED TOTAL:",
  observed_total,
  "\n"
)


if (

  !all(
    final_validation$Match
  ) ||

  observed_total !=
  419
) {

  stop(
    "FINAL SEN SNP VALIDATION FAILED"
  )
}


cat(
  "\n✅ EXACT 419-SNP MATCH ACROSS ALL CLADE × COVERAGE BINS\n"
)


# ============================================================
# FINAL OUTPUT DIRECTORY
# ============================================================

final_outdir <- file.path(

  fga_dir,

  "FINAL_Combined_Split_with_unassigned_vertical_only"
)


dir.create(

  final_outdir,

  showWarnings = FALSE,

  recursive = TRUE
)


# ============================================================
# SAVE FINAL SNP VALIDATION
# ============================================================

write.table(

  final_validation,

  file.path(

    final_outdir,

    "FINAL_SNP_VALIDATION.tsv"
  ),

  sep =
    "\t",

  row.names =
    FALSE,

  quote =
    FALSE
)


# ============================================================
# SAVE FINAL 419-SNP MASTER
# ============================================================

write.csv(

  final_snp,

  file.path(

    final_outdir,

    "FINAL_419_SNP_MASTER.csv"
  ),

  row.names =
    FALSE
)


write.table(

  final_snp,

  file.path(

    final_outdir,

    "FINAL_419_SNP_MASTER.tsv"
  ),

  sep =
    "\t",

  row.names =
    FALSE,

  quote =
    FALSE
)


# ============================================================
# SAVE FINAL COUNTS
# ============================================================

write.csv(

  final_counts,

  file.path(

    final_outdir,

    "counts.csv"
  ),

  row.names =
    FALSE
)


write.table(

  final_counts,

  file.path(

    final_outdir,

    "counts.tsv"
  ),

  sep =
    "\t",

  row.names =
    FALSE,

  quote =
    FALSE
)


# ============================================================
# FINAL PLOT
# ============================================================

final_present_cogs <- sort(

  unique(

    as.character(

      final_counts$COG_label
    )
  )
)


final_present_cogs <-
  final_present_cogs[
    !is.na(
      final_present_cogs
    )
  ]


final_cog_palette <-
  make_distinct_palette(
    final_present_cogs
  )


p_final <- ggplot(

  final_counts,

  aes(

    x =
      Clade,

    y =
      n,

    fill =
      COG_label
  )
) +

  geom_bar(

    stat =
      "identity",

    position =
      "stack",

    width =
      0.7,

    color =
      "black",

    linewidth =
      0.25
  ) +

  coord_flip() +

  facet_grid(

    Impact_Class ~
      Coverage_Bin,

    scales =
      "free_x"
  ) +

  scale_fill_manual(

    values =
      final_cog_palette,

    name =
      "COG",

    drop =
      TRUE
  ) +

  labs(

    x =
      "Clade",

    y =
      "Defining SNP × COG assignments"
  ) +

  theme_minimal(
    base_size =
      20
  ) +

  theme(

    axis.text.x =
      element_text(
        size = 20
      ),

    axis.text.y =
      element_text(
        size = 20
      ),

    axis.title =
      element_text(
        size = 24
      ),

    legend.text =
      element_text(
        size = 20
      ),

    legend.title =
      element_text(
        size = 20
      ),

    legend.key.size =
      unit(
        0.8,
        "cm"
      ),

    strip.text =
      element_text(
        size = 18,
        face = "bold"
      ),

    panel.border =
      element_rect(
        color = "black",
        fill = NA,
        linewidth = 1.2
      ),

    panel.spacing =
      unit(
        0.6,
        "lines"
      )
  )


# ============================================================
# SAVE FINAL PLOTS
# ============================================================

ggsave(

  file.path(

    final_outdir,

    "plot.png"
  ),

  p_final,

  width =
    28,

  height =
    18,

  dpi =
    600
)


ggsave(

  file.path(

    final_outdir,

    "plot.tif"
  ),

  p_final,

  width =
    28,

  height =
    18,

  dpi =
    600,

  compression =
    "lzw"
)


ggsave(

  file.path(

    final_outdir,

    "plot.pdf"
  ),

  p_final,

  width =
    28,

  height =
    18
)


write.table(

  tibble(

    File_Stem =
      "plot",

    Plot_Title =
      paste0(

        "Final COG distribution of vertical clade-defining SNPs for ",

        "Clade 1, Clade 1A, Clade 1B, Clade 2, Clade 3, and Unassigned"
      )
  ),

  file.path(

    final_outdir,

    "plot_TITLE.tsv"
  ),

  sep =
    "\t",

  row.names =
    FALSE,

  quote =
    FALSE
)


cat(
  "\n============================================================\n"
)

cat(
  "✅ ALL 8 ITERATIONS COMPLETE\n"
)

cat(
  "✅ 8 SEPARATE MASTER PACKAGES CREATED\n"
)

cat(
  "✅ MULTI-SNP GENES EXPANDED CORRECTLY\n"
)

cat(
  "✅ FINAL 419-SNP VALIDATION PASSED\n"
)

cat(
  "============================================================\n"
)

cat(
  "\nMaster folder:\n",
  master_root,
  "\n"
)

