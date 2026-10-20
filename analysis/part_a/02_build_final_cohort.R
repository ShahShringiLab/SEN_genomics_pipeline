#!/usr/bin/env Rscript
source(file.path(
  Sys.getenv("SEN_ROOT", unset = getwd()),
  "analysis/part_a/lib/common.R"
))

master_file <- sen_path("data/part_a/genotype/US-SEN_12192025.csv")
key_file <- sen_path("data/part_a/selection/Selected_updated_Isolate.csv")
canonical_file <- sen_path("data/part_a/selection/Finalupdated_SEN_Genomes1.csv")
clade_file <- sen_path("metadata/final_clade_metadata.tsv")

outdir <- ensure_dir(sen_path("results/part_a/cohort"))
out_file <- file.path(outdir, "Finalupdated_SEN_Genomes_rebuilt.csv")

master <- read_sen_master_csv(master_file)
keys <- read_csv(key_file, show_col_types = FALSE, progress = FALSE)
canonical <- read_csv(canonical_file, show_col_types = FALSE, progress = FALSE)
clades <- read_tsv(clade_file, show_col_types = FALSE, progress = FALSE)

required_master <- c(
  "Run", "Strain", "Isolate", "AMR genotypes core", "AMR genotypes",
  "Collection date", "Source type", "Isolation source", "Host",
  "IFSAC category", "Collected by", "Isolate identifiers", "Food origin"
)
assert_columns(master, required_master, "master genotype table")
assert_columns(keys, c("Strain", "Isolate"), "selection key")
assert_columns(
  canonical,
  c(
    "Run", "Isolate", "QR_genotype", "Collection_year", "Source",
    "Origin_group", "Source_group"
  ),
  "canonical final cohort"
)
assert_columns(clades, c("Sample", "Lineage"), "final clade metadata")

if (nrow(master) != 33831L || n_distinct(master$Isolate) != 33831L) {
  stop("Master table must contain 33,831 unique Isolate records.")
}
if (nrow(keys) != 3434L || n_distinct(keys$Isolate) != 3434L) {
  stop("Selection key must contain 3,434 unique Isolate records.")
}
if (nrow(canonical) != 3434L || n_distinct(canonical$Isolate) != 3434L) {
  stop("Canonical final cohort must contain 3,434 unique Isolate records.")
}

keys <- keys %>%
  mutate(across(c(Strain, Isolate), ~ str_squish(as.character(.x))))
master <- master %>%
  mutate(across(c(Strain, Isolate), ~ str_squish(as.character(.x))))

selected <- master %>%
  semi_join(keys, by = c("Strain", "Isolate"))

if (nrow(selected) != 3434L) {
  stop("Strain+Isolate selection produced ", nrow(selected), " rows; expected 3,434.")
}
if (n_distinct(selected$Isolate) != 3434L) {
  stop("Selected cohort contains duplicate Isolate values.")
}
if (n_distinct(selected$Run) != 3434L) {
  stop("Selected cohort contains duplicate or missing Run accessions.")
}

# Rebuild fields from the production rules.
selected$Collection_year_recomputed <- derive_collection_year(
  selected[["Collection date"]]
)

# QR/FQR status is derived from AMRFinderPlus core genotype calls rather than
# the broader CSV field. The latter is vulnerable to a known quoting artifact
# in three source records. This also keeps the binary call tied to curated
# QR/FQR determinant families.
selected$QR_genotype_recomputed <- derive_qr_genotype_core(
  selected[["AMR genotypes core"]]
)
selected$QR_markers <- collapse_qr_core_markers(
  selected[["AMR genotypes core"]]
)
selected$QR_marker_count <- ifelse(
  selected$QR_markers == "",
  0L,
  lengths(strsplit(selected$QR_markers, ";", fixed = TRUE))
)

primary_source <- classify_primary_source(
  selected[["Source type"]],
  selected[["Isolation source"]],
  selected[["Host"]],
  selected[["IFSAC category"]],
  selected[["Collected by"]]
)

origin <- derive_origin_and_source(
  primary_source,
  selected[["Strain"]],
  selected[["Isolate identifiers"]],
  selected[["Food origin"]]
)

selected$Origin_group_recomputed <- origin$Origin_group
selected$Source_group_recomputed <- origin$Source_group
selected$Source_recomputed <- origin$Source

# Compare rebuilt fields against the archived 3,434-genome cohort.
comparison <- selected %>%
  select(
    Run, Isolate, QR_genotype_recomputed, Collection_year_recomputed,
    Source_recomputed, Origin_group_recomputed, Source_group_recomputed
  ) %>%
  inner_join(
    canonical %>%
      select(
        Run, Isolate, QR_genotype, Collection_year, Source,
        Origin_group, Source_group
      ),
    by = c("Run", "Isolate")
  )

if (nrow(comparison) != 3434L) {
  stop("Rebuilt cohort does not match the canonical Run+Isolate set.")
}

comparison <- comparison %>%
  mutate(
    QR_mismatch =
      as.integer(QR_genotype_recomputed) != as.integer(QR_genotype),
    year_mismatch =
      as.integer(Collection_year_recomputed) != as.integer(Collection_year),
    source_mismatch =
      as.character(Source_recomputed) != as.character(Source),
    origin_mismatch =
      as.character(Origin_group_recomputed) != as.character(Origin_group),
    source_group_mismatch =
      as.character(Source_group_recomputed) != as.character(Source_group)
  )

write_csv(
  comparison,
  file.path(outdir, "cohort_rebuild_validation_all.csv")
)

qr_corrections <- comparison %>%
  filter(QR_mismatch) %>%
  arrange(Isolate)

write_csv(
  qr_corrections,
  file.path(outdir, "QR_genotype_corrections.csv")
)

if (any(comparison$year_mismatch, na.rm = TRUE)) {
  stop("Collection-year reconstruction differs from the archived final cohort.")
}
if (
  any(comparison$source_mismatch, na.rm = TRUE) ||
  any(comparison$origin_mismatch, na.rm = TRUE) ||
  any(comparison$source_group_mismatch, na.rm = TRUE)
) {
  stop("Source/origin reconstruction differs from the archived final cohort.")
}

# These three records are true QR/FQR-marker-positive genomes. Their upstream
# AMRFinderPlus core field contains gyrA_D87Y=POINT, while the archived
# QR_genotype value is 0 because the malformed source CSV row shifted the
# broader AMR-genotype field during the historical parse.
expected_qr_corrections <- c(
  "PDT001102834.1",
  "PDT001168391.1",
  "PDT001168656.1"
)
observed_qr_corrections <- sort(as.character(qr_corrections$Isolate))
if (!identical(observed_qr_corrections, sort(expected_qr_corrections))) {
  stop(
    "Unexpected QR-genotype correction set. Observed: ",
    paste(observed_qr_corrections, collapse = ", ")
  )
}

# Preserve the archived QR value for provenance, but use the recomputed
# AMRFinderPlus-core call for all production genotype analyses.
canonical_fields <- canonical %>%
  select(
    Run, Isolate, QR_genotype, Collection_year,
    Source, Origin_group, Source_group
  ) %>%
  rename(
    QR_genotype_archived = QR_genotype,
    Collection_year_canonical = Collection_year,
    Source_canonical = Source,
    Origin_group_canonical = Origin_group,
    Source_group_canonical = Source_group
  )

selected <- selected %>%
  left_join(canonical_fields, by = c("Run", "Isolate")) %>%
  mutate(
    QR_genotype = as.integer(QR_genotype_recomputed),
    FQR_genotype = QR_genotype,
    QR_status = if_else(
      QR_genotype == 1L,
      "QR/FQR marker-positive",
      "QR/FQR marker-negative"
    ),
    QR_corrected = QR_genotype != as.integer(QR_genotype_archived),
    Collection_year = as.integer(Collection_year_canonical),
    Source = as.character(Source_canonical),
    Origin_group = as.character(Origin_group_canonical),
    Source_group = as.character(Source_group_canonical),
    Source3 = case_when(
      str_detect(str_to_lower(Source), "human") ~ "Human",
      str_detect(str_to_lower(Source), "chicken") ~ "Chicken",
      TRUE ~ "Other"
    ),
    NonUS = Source == "International",
    QR_group = if_else(QR_genotype == 1L, "QR", "NON-QR"),
    AMR_pattern = assign_amr_pattern(.data[["AMR genotypes core"]])
  ) %>%
  select(
    -Collection_year_canonical,
    -Source_canonical,
    -Origin_group_canonical,
    -Source_group_canonical
  )

# Attach the authoritative Part B clade labels by SRA Run.
clade_map <- clades %>%
  filter(Sample != "Reference") %>%
  transmute(
    Run = as.character(Sample),
    Clade = as.character(Lineage)
  )

selected <- selected %>%
  left_join(clade_map, by = "Run") %>%
  mutate(
    Final_phylogeny_status = if_else(
      is.na(Clade),
      "Excluded before final phylogeny",
      "Included in final phylogeny"
    )
  )

clade_counts <- selected %>%
  filter(!is.na(Clade)) %>%
  count(Clade, name = "n") %>%
  arrange(Clade)

expected_clades <- c(
  "Clade 1A" = 169L,
  "Clade 1B" = 759L,
  "Clade 2" = 894L,
  "Clade 3" = 1357L,
  "Unassigned" = 127L
)
observed_clades <- setNames(clade_counts$n, clade_counts$Clade)
if (!identical(
  as.integer(observed_clades[names(expected_clades)]),
  as.integer(expected_clades)
)) {
  stop("Final clade counts do not match canonical metadata.")
}
if (sum(!is.na(selected$Clade)) != 3306L) {
  stop("Expected 3,306 selected isolates with final clade labels.")
}
if (sum(selected$QR_genotype == 1L) != 2288L) {
  stop("Expected 2,288 QR/FQR marker-positive isolates after core-field correction.")
}

write_csv(selected, out_file)
write_csv(clade_counts, file.path(outdir, "final_clade_counts.csv"))

summary_table <- tibble(
  check = c(
    "selected_rows",
    "unique_isolates",
    "unique_runs",
    "collection_year_mismatches",
    "source_mismatches",
    "origin_mismatches",
    "source_group_mismatches",
    "QR_genotype_corrections",
    "QR_marker_positive",
    "QR_marker_negative",
    "final_phylogeny_clade_mapped",
    "final_phylogeny_excluded"
  ),
  value = c(
    nrow(selected),
    n_distinct(selected$Isolate),
    n_distinct(selected$Run),
    sum(comparison$year_mismatch, na.rm = TRUE),
    sum(comparison$source_mismatch, na.rm = TRUE),
    sum(comparison$origin_mismatch, na.rm = TRUE),
    sum(comparison$source_group_mismatch, na.rm = TRUE),
    nrow(qr_corrections),
    sum(selected$QR_genotype == 1L),
    sum(selected$QR_genotype == 0L),
    sum(!is.na(selected$Clade)),
    sum(is.na(selected$Clade))
  )
)

write_csv(
  summary_table,
  file.path(outdir, "cohort_rebuild_validation_summary.csv")
)
print(summary_table)
print(clade_counts)

cat(
  "[PASS] Part A2 cohort validated: 3,434 isolates; ",
  "2,288 QR/FQR marker-positive; 3 QR calls corrected from AMRFinderPlus core; ",
  "3,306 final clade labels attached.\n",
  sep = ""
)
