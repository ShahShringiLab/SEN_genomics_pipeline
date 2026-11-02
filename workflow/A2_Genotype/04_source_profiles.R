#!/usr/bin/env Rscript
source(file.path(
  Sys.getenv("SEN_ROOT", unset = getwd()),
  "analysis/part_a/lib/common.R"
))

input <- sen_path("results/part_a/cohort/Finalupdated_SEN_Genomes1_QC_merged.csv")
if (!file.exists(input)) {
  stop("Run python3 scripts/merge_qc_full_metadata.py first.")
}

outdir <- ensure_dir(sen_path("results/part_a/source_profiles"))
sen <- read_csv(input, show_col_types = FALSE, progress = FALSE)

assert_columns(
  sen,
  c(
    "Isolate", "Source", "Source_group", "Isolation source",
    "QR_genotype", "Food origin"
  ),
  "QC-merged 3,434-genome cohort"
)

if (nrow(sen) != 3434L || n_distinct(sen$Isolate) != 3434L) {
  stop("Source profiling requires the validated 3,434-isolate cohort.")
}

source_counts <- sen %>%
  count(Source, name = "n", sort = TRUE) %>%
  mutate(pct = 100 * n / sum(n))

write_csv(source_counts, file.path(outdir, "SEN_Source_counts.csv"))

food_origin_counts <- sen %>%
  mutate(
    Food_origin = str_squish(as.character(.data[["Food origin"]])),
    Food_origin = if_else(
      is.na(Food_origin) | Food_origin == "",
      "Unknown",
      Food_origin
    )
  ) %>%
  count(Food_origin, name = "n", sort = TRUE)

write_csv(
  food_origin_counts,
  file.path(outdir, "SEN_Food_origin_counts.csv")
)

# Refine records whose pre-collapse source group is Other-* using the same
# animal-source vocabulary used in the original Rmd.
other_refined <- sen %>%
  filter(str_detect(as.character(Source_group), "^Other-")) %>%
  transmute(
    Isolate = as.character(Isolate),
    Source_group = as.character(Source_group),
    Other_source_group = classify_animal_source(.data[["Isolation source"]])
  )

other_refined_counts <- other_refined %>%
  count(Source_group, Other_source_group, name = "n") %>%
  group_by(Source_group) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

write_csv(
  other_refined_counts,
  file.path(outdir, "SEN_Other_source_refined_counts.csv")
)

# Reproduce the manuscript-facing International and Other-US animal-source
# decomposition from the final source framework.
source_animal_qr <- sen %>%
  filter(Source %in% c("International", "Other-US")) %>%
  transmute(
    Isolate = as.character(Isolate),
    Group = as.character(Source),
    Animal_source = classify_animal_source(.data[["Isolation source"]]),
    QR_status = if_else(
      as.integer(QR_genotype) == 1L,
      "QR/FQR marker",
      "Q/FQ-susceptible"
    )
  ) %>%
  distinct()

group_totals <- source_animal_qr %>%
  count(Group, name = "n")

write_csv(
  group_totals,
  file.path(outdir, "SEN_International_OtherUS_group_totals.csv")
)

animal_original_counts <- source_animal_qr %>%
  count(Group, Animal_source, name = "n") %>%
  group_by(Group) %>%
  mutate(
    group_total = sum(n),
    pct = 100 * n / group_total
  ) %>%
  ungroup()

write_csv(
  animal_original_counts,
  file.path(outdir, "SEN_International_OtherUS_animal_original_counts.csv")
)

animal_qr_counts <- source_animal_qr %>%
  count(Group, Animal_source, QR_status, name = "n") %>%
  group_by(Group, Animal_source) %>%
  mutate(
    total_source = sum(n),
    pct = 100 * n / total_source
  ) %>%
  ungroup()

write_csv(
  animal_qr_counts,
  file.path(outdir, "SEN_International_OtherUS_animal_QRFQR_counts.csv")
)

audit <- source_animal_qr %>%
  group_by(Group) %>%
  summarise(
    Total_isolates = n_distinct(Isolate),
    QR_FQR_marker = sum(QR_status == "QR/FQR marker"),
    Q_FQ_susceptible = sum(QR_status == "Q/FQ-susceptible"),
    .groups = "drop"
  )

write_csv(
  audit,
  file.path(outdir, "SEN_International_OtherUS_QRFQR_final_audit.csv")
)

expected_source <- c(
  "Human" = 1976L,
  "Chicken-US" = 1059L,
  "Other-US" = 328L,
  "International" = 71L
)

observed_source <- setNames(source_counts$n, source_counts$Source)
if (!identical(
  as.integer(observed_source[names(expected_source)]),
  as.integer(expected_source)
)) {
  stop(
    "Final source counts differ from the release contract. Observed: ",
    paste(names(observed_source), observed_source, sep = "=", collapse = ", ")
  )
}

plot_levels <- c("Chicken", "Turkey", "Cattle", "Sheep/goat", "Swine", "Others")

plot_df <- animal_qr_counts %>%
  mutate(
    Animal_source = factor(Animal_source, levels = plot_levels),
    QR_status = factor(
      QR_status,
      levels = c("Q/FQ-susceptible", "QR/FQR marker")
    )
  )

p <- ggplot(plot_df, aes(Animal_source, n, fill = QR_status)) +
  geom_col(width = 0.75) +
  facet_wrap(~ Group, scales = "free_y") +
  scale_fill_manual(
    values = c(
      "Q/FQ-susceptible" = "green4",
      "QR/FQR marker" = "red3"
    ),
    drop = FALSE
  ) +
  labs(
    x = "Animal source",
    y = "Number of isolates",
    fill = ""
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(face = "bold")
  )

ggsave(
  file.path(outdir, "SEN_International_OtherUS_animal_QRFQR.png"),
  p,
  width = 9,
  height = 5,
  dpi = 600
)
ggsave(
  file.path(outdir, "SEN_International_OtherUS_animal_QRFQR.pdf"),
  p,
  width = 9,
  height = 5
)

cat("[PASS] Part A2 source/origin profiles written to", outdir, "\n")
