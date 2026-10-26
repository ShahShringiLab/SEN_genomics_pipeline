#!/usr/bin/env Rscript
# Faithful original SEN Rmd 3-panel layout; full-master 15,456-positive inputs.
# Reads existing full-master diversity CSVs: does not rerun bootstrap/permutations.
source(file.path(Sys.getenv("SEN_ROOT", unset=getwd()), "analysis/part_a/lib/common.R"))
suppressPackageStartupMessages(library(patchwork))
outdir <- sen_path("results/part_a/master_33831/diversity_full")
div_sum <- readr::read_csv(file.path(outdir,"QR_diversity_q0_q1_q2_bySource_with95CI.csv"),show_col_types=FALSE)
comp <- readr::read_csv(file.path(outdir,"QR_composition_bySource.csv"),show_col_types=FALSE)
source_levels <- c("Human","Chicken-US","Other-US","International")
div_sum$Source <- factor(div_sum$Source, levels=source_levels)
comp$Source <- factor(comp$Source, levels=source_levels)
# Historical plot letters were manually entered without inference; do not assert
# those letters on updated statistics unless independently revalidated.

div_long_plot <- div_sum %>%
  select(
    Source,
    q1, q1_lo, q1_hi,
    q2, q2_lo, q2_hi
  ) %>%
  pivot_longer(
    cols = c(q1, q2),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    lo = if_else(
      metric == "q1",
      q1_lo,
      q2_lo
    ),

    hi = if_else(
      metric == "q1",
      q1_hi,
      q2_hi
    ),

    metric = recode(
      metric,
      q1 = "Exponential Shannon (q1)",
      q2 = "Inverse Simpson (q2)"
    ),

    metric = factor(
      metric,
      levels = c(
        "Exponential Shannon (q1)",
        "Inverse Simpson (q2)"
      )
    )
  ) %>%
  mutate(
    Source = factor(
      as.character(Source),
      levels = source_levels
    )
  )


# =========================================================
# 1) EXPONENTIAL SHANNON
# =========================================================

p_shannon <- div_long_plot %>%
  filter(
    metric == "Exponential Shannon (q1)"
  ) %>%

  ggplot(
    aes(
      x = Source,
      y = value
    )
  ) +

  geom_point(
    size = 2.8
  ) +

  geom_errorbar(
    aes(
      ymin = lo,
      ymax = hi
    ),
    width = 0.15
  ) +

coord_cartesian(
    ylim = c(0, 7)
  ) +

  labs(
    title = "Exponential Shannon (q1)",
    x = "Source",
    y = "Effective number of QR/FQR categories"
  ) +

  theme_bw() +

  theme(
    plot.title = element_text(
      face = "bold",
      size = 11,
      hjust = 0.5
    ),

    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold",
      size = 10
    ),

    axis.text.y = element_text(
      face = "bold",
      size = 10
    ),

    axis.title = element_text(
      face = "bold",
      size = 11
    ),

    plot.margin = margin(
      6, 6, 6, 6
    )
  )


# =========================================================
# 2) INVERSE SIMPSON
# =========================================================

p_simpson <- div_long_plot %>%
  filter(
    metric == "Inverse Simpson (q2)"
  ) %>%

  ggplot(
    aes(
      x = Source,
      y = value
    )
  ) +

  geom_point(
    size = 2.8
  ) +

  geom_errorbar(
    aes(
      ymin = lo,
      ymax = hi
    ),
    width = 0.15
  ) +

coord_cartesian(
    ylim = c(0, 7)
  ) +

  labs(
    title = "Inverse Simpson (q2)",
    x = "Source",
    y = "Effective number of QR/FQR categories"
  ) +

  theme_bw() +

  theme(
    plot.title = element_text(
      face = "bold",
      size = 11,
      hjust = 0.5
    ),

    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold",
      size = 10
    ),

    axis.text.y = element_text(
      face = "bold",
      size = 10
    ),

    axis.title = element_text(
      face = "bold",
      size = 11
    ),

    plot.margin = margin(
      6, 6, 6, 6
    )
  )


# =========================================================
# 3) QR/FQR COMPOSITION
# =========================================================

p_comp <- ggplot(
  comp,
  aes(
    x = Source,
    y = p,
    fill = marker_group
  )
) +

  geom_col(
    width = 0.85
  ) +

  # Freeze the original manuscript color-to-marker mapping (not ggplot's
  # automatic factor-order palette, which can change after a full-cohort rerun).
  scale_fill_manual(
    values = c(
      "gyrA_D87G" = "#F8766D",
      "gyrA_D87C" = "#F8766D",
      "gyrA_D87Y" = "#B79F00",
      "gyrA_S83Y" = "#00BA38",
      "gyrA_D87N" = "#00BFC4",
      "others"   = "#619CFF",
      "qnrB19"   = "#F564E3"
    ),
    breaks = c("gyrA_D87G", "gyrA_D87C", "gyrA_D87Y",
               "gyrA_S83Y", "gyrA_D87N", "others", "qnrB19"),
    drop = TRUE
  ) +

  scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    ),
    expand = expansion(
      mult = c(
        0,
        0.02
      )
    )
  ) +

  labs(
    title = "QR/FQR composition",
    x = "Source",
    y = "Normalized QR/FQR % composition",
    fill = "QR category"
  ) +

  theme_bw() +

  theme(
    plot.title = element_text(
      face = "bold",
      size = 11,
      hjust = 0.5
    ),

    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold",
      size = 9
    ),

    axis.text.y = element_text(
      face = "bold",
      size = 10
    ),

    axis.title = element_text(
      face = "bold",
      size = 11
    ),

    legend.title = element_text(
      face = "bold",
      size = 11
    ),

    legend.text = element_text(
      size = 9
    ),

    legend.key.size = unit(
      0.45,
      "cm"
    ),

    plot.margin = margin(
      6, 6, 6, 6
    )
  )


# =========================================================
# 4) PATCHWORK
#
# SAME HEIGHT for all three panels
#
# Relative widths:
#   Shannon = 1.10
#   Simpson = 1.10
#   Composition = 1.00
#
# Thus Shannon and Simpson are 10% wider.
# =========================================================

combo_side <- (
  p_shannon |
  p_simpson |
  p_comp
) +
  plot_layout(
    widths = c(
      1.10,
      1.10,
      1.00
    )
  )


# =========================================================
# PRINT
# =========================================================

print(
  combo_side
)


# =========================================================
# SAVE
# =========================================================

ggsave(
  file.path(outdir, "QR_diversity_and_composition_side_by_side.png"),
  combo_side,
  width = 18,
  height = 6.2,
  dpi = 600
)

ggsave(
  file.path(outdir, "QR_diversity_and_composition_side_by_side.tif"),
  combo_side,
  width = 18,
  height = 6.2,
  dpi = 600,
  compression = "lzw"
)

cat("[PASS] Original SEN three-panel plot regenerated from full-master diversity tables at",outdir,"\\n")
