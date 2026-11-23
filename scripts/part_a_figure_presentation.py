#!/usr/bin/env python3
"""Apply ONLY publication-requested display formatting to original Part A plots.

All source metadata, mutations, counts, bin boundaries, colours, plot sizes,
geometries, and scientific results stay unchanged. The SHA-frozen original Rmd
is never edited. Patches are applied to execution copies of chunks 6, 19-21.
"""
from __future__ import annotations


def one(src: str, old: str, new: str, count: int = 1, where: str = "") -> str:
    actual = src.count(old)
    if actual != count:
        raise RuntimeError(
            f"Refusing unexpected Figure label patch {where}: "
            f"expected {count} anchor(s), found {actual}: {old[:100]!r}"
        )
    return src.replace(old, new, count)


def patch_chunk(index: int, chunk: str) -> str:
    if index == 6:
        # Figure 2: display ONLY. The count labels remain attached to their
        # original records; original marker frequencies and Type stay intact.
        chunk = one(chunk, "print(p_tokens)", r'''
# Figure 2 manuscript display-only corrections.
stopifnot(sum(plot_df$n_isolates == 18375L) == 1L,
          sum(plot_df$n_isolates == 13486L) == 1L)
p_tokens$data$label_plot <- sub(
  "^18375 \\(", "18,375 (", p_tokens$data$label_plot
)
p_tokens$data$label_plot <- sub(
  "^13486 \\(", "13,486 (", p_tokens$data$label_plot
)
p_tokens <- p_tokens +
  ggplot2::scale_x_discrete(
    labels = function(x) parse(text = sen_marker_plotmath(x, bold = TRUE))
  ) +
  ggplot2::labs(y = expression(paste("Genomes(", italic(n), " (%))")))
print(p_tokens)''', where="Fig2 axis and two counts")
        return chunk

    if index in (19, 20, 21):
        # Facet re-labelling: Human (italic n = 26,296), Chicken—US, etc.
        # Existing Source_lab factor (and source calculations) are untouched.
        chunk = one(
            chunk,
            "~ Source_lab,\n",
            "~ Source_lab,\n    labeller = sen_source_labeller,\n",
            count=2, where=f"chunk {index} source facet labels",
        )

    if index == 19:
        # Figure 3 and companion 5-year plot: italic gene symbols only.
        # Keep the original 30% y-padding and original percent-bin cut points.
        chunk = one(
            chunk,
            "scale_y_discrete(\n    expand = y_pad\n  )",
            "scale_y_discrete(\n    expand = y_pad,\n"
            "    labels = function(x) parse(text = sen_marker_plotmath(x))\n  )",
            count=2, where="Fig3 and 5yr gene labels",
        )
        chunk = one(
            chunk,
            "print(p_year_binned)\nprint(p_5_binned)",
            """# Change only the display of percent-bin labels (1-5 => 1—5).
p_year_binned <- p_year_binned +
  ggplot2::scale_color_discrete(labels = sen_pctbin_labels)
p_5_binned <- p_5_binned +
  ggplot2::scale_color_discrete(labels = sen_pctbin_labels)
print(p_year_binned)
print(p_5_binned)""",
            where="Fig3 percent-bin legends",
        )
        return chunk

    if index == 21:
        # Supplementary expanded-marker count and percent plots. Keep the
        # original two independent scales, seven-category order and colours.
        chunk = one(
            chunk,
            "    breaks = expanded_marker_levels,\n\n    drop = FALSE",
            "    breaks = expanded_marker_levels,\n\n"
            "    labels = function(x) parse(text = sen_marker_plotmath(x)),\n\n"
            "    drop = FALSE",
            count=2, where="S3 count and percent marker legends",
        )

    if index == 22:
        # Three-panel Hill diversity/composition figure. Apply only displayed
        # labels, leaving the bootstrap/permutation calculations unchanged.
        # Patch the three constituent ggplot objects before patchwork joins
        # them; modifying the resulting patchwork object is insufficient.
        chunk = one(
            chunk,
            "combo_side <- (",
            """# Publication typography on each constituent plot.
# Keep Source factor values and order unchanged; display em dash only.
sen_source_display <- function(x) gsub("-", "\\u2014", x, fixed = TRUE)
p_shannon <- p_shannon +
  ggplot2::scale_x_discrete(labels = sen_source_display)
p_simpson <- p_simpson +
  ggplot2::scale_x_discrete(labels = sen_source_display)
p_comp <- p_comp +
  ggplot2::scale_x_discrete(labels = sen_source_display) +
  ggplot2::scale_fill_discrete(
    labels = function(x) parse(text = sen_marker_plotmath(x))
  )
combo_side <- (""",
            where="Hill diversity and composition display labels",
        )
        return chunk

    return chunk
