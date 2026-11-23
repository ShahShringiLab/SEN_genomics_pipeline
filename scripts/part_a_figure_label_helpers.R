# Publication-only label formatting for SEN Part A figures.
# No data recoding, plot geometry, palettes, or statistical calculations.
#
# Figure 2: mixed italic gene/roman mutation labels, and italic n.
# Figure 3, S2, S3: italic n within source facet strips, em dash in
# source names; numeric percentage-bin ranges use em dash.

sen_marker_plotmath <- function(labels, bold = FALSE) {
  vapply(as.character(labels), function(x) {
    if (x %in% c("Q/FQ-susceptible", "Q/FQ-Susceptible", "QS/FQS")) {
      return(if (bold) 'bold("QS/FQS")' else '"QS/FQS"')
    }
    if (x %in% c("others", "Others", "QR/FQR")) {
      return(sprintf('"%s"', x))
    }
    parts <- strsplit(x, "_", fixed = TRUE)[[1]]
    gene <- parts[1]
    # The source tokens contain gene symbols followed by optional mutation
    # suffixes. The mutation suffix itself remains upright.
    if (!grepl("^[A-Za-z][A-Za-z0-9()'-]*$", gene)) {
      stop("Unexpected gene token for manuscript axis: ", x)
    }
    base <- sprintf('%s("%s")',
                    if (bold) "bolditalic" else "italic", gene)
    if (length(parts) == 1L) return(base)
    rest <- paste0("_", paste(parts[-1], collapse = "_"))
    suffix <- sprintf('%s("%s")', if (bold) "bold" else "plain", rest)
    sprintf('paste(%s, %s)', base, suffix)
  }, character(1))
}

sen_pctbin_labels <- function(labels) {
  # Display-only; keep the original numeric cut() boundaries and factor levels.
  gsub("([0-9])-([0-9])", "\\1\u2014\\2",
       as.character(labels), perl = TRUE)
}

sen_source_plotmath <- function(labels) {
  vapply(as.character(labels), function(x) {
    match <- regmatches(x, regexec("^(.+) \\(N=([0-9]+)\\)$", x))[[1]]
    if (length(match) != 3L) {
      stop("Unexpected source facet label: ", x)
    }
    group <- gsub("-", "\u2014", match[2], fixed = TRUE)
    count <- format(as.integer(match[3]), big.mark = ",",
                    trim = TRUE, scientific = FALSE)
    sprintf('paste("%s (", italic(n), " = %s)")', group, count)
  }, character(1))
}

sen_source_labeller <- ggplot2::labeller(
  Source_lab = ggplot2::as_labeller(
    sen_source_plotmath, default = ggplot2::label_parsed
  )
)

# Apply publication punctuation at the export boundary. This operates only on
# visible text in plot specifications, never on sample metadata, grouping,
# breakpoints, factor levels, coordinates, or statistical results.
sen_em_dash_text <- function(x) {
  if (is.character(x)) return(gsub("\u2013", "\u2014", x, fixed = TRUE))
  x
}

sen_em_dash_plot <- function(plot) {
  if (!inherits(plot, "ggplot")) return(plot)
  if (!is.null(plot$labels)) {
    for (key in names(plot$labels)) {
      plot$labels[[key]] <- sen_em_dash_text(plot$labels[[key]])
    }
  }
  for (sc in plot$scales$scales) {
    labels <- sc$labels
    if (is.character(labels)) {
      sc$labels <- sen_em_dash_text(labels)
    } else if (is.function(labels)) {
      original <- labels
      sc$labels <- local({
        fn <- original
        function(...) sen_em_dash_text(fn(...))
      })
    }
  }
  plot
}

# Mask unqualified ggsave() in generated Part A R execution only. The TIFF
# companion covers original plots historically exported as PNG/PDF alone.
# The original explicit TIFF export is still used when present.
sen_ggsave_part_a <- function(filename, plot = ggplot2::last_plot(), ...) {
  displayed <- sen_em_dash_plot(plot)
  result <- ggplot2::ggsave(filename = filename, plot = displayed, ...)
  if (grepl("\\.(png|pdf)$", filename, ignore.case = TRUE)) {
    target <- sub("\\.(png|pdf)$", ".tif", filename, ignore.case = TRUE)
    if (!file.exists(target)) {
      args <- list(...)
      args$dpi <- 600
      args$compression <- "lzw"
      args$device <- "tiff"
      do.call(ggplot2::ggsave,
              c(list(filename = target, plot = displayed), args))
    }
  }
  invisible(result)
}
