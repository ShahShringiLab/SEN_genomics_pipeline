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
