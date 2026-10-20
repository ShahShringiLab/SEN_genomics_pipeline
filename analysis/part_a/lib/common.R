options(stringsAsFactors = FALSE)

required_pkgs <- c("dplyr", "tidyr", "readr", "stringr", "ggplot2", "lubridate")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop(
    "Missing R packages: ", paste(missing_pkgs, collapse = ", "),
    ". Create/update the sen_stats Conda environment first."
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(lubridate)
})

sen_root <- normalizePath(
  Sys.getenv("SEN_ROOT", unset = getwd()),
  winslash = "/",
  mustWork = TRUE
)

sen_path <- function(...) file.path(sen_root, ...)

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

assert_columns <- function(df, cols, label = "input") {
  missing <- setdiff(cols, names(df))
  if (length(missing)) {
    stop(label, " is missing required columns: ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

clean_text <- function(x) {
  x <- stringr::str_squish(as.character(x))
  x[is.na(x)] <- ""
  x
}

# The archived NCBI-derived master CSV contains 24 literal "\\," sequences in
# quoted Isolation-source text and one doubled-quote field around
# "Exudate, Left Heel". The source file is preserved byte-for-byte in data/.
# This function normalizes those parser artifacts only in memory, then requires
# the release contract of 33,831 rows and 51 columns.
read_sen_master_csv <- function(path) {
  if (!file.exists(path)) stop("Missing master genotype input: ", path)

  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  placeholder <- "__SEN_LITERAL_COMMA__"

  if (any(str_detect(lines, fixed(placeholder)))) {
    stop("Reserved parser placeholder already occurs in master input.")
  }

  escaped_comma_n <- sum(str_count(lines, fixed("\\,")))
  malformed_exudate_n <- sum(str_count(
    lines,
    fixed(',""Exudate, Left Heel"","clinical",')
  ))

  lines <- str_replace_all(lines, fixed("\\,"), placeholder)
  lines <- str_replace_all(
    lines,
    fixed(',""Exudate, Left Heel"","clinical",'),
    ',"Exudate, Left Heel","clinical",'
  )

  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(lines, tmp, useBytes = TRUE)

  dat <- read_csv(
    tmp,
    col_types = cols(.default = col_character()),
    show_col_types = FALSE,
    progress = FALSE
  )

  probs <- problems(dat)
  if (nrow(probs)) {
    print(probs)
    stop("Master genotype CSV produced parse problems after in-memory normalization.")
  }

  dat <- dat %>%
    mutate(across(
      where(is.character),
      ~ str_replace_all(.x, fixed(placeholder), ",")
    ))

  if ("Isolation source" %in% names(dat)) {
    dat[["Isolation source"]] <- str_replace(
      dat[["Isolation source"]],
      '""$',
      ""
    )
  }

  if (nrow(dat) != 33831L) {
    stop("Master genotype table rows=", nrow(dat), "; expected 33831.")
  }
  if (ncol(dat) != 51L) {
    stop("Master genotype table columns=", ncol(dat), "; expected 51.")
  }

  message(
    "[INFO] Master CSV parser normalization: ",
    escaped_comma_n, " escaped-comma sequence(s); ",
    malformed_exudate_n, " doubled-quote field(s)."
  )

  dat
}

fqr_marker_regex <- stringr::regex(
  paste0(
    "^(",
    "(gyr[a-b](?:_.+)?|par[c-e](?:_.+)?|grla(?:_.+)?)(?:$)|",
    "(qnr(as|[a-z0-9._-]*))(?:$)|",
    "(qnrvc[a-z0-9._-]*)(?:$)|",
    "(qep[a-z0-9._-]*)(?:$)|",
    "(nfxb(?:_.+)?|emrr(?:_.+)?|norm(?:_.+)?|pmra(?:_.+)?|acrb(?:_.+)?|smef(?:_.+)?)(?:$)|",
    "(aac\\(6'\\)-ib-cr(?:_.+)?|oqx[a-b](?:_.+)?|crpp(?:_.+)?|",
    "acr[a-b](?:_.+)?|tolc(?:_.+)?|ram[ar](?:_.+)?|mar[ar](?:_.+)?|sox[rs](?:_.+)?)(?:$)",
    ")"
  ),
  ignore_case = TRUE
)

# This narrower token rule reproduces the original notebook block that created
# FQR_genotype. Later marker-frequency analyses used the broader rule above.
fqr_genotype_token_regex <- stringr::regex(
  paste0(
    "^(",
    "gyr[a-b]_.+|par[c-e]_.+|",
    "qnr[a-z0-9._-]*|aac\\(6'\\)-ib-cr|qep[a-z0-9._-]*|oqx[a-b]|",
    "qnrvc[a-z0-9._-]*|crpp|",
    "acr[a-b]_.+|tolc_.+|ram[ar]_.+|mar[ar]_.+|sox[rs]_.+",
    ")$"
  ),
  ignore_case = TRUE
)

extract_amr_tokens <- function(amr) {
  txt <- ifelse(is.na(amr), "", as.character(amr))
  stringr::str_extract_all(txt, "[^,=\\s]+(?==)")
}

derive_fqr_genotype <- function(amr) {
  toks <- extract_amr_tokens(amr)
  vapply(
    toks,
    function(v) as.integer(any(str_detect(str_to_lower(v), fqr_genotype_token_regex))),
    integer(1)
  )
}


# Authoritative production QR/FQR calls are derived from the Pathogen Detection
# "AMR genotypes core" field. This avoids downstream-column shifts caused by
# malformed source CSV quoting in a small number of records and restricts the
# binary call to AMRFinderPlus-curated QR/FQR determinant families.
extract_qr_core_markers <- function(amr_core) {
  toks <- extract_amr_tokens(amr_core)
  lapply(
    toks,
    function(v) {
      z <- str_to_lower(v)
      keep <- str_detect(z, fqr_genotype_token_regex)
      unique(v[keep])
    }
  )
}

derive_qr_genotype_core <- function(amr_core) {
  markers <- extract_qr_core_markers(amr_core)
  vapply(markers, function(v) as.integer(length(v) > 0L), integer(1))
}

collapse_qr_core_markers <- function(amr_core) {
  markers <- extract_qr_core_markers(amr_core)
  vapply(
    markers,
    function(v) {
      if (!length(v)) return("")
      paste(sort(unique(pretty_token(v))), collapse = ";")
    },
    character(1)
  )
}

pretty_token <- function(x) {
  x <- str_to_lower(str_replace_all(str_trim(as.character(x)), "\\s+", ""))
  x <- str_replace_all(x, "-", "_")

  replacements <- c(
    "^gyra" = "gyrA", "^gyrb" = "gyrB",
    "^parc" = "parC", "^pard" = "parD", "^pare" = "parE", "^grla" = "grlA",
    "^qnrvc" = "qnrVC", "^qnra" = "qnrA", "^qnrb" = "qnrB",
    "^qnrc" = "qnrC", "^qnrd" = "qnrD", "^qnre" = "qnrE", "^qnrs" = "qnrS",
    "^qepa" = "qepA", "^oqxa" = "oqxA", "^oqxb" = "oqxB", "^crpp" = "crpP",
    "^acra" = "acrA", "^acrb" = "acrB", "^tolc" = "tolC",
    "^rama" = "ramA", "^ramr" = "ramR", "^mara" = "marA", "^marr" = "marR",
    "^soxr" = "soxR", "^soxs" = "soxS"
  )
  for (pat in names(replacements)) {
    x <- str_replace(x, pat, replacements[[pat]])
  }

  x <- str_replace(x, "^aac\\(6'\\)_ib_cr", "aac(6')-Ib-cr")

  parts <- str_split_fixed(x, "_", 2)
  ifelse(
    parts[, 2] == "",
    parts[, 1],
    paste0(parts[, 1], "_", toupper(parts[, 2]))
  )
}

derive_collection_year <- function(x) {
  z <- clean_text(x)
  z <- str_replace_all(z, "_", "-")
  z <- str_replace_all(z, "\\.", "/")

  parsed <- suppressWarnings(lubridate::parse_date_time(
    z,
    orders = c(
      "Ymd", "Y-m-d", "Y/m/d",
      "mdY", "m/d/Y", "m-d-Y",
      "dmY", "d/m/Y", "d-m-Y",
      "b d Y", "d b Y", "B d Y", "d B Y",
      "Y b d", "Y B d", "b Y", "B Y", "Y b", "Y B", "Y"
    ),
    tz = "UTC",
    exact = FALSE
  ))

  yr <- lubridate::year(parsed)
  fallback <- suppressWarnings(as.integer(str_extract(z, "(19|20)\\d{2}")))
  ifelse(is.na(yr), fallback, yr)
}

classify_primary_source <- function(source_type, isolation_source, host, ifsac, collected_by) {
  st <- str_to_lower(clean_text(source_type))
  iso <- str_to_lower(clean_text(isolation_source))

  iso_norm <- iso |>
    str_replace_all("[^a-z]+", " ") |>
    str_replace_all("chicken(?=[a-z])", "chicken ") |>
    str_replace_all("\\bchicekn\\b|\\bchikcen\\b|\\bchcken\\b", "chicken") |>
    str_squish()

  chicken_pat <- paste0(
    "\\bchicken\\b|\\bchick\\b|\\bpoultry\\b|\\bhen\\b|",
    "\\bbroiler\\b|\\blayer\\b|\\bgallus\\b|gallus\\s+gallus|",
    "\\begg\\b|\\beggs\\b|\\begg\\s*product\\b|\\bproduct\\s*eggs\\b|",
    "\\begg\\s*yolk(s)?\\b|\\byolk(s)?\\b|",
    "\\bhatchery\\b|\\bhatching\\b|\\bhatcher(y)?\\b|",
    "\\bbox\\s*liner(s)?\\b|\\bliner(s)?\\b|\\bdrag\\s*swab(s)?\\b"
  )

  src <- case_when(
    str_detect(st, "\\bhuman\\b") ~ "Human",
    str_detect(iso_norm, chicken_pat) ~ "Chicken",
    TRUE ~ "Other"
  )

  host_l <- str_to_lower(clean_text(host))
  ifsac_l <- str_to_lower(clean_text(ifsac))
  collect_l <- str_to_lower(clean_text(collected_by))

  case_when(
    src == "Other" &
      (str_detect(host_l, "chicken") | str_detect(ifsac_l, "chicken")) ~ "Chicken",
    src == "Other" & str_detect(collect_l, "cdc") ~ "Human",
    TRUE ~ src
  )
}

derive_origin_and_source <- function(source, strain, isolate_identifiers, food_origin) {
  usa_pat <- "(\\busa\\b|\\bunited\\s*states\\b|\\bu\\s*\\.\\s*s\\s*\\.\\b|\\bus\\b)"
  id_pat <- "(\\bcfsan\\b|\\bfda\\b)"

  strain_l <- str_to_lower(clean_text(strain))
  id_l <- str_to_lower(clean_text(isolate_identifiers))
  food_l <- str_to_lower(clean_text(food_origin))
  is_fda_like <- str_detect(strain_l, "^fda") | str_detect(id_l, id_pat)

  origin <- case_when(
    !is_fda_like ~ "US",
    food_l == "" ~ "US",
    str_detect(food_l, usa_pat) ~ "US",
    TRUE ~ "Non-US"
  )

  source_group <- case_when(
    str_to_lower(source) == "human" ~ "Human",
    str_to_lower(source) == "chicken" ~ paste0("Chicken-", origin),
    str_to_lower(source) == "other" ~ paste0("Other-", origin),
    TRUE ~ paste0(source, "-", origin)
  )

  final_source <- ifelse(
    source_group %in% c("Chicken-Non-US", "Other-Non-US"),
    "International",
    source_group
  )

  list(
    Origin_group = origin,
    Source_group = source_group,
    Source = final_source
  )
}

classify_animal_source <- function(isolation_source) {
  x <- str_to_lower(clean_text(isolation_source))
  case_when(
    str_detect(x, "\\b(chicken|broiler|hen)\\b") ~ "Chicken",
    str_detect(x, "\\bturkey\\b") ~ "Turkey",
    str_detect(x, "\\b(cattle|bovine|beef|steer|calf|veal|bos taurus)\\b") ~ "Cattle",
    str_detect(x, "\\b(sheep|lamb|ovine|goat|caprine)\\b") ~ "Sheep/goat",
    str_detect(x, "\\b(swine|sow|pork|porcine)\\b") ~ "Swine",
    TRUE ~ "Others"
  )
}

assign_amr_pattern <- function(amr) {
  a <- str_to_lower(ifelse(is.na(amr), "", amr))
  case_when(
    str_detect(a, "gyra[_-]?d87y") ~ "gyrA_D87Y",
    str_detect(a, "gyra[_-]?d87n") ~ "gyrA_D87N",
    str_detect(a, "gyra[_-]?s83y") ~ "gyrA_S83Y",
    str_detect(a, "gyra[_-]?d87g") ~ "gyrA_D87G",
    str_detect(a, "qnrb19") ~ "qnrB19",
    TRUE ~ "others"
  )
}
