#!/usr/bin/env Rscript
# Stage 0 of the original Part A2: normalized 33,831 x 51 raw-data snapshot.
# No derived columns and no changes to data/part_a/genotype/US-SEN_12192025.csv.
source(file.path(Sys.getenv("SEN_ROOT", unset=getwd()), "analysis/part_a/lib/common.R"))
input <- sen_path("data/part_a/genotype/US-SEN_12192025.csv")
out <- ensure_dir(sen_path("results/part_a/master_33831/00_parsed_raw"))
if (!requireNamespace("writexl", quietly=TRUE)) {
  stop("Package writexl missing. Run: conda install -n sen_stats -c conda-forge r-writexl -y")
}
raw <- read_sen_master_csv(input)
stopifnot(nrow(raw)==33831L, ncol(raw)==51L,
          dplyr::n_distinct(raw$Run)==33831L,
          dplyr::n_distinct(raw$Isolate)==33831L)
# Flag suspicious displaced metadata (rather than silently shifting columns).
audit_fields <- function(field, regex) {
  if (!(field %in% names(raw))) return(tibble::tibble(field=field, checked=FALSE, nonmissing=NA_integer_, invalid=NA_integer_))
  x <- trimws(as.character(raw[[field]]))
  valid <- !is.na(x) & x!="" & x!="NA"
  bad <- valid & !grepl(regex,x,perl=TRUE)
  tibble::tibble(field=field, checked=TRUE, nonmissing=sum(valid), invalid=sum(bad))
}
audit <- dplyr::bind_rows(
  audit_fields("BioProject","^PRJ[A-Z]{2}[0-9]+$"),
  audit_fields("BioSample","^SAM[A-Z][0-9]+$"),
  audit_fields("Run","^SRR[0-9]+$")
)
readr::write_csv(audit, file.path(out,"parsed_master_field_alignment_audit.csv"))
if (any(audit$invalid>0,na.rm=TRUE)) {
  print(audit)
  stop("Potentially misaligned accession fields in parsed master. Inspect audit; no snapshot published.")
}
# CSV and XLSX contain the SAME 51 source columns, unchanged after parsing.
csv_path <- file.path(out,"SEN_33831_PARSED_51COL.csv")
xlsx_path <- file.path(out,"SEN_33831_PARSED_51COL.xlsx")
readr::write_csv(raw,csv_path,na="")
writexl::write_xlsx(list(Parsed_master_51col=raw),xlsx_path)
checksums <- tibble::tibble(path=c(input,csv_path,xlsx_path),
                           md5=unname(tools::md5sum(c(input,csv_path,xlsx_path))))
readr::write_csv(checksums,file.path(out,"file_checksums.csv"))
cat("[PASS] Stage 0 parsed master: 33,831 rows x 51 original columns.\n")
cat("[PASS] Field alignment audit and raw XLSX/CSV snapshots: ",out,"\n",sep="")
