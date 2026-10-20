#!/usr/bin/env Rscript
source(file.path(Sys.getenv("SEN_ROOT", unset = getwd()), "analysis/part_a/lib/common.R"))

input <- sen_path("data/part_a/phenotype/FDASEN_AMRPhenotype2.csv")
outdir <- ensure_dir(sen_path("results/part_a/phenotype"))

dat <- read_csv(input, show_col_types = FALSE)
assert_columns(
  dat,
  c("Year","Cip_Resistant","Cip_Tested","DSC_Resistant","DSC_Tested","Source"),
  "phenotype input"
)

dat <- dat %>%
  mutate(across(c(Year,Cip_Resistant,Cip_Tested,DSC_Resistant,DSC_Tested), as.numeric))

both_tested <- !is.na(dat$Cip_Tested) & !is.na(dat$DSC_Tested)
if (any(dat$Cip_Tested[both_tested] != dat$DSC_Tested[both_tested])) {
  stop("Cip_Tested and DSC_Tested differ. Historical shared-denominator analysis cannot be reproduced safely.")
}

trend <- dat %>%
  mutate(
    Year = as.integer(Year),
    Source = factor(as.character(Source), levels = c("Human","Chicken"))
  ) %>%
  filter(!is.na(Year), !Year %in% c(2024, 2025)) %>%
  mutate(
    Tested = pmax(Cip_Tested, DSC_Tested, na.rm = TRUE),
    combined_n = Cip_Resistant + DSC_Resistant
  )

if (any(trend$combined_n > trend$Tested, na.rm = TRUE)) {
  stop("Combined Cip + DSC resistant count exceeds tested denominator.")
}

trend <- trend %>%
  mutate(
    CipR_pct = if_else(Tested > 0, 100 * Cip_Resistant / Tested, NA_real_),
    DSC_pct = if_else(Tested > 0, 100 * DSC_Resistant / Tested, NA_real_),
    QRFQR_pct = if_else(Tested > 0, 100 * combined_n / Tested, NA_real_),
    CipR_se = if_else(Tested > 0, 100 * sqrt((Cip_Resistant/Tested)*(1-Cip_Resistant/Tested)/Tested), NA_real_),
    DSC_se = if_else(Tested > 0, 100 * sqrt((DSC_Resistant/Tested)*(1-DSC_Resistant/Tested)/Tested), NA_real_),
    QRFQR_se = if_else(Tested > 0, 100 * sqrt((combined_n/Tested)*(1-combined_n/Tested)/Tested), NA_real_),
    CipR_lo = pmax(0, CipR_pct - 1.96*CipR_se),
    CipR_hi = CipR_pct + 1.96*CipR_se,
    DSC_lo = pmax(0, DSC_pct - 1.96*DSC_se),
    DSC_hi = DSC_pct + 1.96*DSC_se,
    QRFQR_lo = pmax(0, QRFQR_pct - 1.96*QRFQR_se),
    QRFQR_hi = QRFQR_pct + 1.96*QRFQR_se
  ) %>%
  select(
    Year, Source,
    QRFQR_pct,QRFQR_lo,QRFQR_hi,
    DSC_pct,DSC_lo,DSC_hi,
    CipR_pct,CipR_lo,CipR_hi
  ) %>%
  pivot_longer(
    cols = -c(Year,Source),
    names_to = c("Metric",".value"),
    names_pattern = "^(QRFQR|DSC|CipR)_(pct|lo|hi)$"
  ) %>%
  mutate(
    Metric = recode(Metric, QRFQR="Nal-R/DSC or Cip-R", DSC="Nal-R/DSC", CipR="Cip-R"),
    Metric = factor(Metric, levels=c("Nal-R/DSC or Cip-R","Nal-R/DSC","Cip-R"))
  )

write_csv(trend, file.path(outdir, "FQ_trends_Human_vs_Chicken.csv"))

p <- ggplot(trend, aes(Year, pct, color=Source, fill=Source, linetype=Metric)) +
  geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.15, color=NA, na.rm=TRUE) +
  geom_line(linewidth=0.9, na.rm=TRUE) +
  geom_point(aes(shape=Metric), size=2, na.rm=TRUE) +
  scale_shape_manual(values=c("Nal-R/DSC or Cip-R"=16,"Nal-R/DSC"=15,"Cip-R"=17)) +
  scale_color_manual(values=c("Human"="red","Chicken"="blue")) +
  scale_fill_manual(values=c("Human"="red","Chicken"="blue")) +
  scale_y_continuous(limits=c(0,NA), expand=expansion(mult=c(0,0.05))) +
  labs(
    x="Year of Isolation", y="Percentage (%)",
    color="Source", fill="Source", linetype="Phenotype", shape="Phenotype"
  ) +
  theme_bw() +
  theme(
    axis.title=element_text(face="bold"),
    axis.text.x=element_text(angle=90,vjust=.5,hjust=1)
  )

ggsave(file.path(outdir,"FQ_trends_Human_vs_Chicken.png"),p,width=9,height=4.8,dpi=1200)
ggsave(file.path(outdir,"FQ_trends_Human_vs_Chicken.pdf"),p,width=9,height=4.8)
ggsave(
  file.path(outdir,"FQ_trends_Human_vs_Chicken.tif"),
  p,width=9,height=4.8,dpi=600,compression="lzw"
)

cat("[PASS] Part A1 phenotype trends written to", outdir, "\n")
