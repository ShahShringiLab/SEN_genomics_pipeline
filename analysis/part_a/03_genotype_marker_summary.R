#!/usr/bin/env Rscript
source(file.path(Sys.getenv("SEN_ROOT", unset = getwd()), "analysis/part_a/lib/common.R"))

input <- sen_path("results/part_a/cohort/Finalupdated_SEN_Genomes1_QC_merged.csv")
if (!file.exists(input)) stop("Run python3 scripts/merge_qc_full_metadata.py first.")
outdir <- ensure_dir(sen_path("results/part_a/genotype_summary"))

sen <- read_csv(input, show_col_types=FALSE)
assert_columns(
  sen,
  c("Isolate","AMR genotypes core","QR_genotype","Collection_year","Source","Clade"),
  "QC-merged 3,434-genome cohort"
)

tok <- sen %>%
  transmute(
    isolate_id=as.character(Isolate),
    Source=as.character(Source),
    Collection_year=as.integer(Collection_year),
    token=extract_amr_tokens(.data[["AMR genotypes core"]])
  ) %>%
  unnest_longer(token) %>%
  mutate(
    token_norm=str_to_lower(str_replace_all(str_trim(token),"\\s+",""))
  ) %>%
  filter(str_detect(token_norm, fqr_genotype_token_regex)) %>%
  mutate(token_pretty=pretty_token(token_norm)) %>%
  distinct(isolate_id, Source, Collection_year, token_norm, token_pretty)

N <- n_distinct(sen$Isolate)
positive_ids <- unique(tok$isolate_id)

marker_table <- tok %>%
  distinct(isolate_id, token_pretty) %>%
  count(token_pretty, name="n_isolates") %>%
  mutate(
    Type="QR/FQR marker",
    N=N,
    p=n_isolates/N,
    pct=100*p,
    se=sqrt(p*(1-p)/N),
    lo=pmax(0,100*(p-1.96*se)),
    hi=pmin(100,100*(p+1.96*se))
  )

n_susc <- N-length(positive_ids)
p_susc <- n_susc/N
marker_table <- bind_rows(
  marker_table,
  tibble(
    token_pretty="Q/FQ-susceptible",
    n_isolates=n_susc,
    Type="Q/FQ-susceptible",
    N=N,
    p=p_susc,
    pct=100*p_susc,
    se=sqrt(p_susc*(1-p_susc)/N),
    lo=pmax(0,100*(p_susc-1.96*sqrt(p_susc*(1-p_susc)/N))),
    hi=pmin(100,100*(p_susc+1.96*sqrt(p_susc*(1-p_susc)/N)))
  )
) %>%
  arrange(desc(n_isolates))

write_csv(
  marker_table %>% select(token_pretty,Type,n_isolates,N,pct,lo,hi),
  file.path(outdir,"FQR_marker_frequency_percentage_95CI_with_susceptible.csv")
)

iso_marker <- tok %>%
  distinct(isolate_id, token_pretty) %>%
  group_by(isolate_id) %>%
  summarise(
    n_markers=n(),
    markers=paste(sort(token_pretty),collapse=";"),
    .groups="drop"
  )

write_csv(iso_marker,file.path(outdir,"FQR_marker_burden_per_isolate.csv"))

marker_bins <- iso_marker %>%
  mutate(
    bin = case_when(
      n_markers == 1 ~ "1 marker",
      n_markers == 2 ~ "2 markers",
      n_markers > 2 ~ ">2 markers",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(bin)) %>%
  count(bin, name="n_isolates") %>%
  mutate(
    pct_isolates = 100 * n_isolates / sum(n_isolates),
    bin = factor(bin, levels=c("1 marker","2 markers",">2 markers"))
  ) %>%
  arrange(bin)

write_csv(marker_bins,file.path(outdir,"FQR_marker_burden_bins.csv"))

write_csv(
  iso_marker %>% filter(n_markers>1),
  file.path(outdir,"isolates_more_than_one_FQR_marker.csv")
)

combo_2 <- iso_marker %>%
  filter(n_markers==2) %>%
  count(markers,name="n_isolates") %>%
  mutate(pct=100*n_isolates/sum(n_isolates)) %>%
  arrange(desc(n_isolates),markers)

write_csv(combo_2,file.path(outdir,"FQR_exactly_2_marker_combinations.csv"))

combo_gt2 <- iso_marker %>%
  filter(n_markers>2) %>%
  count(n_markers,markers,name="n_isolates") %>%
  mutate(pct=100*n_isolates/sum(n_isolates)) %>%
  arrange(desc(n_isolates),desc(n_markers),markers)

write_csv(combo_gt2,file.path(outdir,"FQR_greater_than_2_marker_combinations.csv"))

year_stack <- sen %>%
  mutate(
    Collection_year=as.integer(Collection_year),
    QR_status=if_else(as.integer(QR_genotype)==1L,"QR/FQR marker","Q/FQ-susceptible")
  ) %>%
  filter(!is.na(Collection_year),Collection_year>=1990,Collection_year<=2025) %>%
  count(Source,Collection_year,QR_status,name="n") %>%
  group_by(Source,Collection_year) %>%
  mutate(N_year=sum(n),pct=100*n/N_year) %>%
  ungroup()

write_csv(year_stack,file.path(outdir,"QRFQR_by_source_and_year.csv"))

top5 <- tok %>%
  count(token_pretty,sort=TRUE) %>%
  slice_head(n=5) %>%
  pull(token_pretty)

marker_prev <- tok %>%
  mutate(marker_group=if_else(token_pretty %in% top5,token_pretty,"others")) %>%
  distinct(isolate_id,Source,Collection_year,marker_group) %>%
  count(Source,Collection_year,marker_group,name="n") %>%
  left_join(
    sen %>%
      transmute(Source=as.character(Source),Collection_year=as.integer(Collection_year)) %>%
      filter(!is.na(Collection_year)) %>%
      count(Source,Collection_year,name="N"),
    by=c("Source","Collection_year")
  ) %>%
  mutate(pct=100*n/N)

write_csv(
  marker_prev,
  file.path(outdir,"top5_QRFQR_marker_prevalence_by_source_year.csv")
)

plot_df <- marker_table %>%
  slice_head(n=min(100,nrow(marker_table))) %>%
  mutate(
    token_pretty=factor(token_pretty,levels=rev(token_pretty)),
    pct_label=if_else(pct<0.01,"<0.01%",sprintf("%.3f%%",pct)),
    label=paste0(n_isolates," (",pct_label,")")
  )

p_tokens <- ggplot(plot_df,aes(token_pretty,n_isolates,fill=Type)) +
  geom_col() +
  geom_text(aes(label=label),hjust=-0.1,size=3) +
  coord_flip() +
  scale_fill_manual(values=c("Q/FQ-susceptible"="green4","QR/FQR marker"="red3")) +
  scale_y_continuous(expand=expansion(mult=c(0,0.18))) +
  labs(x="QR/FQR marker",y="Number of genomes (n (%))") +
  theme_bw() +
  theme(legend.position="none")

ggsave(
  file.path(outdir,"FQR_marker_counts_top100.png"),
  p_tokens,width=7.5,height=10,dpi=600
)

p_year <- ggplot(year_stack,aes(factor(Collection_year),pct,fill=QR_status)) +
  geom_col() +
  facet_wrap(~Source,ncol=1) +
  scale_fill_manual(values=c("Q/FQ-susceptible"="green4","QR/FQR marker"="red3")) +
  labs(x="Collection year",y="Genomes (%)",fill="") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,vjust=.5,hjust=1,size=7))

ggsave(
  file.path(outdir,"QRFQR_percent_by_source_year.png"),
  p_year,width=8.1,height=6,dpi=600
)

cat("[PASS] Part A2 marker summaries written to",outdir,"\n")
