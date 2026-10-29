#!/usr/bin/env Rscript
source(file.path(Sys.getenv("SEN_ROOT", unset = getwd()), "analysis/part_a/lib/common.R"))

if (!requireNamespace("patchwork", quietly=TRUE)) {
  stop("Missing R package patchwork in sen_stats environment.")
}
suppressPackageStartupMessages(library(patchwork))

input <- sen_path("results/part_a/master_33831/SEN_33831_DERIVED_ANALYSIS_READY.csv")
if (!file.exists(input)) stop("Run analysis/part_a/06_master_33831_analysis.R first.")
outdir <- ensure_dir(sen_path("results/part_a/master_33831/diversity_full"))

sen <- read_csv(input, show_col_types=FALSE)
assert_columns(sen,c("Run","AMR genotypes core","Source"),"derived 33,831-genome master")
stopifnot(nrow(sen) == 33831L, n_distinct(sen$Run) == 33831L)
# Source diversity uses the entire master, independent of Collection_year.
source_levels <- c("Human","Chicken-US","Other-US","International")

tok <- sen %>%
  mutate(
    row_id=row_number(),
    Source=factor(as.character(Source),levels=source_levels),
    token=extract_amr_tokens(.data[["AMR genotypes core"]])
  ) %>%
  select(row_id,Source,token) %>%
  unnest_longer(token) %>%
  mutate(
    token_norm=str_to_lower(str_replace_all(str_trim(token),"\\s+","")),
    token_pretty=pretty_token(token_norm)
  ) %>%
  filter(str_detect(token_norm,fqr_genotype_token_regex))

top5 <- tok %>%
  count(token_pretty,sort=TRUE) %>%
  slice_head(n=5) %>%
  pull(token_pretty)

marker_levels <- c(top5,"others")

marker_long <- tok %>%
  mutate(marker_group=if_else(token_pretty %in% top5,token_pretty,"others")) %>%
  distinct(row_id,Source,marker_group)

iso_mat <- marker_long %>%
  mutate(value=1L) %>%
  pivot_wider(
    id_cols=c(row_id,Source),
    names_from=marker_group,
    values_from=value,
    values_fill=0L
  )

for (m in setdiff(marker_levels,names(iso_mat))) iso_mat[[m]] <- 0L
iso_mat <- iso_mat %>% select(row_id,Source,all_of(marker_levels))

hill_numbers <- function(p) {
  p <- as.numeric(p)
  p <- p[is.finite(p) & p > 0]
  if (!length(p)) return(c(q0=0,q1=0,q2=0))
  c(
    q0=length(p),
    q1=exp(-sum(p*log(p))),
    q2=1/sum(p^2)
  )
}

calc_div <- function(df) {
  df %>%
    group_by(Source) %>%
    summarise(across(all_of(marker_levels),mean),.groups="drop") %>%
    rowwise() %>%
    mutate(
      s=sum(c_across(all_of(marker_levels))),
      p=list(if(s>0)c_across(all_of(marker_levels))/s else rep(0,length(marker_levels))),
      q0=hill_numbers(unlist(p))["q0"],
      q1=hill_numbers(unlist(p))["q1"],
      q2=hill_numbers(unlist(p))["q2"]
    ) %>%
    ungroup() %>%
    select(Source,q0,q1,q2)
}

cat("[AUDIT] QR/FQR marker-positive genomes:", n_distinct(iso_mat$row_id), "\n")
print(iso_mat %>% count(Source, name="QR_positive_N"))
point <- calc_div(iso_mat) %>%
  left_join(iso_mat %>% count(Source,name="N"),by="Source")

set.seed(1)
B <- 1000L

boot_one <- function(df) {
  X <- as.matrix(df[,marker_levels,drop=FALSE])
  n <- nrow(X)
  out <- matrix(NA_real_,B,3,dimnames=list(NULL,c("q0","q1","q2")))
  for (b in seq_len(B)) {
    prev <- colMeans(X[sample.int(n,n,replace=TRUE),,drop=FALSE])
    p <- if(sum(prev)>0) prev/sum(prev) else rep(0,length(prev))
    out[b,] <- hill_numbers(p)
  }
  as.data.frame(out)
}

boot <- split(iso_mat,iso_mat$Source) %>%
  lapply(boot_one) %>%
  Map(function(x,nm){x$Source<-nm;x},.,names(.)) %>%
  bind_rows()

ci <- boot %>%
  group_by(Source) %>%
  summarise(
    q0_lo=quantile(q0,.025),q0_hi=quantile(q0,.975),
    q1_lo=quantile(q1,.025),q1_hi=quantile(q1,.975),
    q2_lo=quantile(q2,.025),q2_hi=quantile(q2,.975),
    .groups="drop"
  )

div_sum <- point %>% left_join(ci,by="Source")
write_csv(
  div_sum,
  file.path(outdir,"QR_diversity_q0_q1_q2_bySource_with95CI.csv")
)

comp <- iso_mat %>%
  group_by(Source) %>%
  summarise(across(all_of(marker_levels),mean),.groups="drop") %>%
  pivot_longer(all_of(marker_levels),names_to="marker_group",values_to="prevalence") %>%
  group_by(Source) %>%
  mutate(p = if (sum(prevalence) > 0) prevalence / sum(prevalence) else rep(0, n())) %>%
  ungroup()

write_csv(comp,file.path(outdir,"QR_composition_bySource.csv"))

global_stat <- function(v) max(v,na.rm=TRUE)-min(v,na.rm=TRUE)

perm_global <- function(metric,nperm=2000L,seed=1L) {
  set.seed(seed)
  obs <- calc_div(iso_mat)
  obs_stat <- global_stat(obs[[metric]])
  perm_stats <- replicate(
    nperm,
    {
      x <- iso_mat
      x$Source <- sample(x$Source)
      global_stat(calc_div(x)[[metric]])
    }
  )
  tibble(
    metric=metric,
    obs_stat=obs_stat,
    p_value=mean(perm_stats>=obs_stat)
  )
}

perm_pair <- function(metric,nperm=2000L,seed=1L) {
  set.seed(seed)
  pairs <- combn(source_levels,2,simplify=FALSE)
  out <- bind_rows(lapply(pairs,function(pp){
    x <- iso_mat %>% filter(as.character(Source)%in%pp)
    obs <- calc_div(x)
    obs_diff <- abs(
      obs[[metric]][as.character(obs$Source)==pp[1]] -
      obs[[metric]][as.character(obs$Source)==pp[2]]
    )
    perm_diff <- replicate(
      nperm,
      {
        z <- x
        z$Source <- sample(z$Source)
        q <- calc_div(z)
        abs(
          q[[metric]][as.character(q$Source)==pp[1]] -
          q[[metric]][as.character(q$Source)==pp[2]]
        )
      }
    )
    tibble(
      metric=metric,
      group1=pp[1],
      group2=pp[2],
      obs_diff=obs_diff,
      p_value=mean(perm_diff>=obs_diff)
    )
  }))
  out %>%
    mutate(p_adj_BH=p.adjust(p_value,"BH")) %>%
    arrange(p_adj_BH)
}

write_csv(perm_global("q1",2000,101),file.path(outdir,"QR_diversity_global_perm_q1.csv"))
write_csv(perm_global("q2",2000,102),file.path(outdir,"QR_diversity_global_perm_q2.csv"))
write_csv(perm_pair("q1",2000,201),file.path(outdir,"QR_diversity_pairwise_perm_q1.csv"))
write_csv(perm_pair("q2",2000,202),file.path(outdir,"QR_diversity_pairwise_perm_q2.csv"))

plot_div <- div_sum %>%
  pivot_longer(c(q1,q2),names_to="metric",values_to="value") %>%
  mutate(
    lo=if_else(metric=="q1",q1_lo,q2_lo),
    hi=if_else(metric=="q1",q1_hi,q2_hi),
    metric=recode(metric,q1="Exponential Shannon (q1)",q2="Inverse Simpson (q2)")
  )

p_div <- ggplot(plot_div,aes(Source,value)) +
  geom_point(size=2.8) +
  geom_errorbar(aes(ymin=lo,ymax=hi),width=.15) +
  facet_wrap(~metric,nrow=1) +
  labs(x="Source",y="Effective number of QR/FQR categories") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45,hjust=1))

p_comp <- ggplot(comp,aes(Source,p,fill=marker_group)) +
  geom_col() +
  scale_y_continuous(labels=scales::percent) +
  labs(x="Source",y="Normalized QR/FQR composition",fill="QR category") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45,hjust=1))

combo <- p_div | p_comp
ggsave(
  file.path(outdir,"QR_diversity_and_composition.png"),
  combo,width=15,height=6,dpi=600
)

cat("[PASS] Full-master QR/FQR diversity/composition outputs written to",outdir,"\n")
