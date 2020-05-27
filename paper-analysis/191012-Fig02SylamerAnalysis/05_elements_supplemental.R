library(tidyverse)
library(ggthemes)

theme_set(theme_tufte(base_family = "Helvetica"))


# load syla files ---------------------------------------------------------

sylamer_results <- list.files("results_data", pattern = "sylamerResult*", full.names = T)

load_syla_file <- function(syl_file) {
  id_name <- basename(syl_file) %>% 
    str_remove("sylamerResult_") %>% 
    str_remove(".tsv") %>% 
    str_remove("_k\\d")
  
  read_tsv(syl_file) %>% 
    mutate(
      utrs_sorted_from = id_name,
      `0` = 0 # sylamer plot should start from 0
    ) %>% 
    gather(key = rank, value = log10pavl, -upper, -utrs_sorted_from) %>% 
    mutate(rank = as.integer(rank)) %>% 
    dplyr::rename(kmer = upper) %>% 
    mutate(
      pvalue = 10**-(abs(log10pavl)),
      rank = as.integer(rank),
      k = str_length(kmer)
    ) %>% 
    separate(utrs_sorted_from, into = c("specie", "ranked_by"), sep = "_")
  
}

sylamer_results <- 
  map_df(sylamer_results, load_syla_file)

sylamer_results <- 
  sylamer_results %>% 
  filter(k==6)

# discard xenopus

sylamer_results <- 
  sylamer_results %>% 
  filter(specie == "fish") %>% 
  mutate(
    ranked_by = factor(ranked_by, levels = c("residual", "observedlog2FC", "predictedStability"))
  )


# kier elements, i go this by trial and error visualiz --------------------

elementos <- bind_rows(
  tibble(kmer=c("GGACTT", "TAGGAC"), motif = "m6a"),
  tibble(kmer=c("TCTATC", "CTATCT", "TATCTA"), motif = "m5c"),
  tibble(kmer=c("GCACTT"), motif = "miR-430")
) %>% 
  inner_join(sylamer_results)

# draw the name

elementos_peak <- elementos %>% 
  group_by(kmer) %>% 
  filter(log10pavl == max(log10pavl))

set.seed(42)
## sample_300 k-mers
mers_sample <- sylamer_results$kmer %>% 
  unique() %>% 
  sample(1000)

sylamer_results %>% 
  filter(!kmer %in% elementos$kmer) %>% 
  filter(kmer %in% mers_sample) %>% 
  ggplot(aes(x=rank, y=log10pavl, group=kmer)) +
  geom_line(position="jitter", color="grey", size=1/10) +
  geom_rangeframe(size=1/3) +
  geom_hline(yintercept = c(-1, 1) * log10(0.001), linetype=3, size=1/3) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 2000, 4000), labels = c("0k", "2k", "4k")) +
  scale_y_continuous(breaks = c(-3, 0, 3, 5, 10)) +
  geom_line(
    data = elementos,
    aes(x=rank, y=log10pavl, group=kmer, color=motif, size=kmer), size=1
  ) +
  ggrepel::geom_text_repel(data = elementos_peak, aes(x=rank, y=log10pavl, label=kmer)) +
  scale_color_manual(values = c("#000000", "#E69F00", "#009E73")) +
  labs(
    x = "3' UTRs sorted by residual values",
    y = "log10 (enrichment P-value)"
  ) +
  facet_grid(~ranked_by)


ggsave("figures/mir430-m6a-m5c-kmerslike-Supplemental.pdf", height = 4, width = 11)
best_pvalues <- elementos %>% 
  group_by(motif, ranked_by) %>% 
  filter(log10pavl == max(log10pavl))

best_pvalues %>% 
  write_csv("results_data/pvalues_comparing_residual_VS_WT.csv")

best_pvalues %>% 
  ungroup() %>% 
  mutate(ranked_by = fct_rev(ranked_by)) %>% 
  ggplot(aes(x=ranked_by, y=log10pavl, group=ranked_by, fill=motif)) +
  geom_bar(stat = "identity", position = 'dodge', color=NA, size=1/4) +
  scale_fill_manual(values = c("#000000", "#E69F00", "#009E73")) +
  coord_flip() +
  theme(legend.position = 'none') +
  facet_grid(motif~.)
ggsave("figures/p-vals-plot.pdf", height = 2, width = 3)
