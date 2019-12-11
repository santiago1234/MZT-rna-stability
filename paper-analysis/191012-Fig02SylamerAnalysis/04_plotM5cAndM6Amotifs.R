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
  filter(k==6, ranked_by == "residual")

# discard xenopus

sylamer_results <- 
  sylamer_results %>% 
  filter(specie == "fish")


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

sylamer_results %>% 
  ggplot(aes(x=rank, y=log10pavl, group=kmer)) +
  geom_line(position="jitter", color="grey", size=1/10) +
  geom_rangeframe(size=1/3) +
  geom_hline(yintercept = c(-1, 1) * log10(0.001), linetype=3, size=1/3) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 2000, 4000), labels = c("0k", "2k", "4k")) +
  geom_line(
    data = elementos,
    aes(x=rank, y=log10pavl, group=kmer, color=motif, size=kmer), size=1
  ) +
  ggrepel::geom_text_repel(data = elementos_peak, aes(x=rank, y=log10pavl, label=kmer)) +
  scale_color_manual(values = c("#000000", "#E69F00", "#009E73")) +
  labs(
    x = "3' UTRs sorted by residual values",
    y = "log10 (enrichment P-value)"
  )

ggsave("figures/mir430-m6a-m5c-kmerslike.pdf", height = 3, width = 6)

