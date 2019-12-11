library(tidyverse)
library(ggthemes)

theme_set(theme_tufte(base_family = 'Helvetica'))


sylamer_results <- list.files("results_data", pattern = "sylamerResult*", full.names = T)

load_syla_file <- function(syl_file) {
  id_name <- basename(syl_file) %>% 
    str_remove("sylamerResult_") %>% 
    str_remove(".tsv") %>% 
    str_remove("_k\\d")
  
  read_tsv(syl_file) %>% 
    mutate(utrs_sorted_from = id_name, `0` = 0) %>% 
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

# define the miR-430 kmers ------------------------------------------------

miR430_kmer <- tibble(
  kmer = c(
    "GCACTT",
    "AGCACT",
    "GCACTTA",
    "AGCACTT",
    "AGCACTTA"
  )
) %>% 
  mutate(
    kmer_name = str_c(str_length(kmer), "-mer ", kmer)
  )
miR430_kmer

# collect the sylamer results for miR-430 ---------------------------------

miR430_kmer <- inner_join(miR430_kmer, sylamer_results)


miR430_kmer %>% 
  ggplot(aes(x=rank, y=log10pavl, color=ranked_by)) +
  geom_hline(yintercept = 0, size=1/5) +
  geom_hline(yintercept = -log10(0.1 / 5) *c(-1, 1), linetype=2, size=1/5) +
  geom_line(size=1) +
  geom_rangeframe(color="black", size=1/5) +
  geom_line(aes(group=ranked_by), size=1/5, color="black") +
  scale_color_viridis_d(option = "B") +
  scale_x_continuous(breaks = c(0, 2000, 4000), labels = c("0k", "2k", "4k") ) +
  facet_grid(specie~kmer_name, scales = "free_x") +
  labs(
    y = "log10 (enrichment P value)",
    x = "sorted 3' UTRs",
    title = "Sylamer enrichment landscape for miR-430 k-mers",
    color = "sorted from:"
  )

ggsave("./figures/SylamerLandscapePlorMir430.pdf", height = 3, width = 7)

