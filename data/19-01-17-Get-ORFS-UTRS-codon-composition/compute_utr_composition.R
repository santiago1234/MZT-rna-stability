library(tidyverse)
library(codonr)

utrs <- read_csv("sequence-data/fish_seq_data_cds_3utr.csv") %>% 
  filter(str_length(`3utr`) > 10)

utr_counts <- 
  utrs %>% 
  mutate(kmer = map(`3utr`, ~count_kmers(., k = 6))) %>% 
  select(-coding) %>% 
  unnest(kmer) %>% 
  spread(key = kmer, value = n) %>% 
  replace(is.na(.), 0)

write_csv(utr_counts, "sequence-data/zfish_3utr6mer_composition.csv")