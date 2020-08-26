library(tidyverse)
library(optimalcodonR)

d <- read_csv("data/position_seqs.csv") %>% 
  filter(nchar(coding) > 9)


# add the codon counts for each position ----------------------------------

d_first <- 
  d %>% 
  select(gene_id, first) %>% 
  rename(coding = first) %>% 
  add_codon_counts() %>% 
  rename_with(function(x) paste0(x, "_first"), starts_with("c_")) %>% 
  select(-coding)

d_second <- 
  d %>% 
  select(gene_id, second) %>% 
  rename(coding = second) %>% 
  add_codon_counts() %>% 
  rename_with(function(x) paste0(x, "_second"), starts_with("c_")) %>% 
  select(-coding)

d_third <- 
  d %>% 
  select(gene_id, third) %>% 
  rename(coding = third) %>% 
  add_codon_counts() %>% 
  rename_with(function(x) paste0(x, "_third"), starts_with("c_")) %>% 
  select(-coding)

# all regions no positions
# the coding colum is the whole cds
d_all <- 
  d %>% 
  select(gene_id, coding) %>% 
  add_codon_counts() %>% 
  rename_with(function(x) paste0(x, "_all"), starts_with("c_")) %>% 
  select(-coding)


# merge the counts in the main table --------------------------------------


d <- d %>% 
  select(-c(coding, first, second, third)) %>% 
  inner_join(d_first) %>% 
  inner_join(d_second) %>% 
  inner_join(d_third) %>% 
  inner_join(d_all)

write_csv(d, "data/counts_by_position.csv")
