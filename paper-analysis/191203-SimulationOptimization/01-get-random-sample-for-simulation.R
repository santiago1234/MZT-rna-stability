library(tidyverse)
library(codonr)



# get a random sample of 250 human genes ---------------------------------
# I am specially interested in this gene

rpl41 <- 
  codonr::train_set %>% 
  filter(specie == "human") %>% 
  filter(gene_id == "ENSG00000229117") %>% 
  select(gene_id, coding) %>% 
  distinct()

set.seed(123)
sample_genes <- 
  codonr::train_set %>% 
  filter(specie == "human") %>% 
  select(gene_id, coding) %>% 
  distinct() %>% 
  sample_n(200)


## append the rpl41 gene

sample_genes <- bind_rows(rpl41, sample_genes)
write_csv(sample_genes, "results-data/sample_genes.csv")
