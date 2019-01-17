library(tidyverse)
library(biomaRt)
library(magrittr)


# download ORFs sequences for fish ----------------------------------------

## I will use a filter from the genes that we quantified

filter_genes <- 
  "../19-01-09-RNAseqProfilesFish/rna-seq-profiles/alpha-amanitin-prolife.csv" %>% 
  read.csv(stringsAsFactors = F) %>% 
  .$Gene_ID %>% 
  unique()

ensembl <- useMart("ensembl", dataset = "drerio_gene_ensembl")


codingseq <- getBM(
  attributes = c("ensembl_gene_id", "coding"),
  filters = "ensembl_gene_id",
  values = filter_genes,
  mart = ensembl
)

utr3_seq <- getBM(
  attributes = c("ensembl_gene_id", "3utr"),
  filters = "ensembl_gene_id",
  values = filter_genes,
  mart = ensembl
)

# aply filters to the sequences -------------------------------------------

codingseq <- 
  codingseq %>% 
  as_tibble() %>% 
  filter(
    !str_detect(coding, "Sequence unavailable"),
    str_length(coding) %% 3 == 0, # coding should be a multiple of 3
    !str_detect(coding, "N") # No Ns is better
  ) %>% 
  group_by(ensembl_gene_id) %>% 
  do(
    arrange(., -str_length(coding)) %>% # get only 1 isophorm (keep longest)
      slice(1)
  ) %>% 
  ungroup()

utr3_seq <- 
  utr3_seq %>% 
  as_tibble() %>% 
  filter(
    !str_detect(`3utr`, "Sequence unavailable"),
    !str_detect(`3utr`, "N")
  ) %>% 
  group_by(ensembl_gene_id) %>% 
  do(
    arrange(., -str_length(`3utr`)) %>% # get only 1 isophorm (keep longest)
      slice(1)
  ) %>% 
  ungroup()

# join table and saves results --------------------------------------------

inner_join(codingseq, utr3_seq) %>% 
  write_csv("sequence-data/fish_seq_data_cds_3utr.csv")
