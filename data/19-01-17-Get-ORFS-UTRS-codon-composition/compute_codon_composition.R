library(tidyverse)
library(codonr) # this is my R package

data_fish <- read_csv("sequence-data/fish_seq_data_cds_3utr.csv")

fish_cc <- codon_composition(data_fish, orf_col = "coding", id_col = "ensembl_gene_id")

write_csv(fish_cc, "sequence-data/zfish_codon_composition.csv")
