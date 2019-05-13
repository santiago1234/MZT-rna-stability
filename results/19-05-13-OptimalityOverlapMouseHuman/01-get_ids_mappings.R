library(biomaRt)
library(tidyverse)


# description -------------------------------------------------------------
# obtain mappings from ensenmble gene ids to Gene names

# get the universe of gene ids --------------------------------------------

ids <- read_csv("../../data/19-04-29-AllDecayProfiles/decay_profiles_and_seqdata.csv")

humangeneids <- 
  ids %>% 
  filter(specie == "human", datatype != "orf-ome") %>% 
  pull(gene_id) %>% 
  unique() %>% 
  tibble(gene_id = .)


# get mappings gene name -> ensembl gene id ->  ---------------------------

human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


res <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters="ensembl_gene_id",
  values=humangeneids$gene_id, mart=human_mart
) %>% 
  as_tibble() %>% 
  rename(
    Name = external_gene_name,
    gene_id = ensembl_gene_id
  )


mappings <- left_join(humangeneids, res)

write_csv(mappings, "results_data/mappings_ensemble_id_to_gene_names.csv")
