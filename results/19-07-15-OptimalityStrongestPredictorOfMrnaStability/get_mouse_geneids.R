library(tidyverse)
library(biomaRt)

# get mappings name -> gene_ids mouse -------------------------------------

slammouse <- readxl::read_excel("../../data/19-04-29-AllDecayProfiles/mouseStemCellsSLAMseqNatureBiotechnology/nmeth.4435-S4.xls") %>% 
  mutate(decay_rate = log(1/2) / `Half-life (h)`)

mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
res <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters="external_gene_name",
  values=slammouse$Name, mart=mouse
) %>% 
  as_tibble() %>% 
  rename(
    Name = external_gene_name,
    gene_id = ensembl_gene_id
  )

write_csv(res, "results_data/mouse_names_to_ensemblGeneIDs.csv")

# I donwloaded mapping for human to transcripts directly from biomart: mart_export.txt