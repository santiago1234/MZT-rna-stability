library(tidyverse)
library(biomaRt)

# get the decay rate data -------------------------------------------------

decay_rates_human <- read_csv("../../../180815-orfome/analysis/18-12-08-CompareAllHalfLifeCSC/data/all_decay_rates.csv") %>% 
  filter(datatype != "orf-ome")


# get the UTR’s for human genes -------------------------------------------


ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

filter_genes <- unique(decay_rates_human$Name)

utr3_seq <- getBM(
  attributes = c("ensembl_gene_id", "3utr"),
  filters = "ensembl_gene_id",
  values = filter_genes,
  mart = ensembl
)

## apply some filter to the data

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

# make predictors data ----------------------------------------------------
# here i count the kmers of interest

predictors <- 
  utr3_seq %>% 
  mutate(
    MiR430 = str_detect(`3utr`, "GCACTT"),
    stable_TCGGCG = str_detect(`3utr`, "TCGGCG"),
    stable_CTCCCG = str_detect(`3utr`, "CTCCCG"),
    stable_CCTGGG = str_detect(`3utr`, "CCTGGG"),
    unstable_AGACGG = str_detect(`3utr`, "AGACGG"),
    unstable_TCCGTA = str_detect(`3utr`, "TCCGTA")
  ) %>% 
  dplyr::select(-`3utr`) %>% 
  rename(Name = ensembl_gene_id)


datp <- decay_rates_human %>% 
  dplyr::select(datatype, cell_type, Name, decay_rate) %>% 
  inner_join(predictors) %>% 
  gather(key = kmer, value = target, -datatype, -cell_type, -Name, -decay_rate)

datp %>% 
  filter(datatype == "slam-seq") %>% 
  ggplot(aes(x = decay_rate, color = target)) +
  stat_ecdf() +
  ggthemes::scale_color_colorblind() +
  facet_grid(cell_type~kmer) +
  coord_cartesian(xlim = c(-.3, 0)) +
  labs(
    title = "SLAM-seq data"
  )
ggsave("kmers-effect-slam-data.pdf", width = 9.5, height = 2.3)

datp %>% 
  filter(datatype != "slam-seq") %>% 
  ggplot(aes(x = decay_rate, color = target)) +
  stat_ecdf() +
  ggthemes::scale_color_colorblind() +
  facet_grid(cell_type~kmer)+
  coord_cartesian(xlim = c(-.2, .2)) +
  labs(
    title = "Endogenous data"
  )
ggsave("kmers-effect-endogenous-data.pdf", width = 9.5, height = 4)
