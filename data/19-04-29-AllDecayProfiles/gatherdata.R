# For data-sets with half-life I convert the half life to
# decay rate
library(tidyverse)
library(biomaRt)

# helper functions --------------------------------------------------------

# embo data ---------------------------------------------------------------

embo2016_fish <- readxl::read_excel("fishxenopus2016/inline-supplementary-material-2.xls", sheet = 1) %>% 
  dplyr::select(gene, Half_life) %>% 
  mutate(
    decay_rate = log(1/2) / Half_life,
    specie = "fish",
    datatype = "aamanitin ribo",
    cell_type = "embryo mzt"
) %>% 
  rename(gene_id = gene) %>% 
  dplyr::select(-Half_life)

embo2016_xenopus <- readxl::read_excel("fishxenopus2016/inline-supplementary-material-2.xls", sheet = 2) %>% 
  dplyr::select(gene, Half_life) %>% 
  mutate(decay_rate = log(1/2) / Half_life) %>% 
  mutate(
    decay_rate = log(1/2) / Half_life,
    specie = "xenopus",
    datatype = "aamanitin ribo",
    cell_type = "embryo mzt"
  ) %>% 
  rename(gene_id = gene) %>% 
  dplyr::select(-Half_life)

# slam-mouse --------------------------------------------------------------

slammouse <- readxl::read_excel("mouseStemCellsSLAMseqNatureBiotechnology/nmeth.4435-S4.xls") %>% 
  mutate(decay_rate = log(1/2) / `Half-life (h)`)

## map ids to ensenmbl ids
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

# in case the same gene-id I will take the median

slammouse <- 
  inner_join(res, slammouse) %>% 
  dplyr::select(gene_id,  decay_rate) %>% 
  group_by(gene_id) %>% 
  summarise(decay_rate = median(decay_rate)) %>% 
  mutate(
    datatype = "slam-seq",
    cell_type = "mES cells",
    specie = "mouse"
  )

# Wu et al 2019 -----------------------------------------------------------
# decay profiles

humanprofiles <- read_csv("../../../180815-orfome/data/human_profiles_decay_rate.csv") %>% 
  dplyr::select(Name, datatype, cell_type, decay_rate) %>% 
  mutate(specie = "human") %>% 
  rename(gene_id = Name)

## add RPE and HELA data

hela_decay <- readxl::read_xlsx("Wu et al_decay_rates.xlsx", sheet = 2)
rpe_decay <-  readxl::read_xlsx("Wu et al_decay_rates.xlsx", sheet = 3)

# NOTE: decay_rate has the opposite sign in these 2 cell types
rpe_and_hela <- bind_rows(
  hela_decay,
  rpe_decay
) %>% 
  mutate(specie = "human", decay_rate = decay_rate * -1) %>% 
  rename(gene_id = Name)

## add data to the other profiles
humanprofiles <- bind_rows(humanprofiles, rpe_and_hela)

# fish Gopal profiles -----------------------------------------------------

fish_new <- read_csv("../../results/19-01-11-GetDecayRateFromTimeCourse/results_data/estimated_decay_rates.csv") %>% 
  filter(term == "beta") %>% 
  dplyr::select(Gene_ID, estimate) %>% 
  mutate(
    specie = "fish",
    datatype = "aamanitin polya",
    cell_type = "embryo mzt"
  ) %>% 
  rename(
    decay_rate = estimate,
    gene_id = Gene_ID
  )

# merge the data sets -----------------------------------------------------

decay_profiles <- 
  bind_rows(
  humanprofiles,
  fish_new,
  embo2016_fish,
  embo2016_xenopus,
  slammouse
)

# - -----------------------------------------------------------------------
# get open reading frames -------------------------------------------------
# - -----------------------------------------------------------------------

get_seq_data <- function(dataset_specie, gene_id_filters) {
  
  # retrive 3' UTR amd orf seqs
  ensemble <- useMart("ensembl", dataset_specie)
  
  seqs_coding <- getBM(
    attributes = c("ensembl_gene_id", "coding"),
    filters="ensembl_gene_id",
    values=gene_id_filters,
    mart=ensemble
  ) %>% 
  as_tibble()
  
  seqs_3utr <- getBM(
    attributes = c("ensembl_gene_id", "3utr"),
    filters="ensembl_gene_id",
    values=gene_id_filters,
    mart=ensemble
  ) %>% 
    as_tibble()
  
  # apply filters
  seqs_coding <- 
    seqs_coding %>% 
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
  
  seqs_3utr <- 
    seqs_3utr %>% 
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
  
  # left join not all orfs contain 3utrs at least for xenopus
  left_join(seqs_coding, seqs_3utr)
}

# get the ORFs and 3UTR seqs
mouse_orfs <- get_seq_data("mmusculus_gene_ensembl",
                           slammouse$gene_id)
xenopus_orfs <- get_seq_data("xtropicalis_gene_ensembl",
                             embo2016_xenopus$gene_id)
fish_orfs <- get_seq_data("drerio_gene_ensembl",
                          unique(c(embo2016_fish$gene_id, fish_new$gene_id)))
human_orfs <- get_seq_data("hsapiens_gene_ensembl",
                           unique(filter(humanprofiles, datatype != 'orf-ome')$gene_id))
# get orf cds sequences
orf_orfs <- readxl::read_excel("orf-data/TRC3_puro_customer informatics lot 11241401MN.xlsx") %>% 
  dplyr::select(DNA_Barcode, ORF_Sequence, `3_Flank`) %>% 
  filter(
    DNA_Barcode != "n/a",
    (str_length(ORF_Sequence) %% 3) == 0
    ) %>% 
  rename(ensembl_gene_id = DNA_Barcode, coding = ORF_Sequence, `3utr` = `3_Flank`)

# bind data ---------------------------------------------------------------

orfs <- bind_rows(
  mouse_orfs,
  xenopus_orfs,
  fish_orfs,
  human_orfs,
  orf_orfs
) %>% 
  rename(gene_id = ensembl_gene_id)

data <- inner_join(decay_profiles, orfs)

write_csv(data, "decay_profiles_and_seqdata.csv")
