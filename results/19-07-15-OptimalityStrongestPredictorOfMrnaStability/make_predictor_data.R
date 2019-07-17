# create a matrix that for each specie contain the pathwas informa --------
# optimality (5 PLS components)
# microRNA targets sites
# 3' UTR length
# GC content
# m6A peaks / targets

library(tidyverse)


# m6a data ----------------------------------------------------------------
# fish, for fish we only consider maternal genes

fish_data <- read_csv("../19-02-24-OverlapPathwaysFig3/results_data/regulatory_pathways_matrix.csv") %>% 
  filter(is_maternal) %>% 
  rename(gene_id = Gene_ID)

fish_m6a <- fish_data %>% 
  select(gene_id, m6A) %>% 
  mutate(m6A = as.numeric(m6A)) %>% 
  mutate(specie = "fish")

human_m6a <- readxl::read_excel("results_data/1-s2.0-S1934590914004512-mmc4.xlsx") %>% 
  group_by(`Ensembl Gene ID`) %>% 
  summarise(m6A = n()) %>% 
  rename(gene_id = `Ensembl Gene ID`) %>% 
  mutate(specie = "human")

mouse_m6a <- readxl::read_excel("results_data/1-s2.0-S1934590914004512-mmc2.xlsx") %>% 
  group_by(`Gene Symbol`) %>% 
  summarise(m6A = n()) %>% 
  rename(Name = `Gene Symbol`)

# map gene names to gene_id, mouse table do not contain ensembl gene ids
idsmaps <- read_csv("results_data/mouse_names_to_ensemblGeneIDs.csv")
mouse_m6a <- inner_join(mouse_m6a, idsmaps) %>% 
  select(-Name) %>% 
  mutate(specie = "mouse")

m6a_data <- bind_rows(mouse_m6a, human_m6a, fish_m6a)

# NOTE: if a gene is not included in the m6a_data table, I assume that this
# gene is not a m6A target, missing values must be replaced with 0s


# microRNAs ---------------------------------------------------------------
# microRNAs targets sites
# fish miR-430

fish_mir <- fish_data %>% 
  select(gene_id, MiR430) %>% 
  mutate(microRNAsites = MiR430) %>% 
  select(-MiR430) %>% 
  mutate(specie = "fish")

mouse_mir <- read_tsv("results_data/TargetScan7.1__miR-291-3p_294-3p_295-3p_302-3p.predicted_targets.txt")

mouse_mir <- 
  mouse_mir %>% 
  select(`Target gene`, `Conserved sites total`) %>% 
  rename(microRNAsites = `Conserved sites total`, Name = `Target gene`) %>% 
  inner_join(idsmaps) %>% 
  select(-Name) %>% 
  mutate(specie = "mouse") %>% 
  filter(microRNAsites != "1*") %>% 
  mutate(microRNAsites = as.numeric(microRNAsites))

idsmapshuman <- read_tsv("results_data/mart_export.txt") %>% 
  select(`Gene stable ID`,  `Transcript stable ID`) %>% 
  rename(gene_id = `Gene stable ID`, `Representative transcript` = `Transcript stable ID`)

human_mir <- read_tsv("results_data/TargetScan7.1__miR-17-3p.predicted_targets.txt") %>% 
  select(`Representative transcript`, `Total sites`) %>% 
  mutate(`Representative transcript` = str_replace(`Representative transcript`, "\\..+$", "")) %>% 
  inner_join(idsmapshuman) %>% 
  select(gene_id, `Total sites`) %>% 
  rename(microRNAsites = `Total sites`) %>% 
  mutate(specie = "human")

mir_data <- bind_rows(human_mir, mouse_mir, fish_mir)

# NOTE: if a gene is not included in the m6a_data table, I assume that this
# gene is not a miR target, missing values must be replaced with 0s

# codon optimality (PLS decomposition) ------------------------------------

plsdecomp <- read_csv("../19-06-19-OptimalityAffectsMrnaHumanLevels/results_data/pls_decomposition_all.csv")

# for mouse and human i will use the slam-seq data

hmpls <- plsdecomp %>% 
  filter(datatype == "slam-seq")

fishpls <- plsdecomp %>% 
  filter(specie == "fish", datatype == "aamanitin polya")

optimality <- bind_rows(hmpls, fishpls) %>% 
  select(gene_id:PLS_9)

pathways <- full_join(optimality, m6a_data, by = c("gene_id", "specie")) %>% 
  full_join(mir_data) %>% 
  replace_na(
    list(m6A = 0, microRNAsites = 0)
  ) %>% 
  select(-datatype)


# get the response data ---------------------------------------------------
# for fish I use the log2FC at 5hrs in polya data for maternal genes
# description stability
# mouse => (slam embryonic stem cells)
# human => (k562 slam)
# fish => (fish polya fold change 6/3 hrs)

fishstability <- read_csv("../../data/19-02-05-FoldChangeData/data/log2FC_earlyVSlate_tidytimecourse.csv") %>% 
  filter(time == 5, sample_condition == "wt_polya", is_maternal) %>% 
  rename(stability = log2FC) %>% 
  select(Gene_ID, stability) %>% 
  rename(gene_id = Gene_ID)


stability <- read_csv("../../data/19-04-29-AllDecayProfiles/decay_profiles_and_seqdata.csv") %>% 
  filter(datatype == "slam-seq")


stability <- stability %>% 
  select(gene_id, decay_rate) %>% 
  rename(stability = decay_rate) %>% 
  bind_rows(fishstability)

fdata <- inner_join(stability, pathways)

# scale the mRNA stability ------------------------------------------------

fdata <- fdata %>% 
  group_by(specie) %>% 
  mutate(stability = as.numeric(scale(stability))) %>% 
  ungroup()

write_csv(fdata, "results_data/pathways_mRNAstability.csv")
