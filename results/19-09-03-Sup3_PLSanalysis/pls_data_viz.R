library(tidyverse)
library(ggthemes)
library(ggpointdensity)

theme_set(theme_tufte(base_family = "Helvetica"))

optimality_counts <- function(seq) {
  # returns a tible with the counts for each codon
  # and a percentage
  seq %>%
    codonr::count_codons() %>%
    dplyr::full_join(codonr::optimality_code_embo2016) %>%
    tidyr::replace_na(list(n = 0)) %>%
    dplyr::group_by(optimality) %>%
    dplyr::summarise(n = sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(percent = n / sum(n))
}


# keep only maternal genes for the analysis
maternalGenes <- 
  read_csv("../19-02-24-OverlapPathwaysFig3/results_data/regulatory_pathways_matrix.csv") %>% 
  filter(is_maternal) %>% 
  pull(Gene_ID)


pls <- read_csv("../19-08-19-OverlapFinal/results_data/pls_species.csv") %>% 
  filter(specie == "fish", gene_id %in% maternalGenes)


# load the fish orfs
orf_fish <- read_csv("../../data/19-01-17-Get-ORFS-UTRS-codon-composition/sequence-data/fish_seq_data_cds_3utr.csv") %>% 
  select(-`3utr`) %>% 
  rename(gene_id = ensembl_gene_id)

datum <- pls %>% 
  inner_join(orf_fish)

datum <- 
  datum %>% 
  mutate(counts = map(coding, optimality_counts)) %>% 
  unnest(counts) %>% 
  filter(optimality == "optimal")


# data-viz ----------------------------------------------------------------


datum %>% 
  ggplot(aes(PLS1, PLS2)) +
  geom_point(shape=16, size=1/3, alpha=.9) +
  geom_rangeframe(size=1/5) +
  labs(
    title = "Projection of codon content (64 vars) in two principal components"
  )
ggsave("figures/PLS_decomposition.pdf", height = 2, width = 2)

datum %>% 
  select(PLS1, PLS2, percent) %>% 
  filter(percent < .6) %>% 
  gather(key = component, value = val, -percent) %>% 
  ggplot(aes(x=percent, y=val)) +
  geom_point(shape=16, size=1/3, alpha=.9) +
  geom_rangeframe(size=1/5) +
  ggpubr::stat_cor(size=3) +
  facet_grid(~component) +
  labs(
    x = "% optimal codons (Bazzini et al 2016)",
    y = "component"
  )
ggsave("figures/componentsVsOptPercent.pdf", height = 2, width = 4)


# compare to fc stabilit --------------------------------------------------

fc <- read_csv("../../data/19-02-05-FoldChangeData/data/log2FC_earlyVSlate_tidytimecourse.csv") %>% 
  filter(time == 6, sample_condition == "wt_ribo") %>% 
  rename(gene_id = Gene_ID) %>% 
  select(log2FC, gene_id) %>% 
  inner_join(datum)

fc %>% 
  select(PLS2, percent, log2FC) %>% 
  filter(percent < 0.6) %>% 
  gather(key = per, value = val, -log2FC) %>% 
  ggplot(aes(y = log2FC, x=val)) +
  geom_point() +
  ggpubr::stat_cor() +
  facet_wrap(per ~., scales = "free_x", nrow = 2)


# stabilit profiles -------------------------------------------------------

pls <- read_csv("../19-08-19-OverlapFinal/results_data/pls_species.csv")
testing <- read_csv("../19-04-30-PredictiveModelDecayAllSpecies/19-04-30-EDA/results_data/validation_set.csv")

pls %>% 
  inner_join(testing) %>% 
  filter(!cell_type %in% c("293t", "k562", "RPE"), specie != "xenopus", datatype != "aamanitin ribo") %>% 
  select(PLS1, PLS2, decay_rate, specie) %>% 
  gather(key = component, value = val, -decay_rate, -specie) %>% 
  ggplot(aes(val, decay_rate)) +
  geom_point(shape='.', alpha=.99, size=1) +
  ggpubr::stat_cor(size=2) +
  geom_rangeframe(size=1/5) +
  scale_color_viridis_c(option = "E") +
  facet_grid(.~component) +
  theme(axis.ticks = element_line(size=1/5))

ggsave("figures/componentsVsOutcome.pdf", height = 2, width = 3.5)
