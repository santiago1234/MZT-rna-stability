library(tidyverse)
library(ggthemes)

theme_set(theme_tufte())


# data for overlap m6a ----------------------------------------------------

m6Agenes <- read_csv("../19-05-13-OptimalityOverlapMouseHuman/results_data/mmc5.csv") %>% 
  filter(Region == "3' UTR") %>% 
  pull(Gene) %>% 
  unique()

mappings <- read_csv("../19-05-13-OptimalityOverlapMouseHuman/results_data/mappings_ensemble_id_to_gene_names.csv")


m6a_targets <- mappings %>% 
  mutate(m6a = Name %in% m6Agenes) %>% 
  filter(m6a) %>% 
  pull(gene_id)

humandata <- read_csv("../19-05-13-OptimalityOverlapMouseHuman/results_data/pls_components_human.csv") %>% 
  mutate(m6A = gene_id %in% m6a_targets)

# fish data m6a -----------------------------------------------------------

programs <- read_csv("../19-02-24-OverlapPathwaysFig3/results_data/regulatory_pathways_matrix.csv") %>% 
  filter(is_maternal == TRUE) %>% 
  select(Gene_ID, m6A, MiR430, PLS1, PLS2)

# load log2 fc  -----------------------------------------------------------

log2fc <- read_csv("../../data/19-02-05-FoldChangeData/data/log2FC_earlyVSlate_tidytimecourse.csv") %>% 
  filter(is_maternal, sample_condition == "wt_ribo") %>% 
  filter(time == 6)

fish_data <- inner_join(programs, log2fc)

# plot m6a ----------------------------------------------------------------
# let's plot m6a for human

hum_m6a <- 
  humandata %>% 
  filter(datatype == "slam-seq") %>% 
  rename(stability = decay_rate, Gene_ID = gene_id) %>% 
  mutate(specie = "human") %>% 
  filter(stability > -.4)

fish_m6a <- fish_data  %>% 
  mutate(specie = "fish") %>% 
  rename(
    stability = log2FC
  ) %>% 
  filter(stability < 3, stability > -6) # outliers

data_m6a <- bind_rows(hum_m6a, fish_m6a) %>% 
  group_by(specie) %>% 
  mutate(optimality = ntile(-PLS1, 5)) %>% 
  filter(optimality %in% c(1, 5)) %>% 
  mutate(optimality = as.character(optimality)) %>% 
  mutate(
    m6A = if_else(m6A, "m6A targets", "no targets"),
    m6A = factor(m6A, levels = c("no targets", "m6A targets"))
  ) %>% 
  ungroup() %>% 
  group_by(specie) %>% 
  mutate(
    stability = as.numeric(scale(stability))
  )

data_m6a %>% 
  ggplot(aes(m6A, stability, fill=optimality)) +
  geom_boxplot(alpha=2/3, outlier.shape = NA, size=1/4) +
  facet_wrap(specie~.) +
  scale_fill_manual(values = c("red", "blue"),labels = c("optimal genes", "non-optimal genes")) +
  labs(
    y = "mRNA stability\n(scaled)",
    title = "m6A overalp with optimality",
    x = NULL
  )
ggsave("overlap_m6a_human_and_fish.pdf", width = 5, height = 2.5)

# miR-430 -----------------------------------------------------------------

xenopus_data <- read_csv("results_data/xenopus_verlap_data.csv") %>% 
  select(-jgi_id, -external_gene_name, -sample_condition, -in_top_4k_expressed_at_0hrs) %>% 
  mutate(specie = "xenopus")


# human mir2 --------------------------------------------------------------

quantifications <- read_csv("../../data/19-05-08-MammalianMircroRNAs/rna-seq-quants-mammlian-micro-rnas.csv")

micro_rnas <- read_csv("../../data/19-05-08-MammalianMircroRNAs/mircro_rnas_tidy.csv")

pls <- read_csv("../19-05-13-OptimalityOverlapMouseHuman/results_data/pls_components_human.csv") %>% 
  filter(cell_type == "hela")

mir1 <- quantifications %>% 
  filter(
    sample_name %in% c("mrna_mir1_32hr")
  ) %>% 
  left_join(filter(micro_rnas, `Representative miRNA` != "hsa-miR-155-5p")) %>% 
  replace_na(list(`Conserved sites total` = 0)) %>% 
  inner_join(pls)


human_mir1 <- mir1 %>% 
  select(gene_id, expression_level, PLS2, `Conserved sites total`) %>% 
  mutate(
    optimality = ntile(-PLS2, 5),
    mirSites = `Conserved sites total`,
    stability = log(expression_level + 1),
    micro = "miR-1 human"
  )

fish_mir <- fish_data %>% 
  mutate(
    stability = log2FC,
    mirSites = MiR430,
    micro = "miR-430 fish",
    optimality = ntile(-PLS1, 5)
  )

xen_mir <- xenopus_data %>% 
  mutate(
    stability = log2FC,
    optimality = ntile(-PLS1, 5),
    mirSites = MiR430,
    micro = "miR-427 xenopus"
  )

# make plot ---------------------------------------------------------------

data_mirs <- bind_rows(xen_mir, human_mir1, fish_mir) %>% 
  mutate(
    mirSites = map_chr(mirSites, ~if_else(. >= 1, ">=1", as.character(.))),
    mirSites = factor(mirSites, levels = c("0", ">=1"))
  ) %>% 
  group_by(micro) %>% 
  mutate(stability = as.numeric(scale(stability))) %>% 
  ungroup() %>% 
  filter(optimality %in% c(1, 5)) %>% 
  mutate(
    optimality = as.character(optimality),
    micro = factor(micro, levels = c("miR-430 fish", "miR-427 xenopus", "miR-1 human"))
  ) 

data_mirs %>% 
  ggplot(aes(x=mirSites, y = stability, fill=optimality)) +
  geom_boxplot(alpha=2/3, outlier.shape = NA, size=1/4) +
  scale_fill_manual(values = c("red", "blue"),labels = c("optimal genes", "non-optimal genes")) +
  facet_wrap(micro~., nrow = 1) +
  labs(
    x = "micro-RNA target sites",
    title = "Overlap optimality with microRNAs",
    y = "mRNA\nstability/level (scaled)"
  ) +
  coord_cartesian(ylim = c(-3.5, 3))
ggsave("overlap_micro_human_and_fish.pdf", width = 6, height = 2.5)
