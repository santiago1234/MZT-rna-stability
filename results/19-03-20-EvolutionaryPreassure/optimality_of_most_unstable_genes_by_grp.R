library(tidyverse)
library(ggthemes)

theme_set(theme_tufte(base_family = "Helvetica"))
# get the optimality for fish and xenopus\ --------------------------------
top_most_unstable <- 1000


fish_data <- read_csv("../19-02-24-OverlapPathwaysFig3/results_data/regulatory_pathways_matrix.csv")

## get log fold change for fish (6 hrs ribo)

fc_fish <- read_csv("../../data/19-02-05-FoldChangeData/data/log2FC_earlyVSlate_tidytimecourse.csv") %>% 
  filter(time == 6, sample_condition == "wt_ribo", is_maternal) %>% 
  select(-sample_condition, -time, -is_maternal)

# tidy the data for visualization -----------------------------------------
xen_time_course <- read_csv("../../data/19-06-12-XenopusData/time_course_xenopus.csv")

xenlog <- 
  xen_time_course %>% 
  filter(time %in% c(2, 8), sample_condition == "wt_polya") %>% 
  spread(key = time, value = expression_quantification) %>% 
  mutate(log2FC = log(`8`/`2`)) %>% 
  select(ensembl_gene_id, log2FC) %>% 
  rename(Gene_ID = ensembl_gene_id)

optimality_xen <- read_csv("../../data/19-06-12-XenopusData/regulatory_pathways_xenopus.csv")

xendat <- 
  inner_join(xenlog, optimality_xen) %>% 
  select(-in_top_4k_expressed_at_0hrs) %>% 
  filter(!is.na(MiR430)) %>% 
  select(Gene_ID, log2FC, PLS1, MiR430) %>% 
  mutate(
    miR430 = MiR430 > 0,
    no_utr_element = MiR430 == 0,
    specie = "xenopus"
  ) %>% 
  filter(!is.infinite(log2FC)) %>% 
  arrange(log2FC) %>% 
  slice(1:top_most_unstable) %>% 
  select(-MiR430) %>% 
  gather(key = grp, value = in_grp, -Gene_ID, -specie, -log2FC, -PLS1) %>% 
  filter(in_grp)

fish_data <- fish_data %>% 
  select(Gene_ID, m6A, MiR430, PLS1) %>% 
  inner_join(fc_fish) %>% 
  mutate(
    specie = "fish",
    miR430 = MiR430 > 0,
    m6A_and_miR430 = miR430 & m6A,
    no_utr_element = !miR430 & !m6A
  ) %>% 
  arrange(log2FC) %>% 
  slice(1:top_most_unstable) %>% 
  select(-MiR430) %>% 
  gather(key = grp, value = in_grp, -Gene_ID, -specie, -log2FC, -PLS1) %>% 
  filter(in_grp)

dset <- bind_rows(xendat, fish_data) %>% 
  mutate(
    grp = factor(grp, levels = c("no_utr_element", "miR430", "m6A", "m6A_and_miR430"))
  )

# make plot ---------------------------------------------------------------

dset %>% 
  ggplot(aes(x=grp, y=PLS1, fill=grp)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~specie, scales="free_y") +
  scale_fill_manual(values = c("grey40", "goldenrod3", "forestgreen", "purple")) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1)) +
  geom_rangeframe() +
  labs(
    x=NULL,
    y = "optimality"
  )

ggsave("figures/optimality_of_most_unstable_genes.pdf", height = 3, width = 4)
# compute p value for xenopus ---------------------------------------------

lm(PLS1 ~ grp, data = xendat) %>% summary()

