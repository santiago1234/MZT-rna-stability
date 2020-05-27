library(tidyverse)
library(ggthemes)

# GOAL:
# Make a plot where you only show the MA genes

theme_set(theme_tufte(base_family = "Helvetica"))

fish_fc <- read_csv("../../data/19-02-05-FoldChangeData/data/log2FC_earlyVSlate_tidytimecourse.csv") %>% 
  filter(
    sample_condition == "wt_ribo",
    is_maternal,
    time == 6
  ) %>% 
  select(-is_maternal, -sample_condition, -time)


optimality <- read_csv("results_data/pls_species.csv") %>% 
  rename(Gene_ID = gene_id)

data <- inner_join(fish_fc, optimality)


# load the definition of MA etc. ------------------------------------------

mir430codes <- read_csv("../../data/19-01-17-EMBO2016DATA/datasets/Groups/mir430_genes_targets.csv") %>% 
  rename(Gene_ID = gene_id)

only_ma <- inner_join(data, mir430codes) %>% 
  filter(MicroArray != "No")
  
only_ma %>% 
  ggplot(aes(x=PLS1, y=log2FC)) +
  geom_point(shape=16, alpha=.98, size=1/2) +
  ggpubr::stat_cor() +
  geom_smooth(method = MASS::rlm) +
  geom_rangeframe(size=1/4) +
  labs(
    title = "miR-430 targets (Microarray detected)",
    x = "PLS1\nCodon Optimality level",
    y = "log2 Fold Change"
  ) +
  theme(axis.ticks = element_line(size=1/5))

ggsave("figures/mzt_mir430_overlap-MA.pdf", height = 2, width = 3)
