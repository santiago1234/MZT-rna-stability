library(tidyverse)
library(ggthemes)

theme_set(theme_tufte(base_family = "Helvetica"))

optimality <- read_csv("../../results/19-08-19-OverlapFinal/results_data/pls_species.csv")
fc_timepoints <- read_csv("results-data/predicted-log2fold-timepoints.csv")

d <- inner_join(fc_timepoints, optimality) %>% 
  select(gene_id, time, PLS1, PLS2, log2FC, specie)


# which genes do contain miR-430 sites ------------------------------------

xen <- read_csv("../../data/19-06-12-XenopusData/regulatory_pathways_xenopus.csv") %>% 
  select(Gene_ID, MiR430) %>% 
  rename(gene_id = Gene_ID)

fish <- read_csv("../../results/19-02-24-OverlapPathwaysFig3/results_data/regulatory_pathways_matrix.csv") %>% 
  select(Gene_ID, MiR430) %>% 
  rename(gene_id = Gene_ID)


d <- 
  bind_rows(xen, fish) %>% 
  mutate(miR430 = MiR430 > 0) %>% 
  inner_join(d)


d_xen <- d %>% 
  filter(
    time %in% 8:12, specie == "xenopus",
    !is.na(miR430)
    ) 

mean_log2fc <- d_xen %>% 
  group_by(time, miR430) %>% 
  summarise(mlf = mean(log2FC))


# make time a factor for ploting ------------------------------------------


d_xen %>% 
  ggplot(aes(x = PLS1, y = log2FC, color = miR430)) +
  geom_point(shape = 16, size = .5, alpha = .9) +
  geom_smooth(method = "lm", se = F, color = "blue",alpha = .5, size = .5) +
  geom_hline(data = mean_log2fc, aes(yintercept = mlf), linetype = 2, size = .4) +
  scale_color_manual(values = c("grey", '#04591A')) +
  facet_grid(time~miR430) +
  coord_cartesian(ylim = c(-7, 6)) +
  theme(axis.line = element_line(size = .2), legend.position = "none") +
  labs(
    x = "Codon optimality level",
    subtitle = "Xenopus maternal genes"
  )

ggsave("figures/xen-mir-optimality-overlap-by-timepoint.pdf", height = 5, width = 2.4)



# do the same plot in fish data -------------------------------------------

d_fish <- d %>% 
  filter(
    time %in% 5:10, specie == "fish",
    !is.na(miR430)
  ) 

mean_log2fc <- d_fish %>% 
  group_by(time, miR430) %>% 
  summarise(mlf = mean(log2FC))


# make time a factor for ploting ------------------------------------------


d_fish %>% 
  ggplot(aes(x = PLS1, y = log2FC, color = miR430)) +
  geom_point(shape = 16, size = .5, alpha = .9) +
  geom_smooth(method = "lm", se = F, alpha = .5, color = "blue", size = .5) +
  geom_hline(data = mean_log2fc, aes(yintercept = mlf), linetype = 2, size = .4) +
  ggpubr::stat_cor(size = 1.5, color = "black") +
  coord_cartesian(ylim = c(-7, 6)) +
  scale_color_manual(values = c("grey", '#04591A')) +
  facet_grid(time~miR430) +
  theme(axis.line = element_line(size = .2), legend.position = "none") +
  labs(
    x = "Codon optimality level",
    subtitle = "Zebrafish maternal genes"
  )

ggsave("figures/fish-mir-optimality-overlap-by-timepoint.pdf", height = 4, width = 2.4)

