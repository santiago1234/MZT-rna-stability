library(tidyverse)
library(ggthemes)
library(gridExtra)

theme_set(theme_tufte(base_family = "Helvetica"))

# plot MZT programs fish --------------------------------------------------

pathways <- read_csv("../19-02-24-OverlapPathwaysFig3/results_data/regulatory_pathways_matrix.csv")

pathways <- pathways %>% 
  filter(is_maternal) %>% 
  select(Gene_ID, m6A, MiR430) %>% 
  rename(gene_id = Gene_ID)

optimality <- read_csv("results_data/pls_species.csv")

pathways <- inner_join(pathways, optimality)

# get fold change data ----------------------------------------------------

fc_data <- read_csv("../../data/19-02-05-FoldChangeData/data/log2FC_earlyVSlate_tidytimecourse.csv") %>% 
  filter(
    time == 6,
    is_maternal,
    sample_condition == "wt_ribo"
  ) %>% 
  select(-is_maternal, -sample_condition, -time) %>% 
  rename(gene_id = Gene_ID)

datum <- inner_join(pathways, fc_data)

p1 <- datum %>% 
  ggplot(aes(x=log2FC, color=m6A)) +
  stat_ecdf() +
  scale_color_manual(values = c("grey30", "darkgoldenrod3")) +
  coord_cartesian(xlim = c(-5, 3)) +
  scale_x_continuous(breaks = c(-2, 0, 2)) +
  theme(legend.position  = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  labs(
    y = NULL,
    title = "m6A program",
    x = NULL
  )

p2 <- datum %>% 
  mutate(
    MiR430 = map_chr(MiR430, ~if_else(. > 0, ">1", as.character(.)))
  ) %>% 
  ggplot(aes(x=log2FC, color=MiR430, group=MiR430)) +
  stat_ecdf() +
  scale_color_manual(values = c("forestgreen", "grey30")) +
  scale_x_continuous(breaks = c(-2, 0, 2)) +
  coord_cartesian(xlim = c(-5, 3)) +
  theme(legend.position  = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  labs(
    y = NULL,
    x = NULL,
    title = "miR-430 program"
  )

p3 <- datum %>% 
  mutate(opt = ntile(PLS1, n=5)) %>% 
  ggplot(aes(x=log2FC, color=opt, group=opt)) +
  stat_ecdf() +
  scale_x_continuous(breaks = c(-2, 0, 2)) +
  coord_cartesian(xlim = c(-5, 3)) +
  scale_color_gradient2(low = 'blue', high = "red", mid = "grey30", midpoint = 3) +
  theme(legend.position  = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  labs(
    y = NULL,
    x = NULL,
    title = "codon optimality program"
  )


# mouse programs ----------------------------------------------------------

tr <- read_csv("../19-07-18-PredictiveModelWithM6AandMicroRNAs/results_data/training_set.csv")
te <- read_csv("../19-07-18-PredictiveModelWithM6AandMicroRNAs/results_data/validation_set.csv")
optm <- read_csv("results_data/pls_species.csv")

dtum <- 
  bind_rows(tr, te) %>% 
  filter(specie == "mouse") %>% 
  select(-specie, -cell_type, -datatype, -utrlenlog, -cdslenlog) %>% 
  inner_join(optm)

p1m <- dtum %>% 
  ggplot(aes(x=decay_rate, color=m6A > 0)) +
  stat_ecdf() +
  scale_color_manual(values = c("grey30", "darkgoldenrod3")) +
  coord_cartesian(xlim = c(-3, 2)) +
  scale_x_continuous(breaks = c(-2, 0, 2)) +
  theme(legend.position  = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  labs(
    y = NULL,
    title = "m6A program",
    x = NULL
  )

p2m <- dtum %>% 
  mutate(
    MiR430 = map_chr(microRNAsites, ~if_else(. > 0, ">1", as.character(.)))
  ) %>% 
  ggplot(aes(x=decay_rate, color=MiR430, group=MiR430)) +
  stat_ecdf() +
  scale_color_manual(values = c("forestgreen", "grey30")) +
  scale_x_continuous(breaks = c(-2, 0, 2)) +
  coord_cartesian(xlim = c(-3, 2)) +
  theme(legend.position  = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  labs(
    y = NULL,
    x = NULL,
    title = "miR-291a program"
  )


p3m <- 
  dtum %>% 
  mutate(
    opt = ntile(PLS1, n=5)
  ) %>% 
  ggplot(aes(x=decay_rate, color=opt, group=opt)) +
  stat_ecdf() +
  scale_color_gradient2(low = 'blue', high = "red", mid = "grey30", midpoint = 3) +
  scale_x_continuous(breaks = c(-2, 0, 2)) +
  coord_cartesian(xlim = c(-3, 2)) +
  theme(legend.position  = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  labs(
    y = NULL,
    x = NULL,
    title = "codon optimality program"
  )

pdf("figures//programs_mzt_emsMouse.pdf", width = 6.5, height = 3)
grid.arrange(p1, p2, p3, p1m, p2m, p3m, nrow=2)
dev.off()
