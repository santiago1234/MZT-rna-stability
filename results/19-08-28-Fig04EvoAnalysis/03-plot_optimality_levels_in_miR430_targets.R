# plot the optimility levels as measured by PLS1 of miR-430 targets
# and not targets in fish and xenopus
library(tidyverse)
library(ggthemes)
library(ggforce)

theme_set(theme_tufte(base_family = "Helvetica"))

opt <- read_csv("results_data/optimality_information.csv")
foldchange <- read_csv("results_data/fold_change_late_vs_early.csv")

datum <- inner_join(opt, foldchange) %>% 
  mutate(utrlen = str_length(`3utr`)) %>% 
  filter(!is.na(`3utr`), !is.na(coding), !is.infinite(log2FC))

mir_seed_type <- function(utr_seq) {
  if (str_detect(utr_seq, pattern = "AGCACTTA")) {
    return("8-mer")
  }
  
  if (str_detect(utr_seq, pattern = "AGCACTT") | str_detect(utr_seq, pattern = "GCACTTA")) {
    return("7-mer")
  }
  
  if (str_detect(utr_seq, pattern = "GCACTT") | str_detect(utr_seq, pattern = "AGCACT")) {
    return("6-mer")
  }
  
  "no seed"
  
}

# tidy data for plotting
datum_plot <- 
  datum %>% 
  group_by(specie) %>% 
  mutate(
    `miR-430 seed strength` = map_chr(`3utr`, mir_seed_type),
    `number of miR-430 seeds` = map_dbl(miR430, ~if_else(. > 3, 3, .)),
    stability = ntile(x = log2FC, n = 4),
  ) %>% 
  filter(stability == 1) %>% 
  select(gene_id, specie, PLS1,`miR-430 seed strength`, `number of miR-430 seeds`) %>% 
  gather(key = reg, value = typer, -gene_id, -specie, -PLS1) %>% 
  mutate(
    reg = factor(reg, levels = c("number of miR-430 seeds", "miR-430 seed strength")),
    typer = factor(typer, levels = c("0", "1", "2", "3", "no seed", "6-mer", "7-mer", "8-mer"))
  )


dp_median <- 
  datum_plot %>% 
  group_by(reg, typer, specie) %>% 
  summarise(median_pls = median(PLS1)) %>% 
  filter(!is.na(median_pls))

plt_colors <- c("grey", rep("grey20", 2), "grey", rep("grey20", 3))

datum_plot %>% 
  ggplot(aes(x=typer, y=PLS1, color=typer)) +
  geom_sina(shape=16, size=1/3, alpha=.99) +
  geom_rangeframe(size=1/5, color="black") +
  geom_errorbar(
    data = dp_median,
    aes(ymin = median_pls, ymax = median_pls, y = median_pls),
        color = "black",
    size = 1/5
  ) +
  facet_grid(specie~reg, scales = "free_x", space = "free_x") +
  scale_color_manual(values = c("#bf812d", "#35978f", "#01665e","#003c30", "#bf812d", "#35978f", "#01665e","#003c30")) +
  
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 9, angle = 30, hjust = 1)
  ) +
  labs(
    y = "optimality level: PLS1",
    x = NULL
  )

ggsave("figures/03_optimality_level_miR430targets.pdf", height = 3, width = 3)


# statistics --------------------------------------------------------------

datum %>% 
  group_by(specie) %>% 
  mutate(
    stability = ntile(x = log2FC, n = 4),
    miR430 = map_dbl(miR430, ~if_else(. > 1, 2, .))
  ) %>% 
  filter(stability == 1) %>% 
  nest() %>% 
  mutate(
    fit = map(data, ~lm(PLS1 ~ miR430, data = .)),
    tfit = map(fit, broom::tidy)
  ) %>% 
  unnest(tfit) %>% 
  select(-data, -fit) %>% 
  knitr::kable()
