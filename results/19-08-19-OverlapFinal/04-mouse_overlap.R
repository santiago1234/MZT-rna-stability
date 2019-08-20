library(tidyverse)
library(ggthemes)
library(scales)

theme_set(theme_tufte(base_family = "Helvetica"))

# mouse analysis ----------------------------------------------------------

data <- bind_rows(
  read_csv("../19-07-18-PredictiveModelWithM6AandMicroRNAs/results_data/training_set.csv"),
  read_csv("../19-07-18-PredictiveModelWithM6AandMicroRNAs/results_data/validation_set.csv")
) %>% 
  filter(specie == "mouse")

mouse_optimality <- read_csv("results_data/pls_species.csv") %>% 
  filter(specie == "mouse")

data <- data %>% 
  select(gene_id, decay_rate, m6A, microRNAsites) %>% 
  gather(key = pathway, value = sites, -gene_id, -decay_rate) %>% 
  mutate(target = sites > 0) %>% 
  inner_join(mouse_optimality)
  

data <- data %>% 
  group_by(pathway, target) %>% 
  mutate(optimality = ntile(-PLS1, 4)) %>% 
  ungroup()

data %>% 
  filter(optimality %in% c(1, 4)) %>% 
  ggplot(aes(x=target, y=decay_rate, color=as.character(optimality))) +
  geom_tufteboxplot() +
  geom_rangeframe(sides = "l", color="black", alpha=2/3) +
  facet_grid(~pathway) +
  scale_x_discrete(labels = c("no targets", "targets")) +
  scale_color_manual(values = c("red", "blue")) +
  theme(
    legend.position = "none"
  ) +
  labs(
    y = "mRNA stability\ndecay rate scaled",
    title = "combinatorial code in mESCs",
    x = NULL
  )

ggsave("figures/mouse.pdf", height = 2, width = 4)

data %>% 
  group_by(pathway, target) %>% 
  nest() %>% 
  mutate(
    fit = map(data, ~lm(decay_rate ~ PLS1, data = .)),
    tf = map(fit, broom::tidy)
  ) %>% 
  unnest(tf) %>% 
  filter(term == "PLS1")

