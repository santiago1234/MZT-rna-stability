library(tidyverse)
library(ggthemes)
library(scales)
library(ggforce)
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
  geom_sina(size=1/5, shape=16, alpha=.99) +
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

# minimalist plot ---------------------------------------------------------

summary_stability <- 
  data %>% 
  mutate(optimality = str_c("q", optimality)) %>%
  group_by(pathway, target, optimality) %>% 
  summarise(
    mediana = median(decay_rate),
    n = n()
  )


data %>% 
  mutate(optimality = str_c("q", optimality)) %>% 
  ggplot(aes(x=optimality, y=decay_rate, color=optimality)) +
  geom_hline(yintercept = 1, size = .2, color = "grey") +
  geom_sina(shape=16, alpha=.99, size=1/4) +
  geom_rangeframe(sides = "l", color="black", size=1/5) +
  geom_errorbar(
    data = summary_stability,
    aes(ymin=mediana, ymax=mediana, y=mediana, x=optimality),
    color="black",
    size=1/7
  ) +
  geom_text(
    data = summary_stability,
    aes(x = optimality, y = 2, label = paste0("n=", n)),
    color = "grey",
    size=1
  ) +
  scale_color_manual(values = c("#ca0020", "#f4a582", "#92c5de", "#0571b0")) +
  facet_grid(~pathway + target) +
  theme(
    legend.position = "none"
  ) +
  labs(
    y = "mRNA stability\ndecay rate scaled",
    title = "combinatorial code in mESCs",
    x = NULL
  )
ggsave("figures/mouse_minimal.pdf", height = 2, width = 5)

