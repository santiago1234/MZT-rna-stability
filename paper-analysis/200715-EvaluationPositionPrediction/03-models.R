library(tidyverse)

d <- read_csv("data/counts_by_position.csv")

test <- filter(d, allocation == "testing") %>% 
  select(-allocation)

train <- filter(d, allocation == "training") %>% 
  select(-allocation)


# fit the positional models -----------------------------------------------


d_first <- select(train, decay_rate, contains("_first"))
first_mdl <- lm(decay_rate ~ ., data = d_first)


d_second <- select(train, decay_rate, contains("_second"))
second_mdl <- lm(decay_rate ~ ., data = d_second)

d_third <- select(train, decay_rate, contains("_third"))
third_mdl <- lm(decay_rate ~ ., data = d_third)

d_positional <- select(train, decay_rate, contains("_third"), contains("_first"), contains("_second"))
positional_mdl <- lm(decay_rate ~ ., data = d_positional)

d_all <- select(train, decay_rate, contains("_all"))
all_mdl <- lm(decay_rate ~ ., data = d_all)


# make a tibble with all the models ---------------------------------------

resultados <- 
  tibble(
  model = c("first", "second", "third", "positional", "all"),
  fit = list(first_mdl, second_mdl, third_mdl, positional_mdl, all_mdl)
) %>% 
  mutate(
    tf = map(fit, broom::glance)
  ) %>% 
  unnest(tf)

resultados %>% 
  ggplot(aes(x = reorder(model, adj.r.squared), y = adj.r.squared)) +
  geom_text(aes(label = round(adj.r.squared, 3)), color = "grey30") +
  geom_errorbar(aes(ymin = 0, ymax = adj.r.squared), width = 0, color = "steelblue", alpha =.7) +
  labs(
    x = "model",
    title = "zebrafish data"
  ) +
  theme_light()
ggsave("model-results.pdf", height = 2, width = 4)
