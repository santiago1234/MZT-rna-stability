library(tidyverse)

# do mir430 interacts with optimality? ------------------------------------

optimality <- read_csv("results_data/optimality_information.csv") %>% 
  mutate(utr_len = str_length(`3utr`))
fc <- read_csv("results_data/fold_change_late_vs_early.csv") %>% 
  filter(!is.na(log2FC), !is.infinite(log2FC))

datum <- inner_join(optimality, fc) %>% 
  filter(specie == "fish")

fit <- lm(log2FC ~ PLS2*miR430 + PLS1*miR430 + utr_len, data = datum)


datum2 <- tibble(
  PLS2 = seq(-5, 5, length.out = 50),
  PLS1 = 0,
  utr_len = median(datum$utr_len)
) %>% 
  crossing(tibble(miR430 = 0:4))

datum2$prediction <- predict(fit, newdata = datum2)

datum2 %>% 
  ggplot(aes(x=PLS2, y = prediction, color=miR430)) +
  geom_point()

fit %>% 
  summary()

# the data does not support an interacting effect, but more an additive effect
# p.value = 0.0439