library(tidyverse)
library(ggthemes)
library(ggpointdensity)

theme_set(theme_tufte(base_family = "Helvetica"))


# predictions test data ---------------------------------------------------

testdata <- read_csv("results-data/test_data_predictions.csv")

testdata %>% 
  ggplot(aes(x=predicted, y=observed)) +
  geom_pointdensity(size=1/3, shape=16, alpha=.99) +
  scale_color_viridis_c(option = "E") +
  geom_rangeframe() +
  ggpubr::stat_cor() +
  labs(
    title = paste0("predictions in test set genes\nn = ", nrow(testdata)),
    x = "predicted stability",
    y = "observed stability"
  )
ggsave("figures/predictions_test_data.pdf", width = 4, height = 2.5)




testdata %>% 
  ggplot(aes(x=predicted, y=observed)) +
  geom_point(alpha=.1, shape='.') +
  geom_rangeframe(size=1/3) +
  labs(
    x = "predicted stability",
    y = "observed stability"
  ) +
  facet_wrap(~specie) +
  theme(axis.ticks = element_line(size=1/5)) +
  ggpubr::stat_cor(size=2)
ggsave("figures/predictions_test_data_by_specie.pdf", width = 2, height = 2)
