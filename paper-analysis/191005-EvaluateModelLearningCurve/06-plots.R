library(tidyverse)
library(ggthemes)
library(ggpointdensity)

theme_set(theme_tufte(base_family = "Helvetica"))


# predictions test data ---------------------------------------------------

testdata <- read_csv("results-data/test_data_predictions.csv")

testdata %>% 
  ggplot(aes(x=predicted, y=observed)) +
  geom_pointdensity(size=1/3, shape=16, alpha=.99) +
  scale_color_viridis_c(option = "A") +
  geom_rangeframe() +
  ggpubr::stat_cor() +
  geom_abline(linetype=2, size=1/5) +
  labs(
    title = paste0("predictions in test set genes\nn = ", nrow(testdata)),
    x = "predicted stability",
    y = "observed stability"
  )
ggsave("figures/predictions_test_data.pdf", width = 5, height = 2.5)


# break the data by specie ------------------------------------------------


orders <- c(
  "human hela",
  "human RPE",
  "human k562",
  "human 293t",
  "mouse mES cells",
  "fish embryo mzt",
  "xenopus embryo mzt"
)

testdata %>% 
  mutate(
    id = paste0(specie, " ", cell_type),
    id = factor(id, levels = orders)
  ) %>% 
  ggplot(aes(x=predicted, y=observed)) +
  geom_pointdensity(size=1/3, shape=16, alpha=.99) +
  scale_color_viridis_c(option = "A") +
  geom_rangeframe(size=1/3) +
  ggpubr::stat_cor(size=2) +
  labs(
    x = "predicted stability",
    y = "observed stability"
  ) +
  facet_grid(~id)
ggsave("figures/predictions_test_data_by_specie.pdf", width = 14, height = 2)


# r2 boostrap intervals by profile ----------------------------------------

r2_profiles <- read_csv("results-data/boostrap_r2_byProfile.csv")

r2_profiles %>%
  mutate(
    id = paste(cell_type, datatype),
    specie = factor(specie, levels = c("human", "mouse", "fish", "xenopus"))
  ) %>% 
  ggplot(aes(x=r2_estimate, y=reorder(id, r2_estimate))) +
  geom_point(shape=1, alpha=.99) +
  geom_errorbarh(aes(xmin=r2_ci_low, xmax=r2_ci_high), height=0) +
  geom_rangeframe(sides = "l", size=1/3) +
  coord_cartesian(xlim = c(0, .35)) +
  facet_grid(specie~., scales = "free_y", space = "free_y") +
  labs(
    x = "R2",
    y = "mRNA stability profile"
  )
