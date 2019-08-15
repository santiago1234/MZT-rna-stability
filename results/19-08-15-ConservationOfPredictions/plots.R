library(tidyverse)
library(ggthemes)

theme_set(theme_tufte(base_family = "Helvetica"))


# plot predictions accross species ----------------------------------------

preds <- read_csv("results_data/orfs_prdictions.csv")

preds %>% 
  ggplot(aes(x=mouse, y=fish)) +
  geom_point(shape=".") +
  ggpubr::stat_cor()


yvar <- preds %>% 
  gather(key = "specie", value = "yvar",-coding)

xvar <- preds %>% 
  gather(key = "specie", value = "xvar",-coding)

inner_join(yvar, xvar, by="coding") %>% 
  ggplot(aes(x=xvar, y=yvar)) +
  geom_point(shape=16, size=1/8, alpha=.8) +
  facet_grid(specie.x~specie.y) +
  ggpubr::stat_cor(size=1.5) +
  labs(
    x = "predicted mRNA stability",
    y = "predicted mRNA stability"
  )

ggsave("preds_across_species.pdf", width = 5, height = 5)
