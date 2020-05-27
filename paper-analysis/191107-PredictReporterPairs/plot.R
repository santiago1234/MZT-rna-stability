library(tidyverse)
library(ggthemes)

theme_set(theme_tufte(base_family = "Helvetica"))

resutls <- read_csv("reporters_with_predicted_stability_all_sepcies.csv") %>% 
  separate(gene_id, into = c("reporter", "optimality"), sep = "\\|")

# filter the data for the plot --------------------------------------------

fish <- resutls %>% 
  filter(specie == "fish", datatype == "aamanitin polya")

human <- resutls %>% 
  filter(specie == "human", cell_type == "RPE")

datos <- bind_rows(fish, human)


# plot --------------------------------------------------------------------

datos %>% 
  ggplot(aes(x=predicted_stability, y=reorder(reporter, predicted_stability, median), color = optimality)) +
  geom_rangeframe(color="black", size=1/5) +
  geom_line(color="black", linetype=3, size=1/5) +
  geom_point(shape=16, alpha=.9) +
  scale_color_manual(values = c("blue", "red")) +
  facet_grid(~specie)

ggsave("predictions_reporters.pdf", height = 2, width = 6)


# slope plot --------------------------------------------------------------

tags <- datos %>% 
  filter(reporter != "embo2016-B", optimality == "non-optimal")

datos %>% 
  filter(reporter != "embo2016-B") %>% 
#  select(reporter, optimality, specie, predicted_stability) %>% 
  ggplot(aes(x=optimality, y=predicted_stability, group=reporter)) +
  geom_point(shape=16, alpha=.9) +
  ggrepel::geom_text_repel(aes(label=round(predicted_stability, 2)), size=4, color="grey") +
  geom_text(data=tags, aes(label=reporter)) +
  geom_line() +
  facet_grid(~specie) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(y=NULL, x=NULL)

ggsave("predictions_reporters2.pdf", height = 3, width = 4)
