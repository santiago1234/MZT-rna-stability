library(tidyverse)
library(ggthemes)
library(ggforce)

theme_set(theme_tufte(base_family = "Helvetica"))

mztdatos <- read_csv("results-data/mzt_predictionsResidualsLog2Fc.csv")


# classify the 3utr according to the kmer type ----------------------------

classify_utr <- function(utr_seq) {
  if (str_detect(utr_seq, "AGCACTTA")) return("8-mer")
  if (str_detect(utr_seq, "(AGCACTT|GCACTTA)")) return("7-mer")
  if (str_detect(utr_seq, "(AGCACT|GCACTT)")) return("6-mer")
  return("no miR seed")
  
}

mztdatos %>% 
  filter(!is.na(`3utr`)) %>% 
  mutate(
    miR430 = map_chr(`3utr`, classify_utr),
    miR430 = factor(miR430, levels = c("no miR seed", "6-mer", "7-mer", "8-mer"))
  ) -> mztdatos



# plot --------------------------------------------------------------------

# get the median to draw a line in the sina plot

residual_median <- 
  mztdatos %>% 
  group_by(specie, miR430) %>% 
  summarise(mediana = median(resid))

mztdatos %>% 
  ggplot(aes(x=miR430, y=resid, color=miR430)) +
  geom_sina(shape=16, size=1/3) +
  geom_errorbar(data=residual_median, aes(y=mediana, ymin=mediana, ymax=mediana), color="grey",  size=1/3) +
  facet_grid(specie~.) +
  coord_cartesian(ylim = c(-4.5, 4.5)) +
  scale_color_viridis_d(option = "C") +
  geom_rangeframe(color="black", size=1/4) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("figures/residualsMir430.pdf", height = 4, width = 2.5)
