library(tidyverse)
library(ggthemes)
library(scales)

optimization <- 
  list.files("results-data/", pattern = "ENS.*-optimization.csv", full.names = T) %>% 
  map_df(read_csv) %>% 
  select(gene_id, sequences, iteration, fitness, optimization)



optimization %>% 
  filter(iteration == 0, optimization == "minimization") %>% 
  group_by(gene_id) %>% 
  slice(1:1) %>% 
  ungroup() %>% 
  mutate(optimality = ntile(fitness, 5))  %>% 
  select(gene_id, optimality) %>% 
  ungroup() -> level_optimal

optimization <- 
  optimization %>% 
  inner_join(level_optimal)

optimization %>% 
  mutate(lenseq = str_length(sequences), id=paste(gene_id, optimization)) %>% 
  filter(optimality %in% c(5)) %>% 
  group_by(id, optimization, iteration) %>% 
  slice(1:1) %>% 
  ggplot(aes(x=iteration+1, y=fitness, color=lenseq, group=id)) +
  geom_line(size=1/2, position="jitter") +
  geom_rangeframe(color="black", size=1/5) +
  theme_tufte() +
  scale_x_continuous(breaks = scales::breaks_width(width = 25)) +
  scale_y_continuous(breaks = c(-2.5, 0, 2.5)) +
  scale_color_viridis_c(
    option = "D",
    trans="log10",
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  coord_cartesian(ylim = c(-3, 2.5)) +
  labs(
    y = "predicted stability",
    x = "Iteration"
  )

ggsave("optimization-genes.pdf", height = 2, width = 3)


# get the plot for the distribution of endegenous genes -------------------

humangenes_endoStability <- codonr::train_set %>% 
  filter(specie == "human") %>% 
  group_by(gene_id) %>% 
  slice(1:1)

humangenes_endoStability %>%
  ggplot(aes(x=-decay_rate)) +
  geom_histogram(bins=50, fill="grey", color="white", size=1/3) +
  theme_tufte() +
  geom_hline(yintercept = c(200, 400, 600), color="white", size=1/5) +
  scale_x_continuous(expand = c(0,0), breaks = c(-2.5, 0, 2.5)) +
  coord_cartesian(xlim = c(-3, 2.5))
ggsave("endo-genes-distribution.pdf", height = 1, width = 2)
