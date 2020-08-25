library(tidyverse)
library(ggthemes)
library(ggforce)

theme_set(theme_tufte(base_family = "Helvetica"))

d <- read_csv("../191010-PredictStabilityInMZT/results-data/mzt_predictionsResidualsLog2Fc.csv")
output_plot <- snakemake@output$plot
output_stats <- snakemake@output$stats

m6a_motif <- "GGACT" # Extended Data Figure 3, 

d <- 
  d %>% 
  select(coding, gene_id, `3utr`,resid, log2FC, specie)


# add m6a motif counts ----------------------------------------------------

d <- 
  d %>% 
  mutate(
    m6a_coding = str_count(coding, m6a_motif),
    m6a_3utr = str_count(`3utr`, m6a_motif),
  ) %>% 
  select(-coding, -`3utr`)


d <- d %>% 
  pivot_longer(c(m6a_coding, m6a_3utr), names_to = "position", values_to = "n")


# compute the p-values ----------------------------------------------------

d %>% 
  group_by(specie, position) %>% 
  nest() %>% 
  mutate(
    fit = map(data, ~lm(resid ~ n, data = .)),
    tfit = map(fit, broom::tidy)
  ) %>% 
  select(-data, -fit) %>% 
  unnest(tfit) %>% 
  filter(term == "n") %>% 
  select(specie, position, p.value) %>% 
  write_csv(output_stats)



d <- d %>% 
  filter(!is.na(n)) %>% 
  mutate(
    n_sites = map_chr(n, ~if_else(. > 1, ">1", as.character(.))),
    n_sites = factor(n_sites, levels = c("0", "1", ">1")),
    position = factor(position, levels = c("m6a_coding", "m6a_3utr"), labels = c("coding", "3' UTR"))
  )

# compute the median to draw line

resid_median <- d  %>% 
  group_by(specie, position, n_sites) %>% 
  summarise(mediana_r = median(resid), n=n())


d %>% 
  ggplot(aes(x = n_sites, y = resid, color = n_sites)) +
  geom_sina(size = .01/2, shape=16, alpha=.9) +
  geom_errorbar(data = resid_median,
                aes(y=mediana_r, x=n_sites, ymin=mediana_r, ymax=mediana_r),
                color="black",
                size=1/5) +
  scale_color_manual(values = c("grey", "#E69F00", "#E69F00")) +
  
  facet_grid(specie~position, scales = "free_y") +
  theme(
    axis.line = element_line(colour = "black", size = .1),
    legend.position = "none"
  )

ggsave(output_plot, height = 3, width = 2.5)
