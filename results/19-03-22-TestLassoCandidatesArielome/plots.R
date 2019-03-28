library(tidyverse)


# get data ----------------------------------------------------------------

arielome_results <- list.files("results_data/", pattern = "mdlfit*", full.names = T) %>% 
  map_df(read_csv) %>% 
  filter(parameter == "b_latetimepointTRUE")

lassores <- 
  read_csv("../19-03-18-LassoFindKmers/results_data/candidate_elements_laso_analysis.csv") %>% 
  mutate(
    ef = sign(effect),
    ef = map_chr(ef, ~if_else(. == 1, "stabalizing", "destabalizing"))
  )
  
# this line is to plot the data in the right order
ranking <- 
  lassores %>% 
  group_by(vars) %>% 
  summarise(me = median(effect)) %>% 
  arrange(me) %>% 
  pull(vars)

# append the mir kmer
ranking <- c(ranking, c("AGCACTT", "GCACTTA", "AGCACTTA"))
results <- 
  lassores %>% 
  select(vars, ef) %>% 
  filter(!str_detect(vars, "PLS")) %>% 
  rename(kmer = vars) %>% 
  full_join(arielome_results) %>% 
  replace_na(list(ef = "MiR-430"))


# plots -------------------------------------------------------------------

results %>% 
  mutate(
    ef = factor(ef, levels = c("stabalizing", "destabalizing", "MiR-430")),
    kmer = factor(kmer, levels = ranking)
  ) %>% 
  ggplot(aes(y = kmer, x = Estimate / log(2), color = Estimate / log(2))) +
  geom_point(alpha=.9) +
  geom_point(shape=1, color="black", alpha=.92) +
  geom_vline(xintercept = 0) +
  scale_color_viridis_c() +
  geom_errorbarh(aes(xmin=Q2.5 / log(2), xmax=Q97.5 / log(2)), height = 0) +
  facet_grid(ef~., space = "free_y", scales = "free_y") +
  labs(
    x = "enrichment/depletion (log2 fold change)",
    title = "K-mer effect in reporter library"
  )

ggsave("figures/arielome-effects.pdf", width = 5, height = 7)