library(tidyverse)

# put together the results ------------------------------------------------

marg_efcts <- list.files("results_data", pattern = "marginal-effec", full.names = T) %>% 
  map_df(read_csv)

## load previous results
lassores <- 
  read_csv("../19-03-18-LassoFindKmers/results_data/candidate_elements_laso_analysis.csv") %>% 
  mutate(
    ef = sign(effect),
    ef = map_chr(ef, ~if_else(. == 1, "stabalizing", "destabalizing")),
    ef = factor(ef, levels = c("stabalizing", "destabalizing"))
  ) 

ranking <- 
  lassores %>% 
  group_by(vars) %>% 
  summarise(me = median(effect)) %>% 
  arrange(me) %>% 
  pull(vars)


marg_efcts <- lassores %>% 
  select(vars, ef) %>% 
  rename(kmer = vars) %>% 
  inner_join(marg_efcts) %>% 
  mutate(kmer = factor(kmer, levels = ranking))
  
# add better names --------------------------------------------------------


meaning_name <- function(x) {
  if (x %in% c("wt1_5h_polyfc", "wt2_6h_polyfc", "wt1_6h_ribofc", "mzdicerDuplex_6h_polyfc")) return("WT")
  if (x %in% c("mzdicer_6h_polyfc", "lna_6h_ribofc")) return("DIZER")
  else ("AA")
}

is_sig <- function(x, y) {
  if (sign(x) == sign(y)) return("*")
  else return("")
}

order_x <- c( "wt1_5h_polyfc", "wt2_6h_polyfc", "wt1_6h_ribofc",
              "mzdicerDuplex_6h_polyfc","mzdicer_6h_polyfc",
              "lna_6h_ribofc", "mcounts")

marg_efcts <- 
  marg_efcts %>% 
  mutate(
    condition = map_chr(rnasample, meaning_name),
    condition = factor(condition, levels = c("WT", "DIZER", "AA")),
    rnasample = factor(rnasample, levels = order_x),
    sig = map2_chr(ci_l, ci_u, is_sig)
  )

marg_efcts %>% 
  ggplot(aes(x=rnasample, y = kmer, fill = estimate)) +
  geom_tile() +
  geom_text(aes(label = sig)) +
  scale_fill_gradient2() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  facet_grid(ef  ~ condition, scales = "free", space = "free") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
    strip.text.x = element_text(size = 7)
  ) +
  labs(
    title = "k-mer effect in log2Fold Change\nlate/early",
    subtitle = "validation data",
    x = NULL,
    y = NULL
  )

ggsave("figures/log2fc_val_data_candidates.pdf", height = 7, width = 3.5)
ggsave("figures/log2fc_val_data_candidates.png", height = 7, width = 3.5)


# save the results --------------------------------------------------------

write_csv(marg_efcts, "results_data/log2fc_validation_test.csv")
