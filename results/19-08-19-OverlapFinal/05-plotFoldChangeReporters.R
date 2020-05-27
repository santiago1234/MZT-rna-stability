library(tidyverse)
library(broom)
library(ggthemes)
library(scales)
library(ggforce)
library(optimalcodonR)
library(cowplot)

theme_set(tidybayes::theme_tidybayes())


# get the predicted stability ---------------------------------------------
optimal_Seq <- 'GACATCTTTGGCTTTGAGAACTTTGAGGTCAACCGCTTTGAGCAGTTCAACATTAACTATGCAAACGAGAAGCTTCAGGAGTATTTCAACAAGCACATTTTCTCACTGGAGCAGCTTGAGTTCAGGAAGGTGCAGCATGAGCTGGAGGAGGCTCAGGAGAGAGCTGACATCGCCGAGTCCCAGGTCAACAAGCTCAGAGCTAAAAGCCGTGAATTTGGAAAGGGTAAAGAGGCTGAGGAGGCTGACTCCTTCGACTATAAGAGCTTCTTCGCCAAGGTTGGGCTGTCCGCCAAGACTCCTGATGACATCAAGAAGGCTTTTGCTGTCATTGACCAGGACAAGAGCGGCTTCATTGAGGAGGATGTGGAGGACTCCCTCTGTGAGGCCAAAGAGCTGTTCATCAAGACAGTCAAGCACTTCGGTGAGGACGCTGATAAGATGCAGCCTGATGAGTTCTTTGGGATTTTCGACCAGTTCTTGCGTATCCCCAAGGAGCAGGGCTTCCTGTCGTTCTGGAGAGGAAACTTGGCCAACGTCATCAGATACTTCCCCACACAGGCCCTCAACTTTGCTTTCAAGGACAAGTACAAGAAGGTCTTCGACATCACAGACAAGCTGGAGAACGAGCTGGCCAATAAGGAGGCTTTCCTCAGACAGATGGAGGAGAAGAACAGGCAGTTGCAGGAGCGGCTTGAGTTGGCAGAGCAGAAGCTCCAGCAG'
nonOptimal_Seq <- 'ACATCTTTGGCTTTGAGAACTTTGAGGTCAACCGCTTTGAGCAGTTCAACATTAACTATGCAAACGAGAAGCTTCAGGAGTATTTCAACAAGCACATTTTCTCACTGGAGCAGCTTGAGTTCAGGAAGGTGCAGCATGAGCTGGAGGAGGCTCAGGAGAGAGCTGACATCGCCGAGTCCCAGGTCAACAAGCTCAGAGCTAAAAGCCGTGAATTTGGAAAGGGTAAAGAGGCTGAGGAGGCTGACTCCTTCGACTATAAGAGCTTCTTCGCCAAGGTTGGGCTGTCCGCCAAGACTCCTGATGACATCAAGAAGGCTTTTGCTGTCATTGACCAGGACAAGAGCGGCTTCATTGAGGAGGATGTGGAGGACTCCCTCTGTGAGGCCAAAGAGCTGTTCATCAAGACAGTCAAGCACTTCGGTGAGGACGCTGATAAGATGCAGCCTGATGAGTTCTTTGGGATTTTCGACCAGTTCTTGCGTATCCCCAAGGAGCAGGGCTTCCTGTCGTTCTGGAGAGGAAACTTGGCCAACGTCATCAGATACTTCCCCACACAGGCCCTCAACTTTGCTTTCAAGGACAAGTACAAGAAGGTCTTCGACATCACAGACAAGCTGGAGAACGAGCTGGCCAATAAGGAGGCTTTCCTCAGACAGATGGAGGAGAAGAACAGGCAGTTGCAGGAGCGGCTTGAGTTGGCAGAGCAGAAGCTCCAGCAGG'
predictor <- optimalcodonR::predict_stability('human')

datos <- tibble(
  Optimality = c('Optimal', 'Non-optimal'),
  seq = c(optimal_Seq, nonOptimal_Seq)
) %>% 
  mutate(
    predicted_stability = map_dbl(seq, predictor) %>% unscale_decay_to_mouse()
  )


testing %>% 
#  filter(specie == 'human') %>% 
  mutate(stability = unscale_decay_to_mouse(decay_rate)) %>% 
  select(stability) %>% 
  filter(between(stability, 1.7, 20)) %>% 
  ggplot(aes(x = stability)) +
  geom_density(fill = 'grey', color = NA) +
  geom_vline(xintercept = datos$predicted_stability[1], color = 'red') +
  geom_vline(xintercept = datos$predicted_stability[2], color = 'blue')
  
ggsave('figures/fold-change-hist.pdf', height = 1.5, width = 4)
datum <- readxl::read_excel("results_data/Human_final_table.xlsx")

datum <- 
  datum %>% 
  mutate(
    Optimality = map_chr(Optimality, ~if_else(. == 1, "Optimal", "Non-optimal")),
    Optimality = factor(Optimality, levels = c("Optimal", "Non-optimal")),
    miR = map_chr(miR, ~if_else(. == 1, "Seed", "no Seed")),
    UTR = str_detect(ID, "205") %>% map_chr(~if_else(., "weak seed", "strong seed")),
    Replicate = map_chr(Replicate, ~if_else(. == 1, "replicate A", "replicate B"))
  )


resultados <- 
  datum %>% 
  group_by(UTR, Optimality) %>% 
  nest() %>% 
  mutate(
    fit = map(data, ~lm(log2(value) ~ Replicate + miR, data = .)),
    tidy_fit = map(fit, tidy, conf.int = T)
  ) %>% 
  unnest(tidy_fit)


resultados %>% 
  inner_join(datos) %>% 
  filter(term == 'miRSeed') %>% 
  ggplot(aes(x = predicted_stability, y = estimate, color = Optimality)) +
  geom_point(shape = 16, alpha = .9, size = 1) +
  scale_color_manual(values = c('blue', 'red')) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = .1) +
  geom_hline(yintercept = 0, linetype = 3) +
  scale_x_continuous(breaks = round(datos$predicted_stability, 2)) +
  labs(
    y = 'log2 fold change'
  ) +
  facet_wrap(. ~ UTR) +
  theme(legend.position = 'none', axis.text.x = element_text(angle = 20))

ggsave('figures/fold-change-optimalVSnon-optimal_majo.pdf', height = 1.5, width = 2)
