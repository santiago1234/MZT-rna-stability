library(tidyverse)
library(ggthemes)
theme_set(theme_tufte(base_family = 'Helvetica'))


eGFPseq <- "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAG"
egfp_path <- read_csv("results_data/eGFP_optimization_path.csv") %>% 
  select(-X1)


egfp_path <- 
  egfp_path %>% 
  select(-fitness) %>% 
  gather(key = specie, value = stability, -iteration, -sequences, -optimization)

eGFPseq <- egfp_path %>% 
  filter(sequences == eGFPseq, optimization == "maximization") %>% 
  distinct()

egfp_path %>%
  ggplot(aes(x=iteration, y = stability)) +
  geom_line(color="grey") +
  geom_point(
    data = eGFPseq,
    color = "red"
  ) +
  ggrepel::geom_text_repel(
    data = eGFPseq,
    color = "red",
        aes(label="eGFP")
  ) +
  facet_grid(~specie) +
  geom_rangeframe()

ggsave("figures/gfp_optimization.pdf", height = 3, width = 6)
# print the best GFP ------------------------------------------------------

cat(">eGFP\n")
cat(eGFPseq$sequences[1])

cat("\n>supremaGFP\n")
egfp_path %>% 
  filter(specie == "fish | embryo mzt | aamanitin polya") %>% 
  arrange(-stability) %>% 
  slice(1:1) %>% 
  pull(sequences) %>% 
  cat()

cat("\n>infimaGFP\n")
egfp_path %>% 
  filter(specie == "fish | embryo mzt | aamanitin polya") %>% 
  arrange(stability) %>% 
  slice(1:1) %>% 
  pull(sequences) %>% 
  cat()
cat("\n")


# plot predictions for silent reporters -----------------------------------

reporters <- read_csv("results_data/predict_silent_reporters.csv")

reporters %>% 
  ggplot(aes(x=predicted_stability, y=protein, color=p_optimal)) +
  geom_point() +
  geom_rangeframe(color="black") +
  scale_color_gradient(low="blue", high = "red") +
  ggrepel::geom_text_repel(aes(label=p_optimal)) +
  labs(
    y = "Fluorescent Intensity",
    x = "predicted mRNA stability",
    title = "Synonimous Reporters"
  ) +
  theme(legend.position = "none")
ggsave("figures/synonimous_reporters.pdf", height = 2, width = 2.5)

