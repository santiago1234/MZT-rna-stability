library(tidyverse)
library(ggthemes)
library(gridExtra)

theme_set(theme_tufte(base_family = "Helvetica"))


results <- read_csv("results-data/lasso-coefficients.csv") %>% 
  filter(predictor != "(Intercept)") %>% 
  filter(!predictor %in% c("utrlenlog", "len_log10_coding"))



# filter only one specie --------------------------------------------------

# option A
species <- 
  bind_rows(
  filter(results, datatype == "slam-seq"), # human and mouse
  filter(results, specie == "fish", datatype == "aamanitin polya"),
  filter(results, specie == "xenopus"),
) %>% 
  select(-cell_type, -datatype, -data)

# option B

species <- results %>% 
  select(-data) %>% 
  rename(codon = predictor) %>% 
  group_by(specie, codon) %>% 
  summarise(optimality = median(coef))

#  add the genetic code ---------------------------------------------------

stop_codons <- c("TAG", "TAA", "TGA")

genetic_code <- 
  tibble(
    amino = Biostrings::AMINO_ACID_CODE %>% names(),
    amino_3 = Biostrings::AMINO_ACID_CODE
  ) %>% 
  inner_join(tibble(codon = Biostrings::GENETIC_CODE %>% names(), amino = Biostrings::GENETIC_CODE)) %>% 
  select(-amino) %>% 
  bind_rows(tibble(codon = stop_codons, amino_3 = "Stop"))

species <- species %>% 
  inner_join(genetic_code)


# extract the letter position for plotting --------------------------------


order_nucs <- c("T", "C", "A", "G")


species <- 
  species %>% 
  mutate(
    first = str_sub(codon, 1, 1),
    second = str_sub(codon, 2, 2),
    third = str_sub(codon, 3, 3),
    third = factor(third, levels = c("G", "A", "C", "T")),
    second = factor(second, levels = order_nucs),
    first = factor(first, levels = order_nucs),
    # stop codons should be zero optimality
    optimality = map2_dbl(codon, optimality, function(x, y) if_else(condition = x %in% stop_codons, 0.0, y)),
    codon = paste0(codon, " ", amino_3)
)


# draw plot ---------------------------------------------------------------

draw_genetic_code <- function(specie_) {
  species %>% 
    filter(specie == specie_) %>% 
    ggplot(aes(y=third, label=codon, x="A", fill=optimality)) +
    geom_tile() +
    geom_text(nudge_x = -.1, size=2) +
    scale_fill_gradient2(low = "blue3", high = "red2", mid = "white",limits = c(-.06, .06), oob = scales::squish) +
    facet_grid(first~second) +
    theme_linedraw() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 5),
      axis.title.x = element_blank(),
      panel.spacing = unit(0, "lines"),
      axis.ticks = element_blank(),
      strip.text.y = element_text(angle = 0),axis.line = element_line(size=1/10), panel.border = element_rect(size = 1/10)
      
    ) +
    labs(title = specie_)
  
}


gene_plots <- species$specie %>% 
  unique() %>% 
  map(draw_genetic_code)


pdf("genetic-code-optimality.pdf", height = 5, width = 8)
grid.arrange(gene_plots[[1]], gene_plots[[3]], gene_plots[[4]], gene_plots[[2]], top = "Codon optimality code in vertebrates")
dev.off()


# heatmap version by amino acid -------------------------------------------

average_effect_amino_order <- 
  species %>% 
  ungroup() %>% 
  group_by(amino_3) %>% 
  summarise(me = mean(optimality)) %>% 
  arrange(me) %>% 
  pull(amino_3)

species %>% 
  ungroup() %>% 
  mutate(
    amino_3 = factor(amino_3, levels = average_effect_amino_order),
    specie = factor(specie, levels = c("xenopus", "fish", "mouse", "human"))
  ) %>% 
  ggplot(aes(x=str_sub(codon, 1, 3), y=specie, fill = optimality)) +
  geom_tile() +
  scale_fill_viridis_c(option = "B",limits = c(-.03, .03), oob = scales::squish, breaks = c(-.02, 0, .02)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  facet_grid(~amino_3, scales="free_x", space = "free_x") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggsave("lasso-coefs-by-amino.pdf", height = 1.7, width = 15)
