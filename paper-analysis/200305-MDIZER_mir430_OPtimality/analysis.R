library(tidyverse)
library(forcats)
library(optimalcodonR)
library(ggthemes)



fish_genes <- 
  bind_rows(training, testing) %>% 
  filter(specie == "fish") %>% 
  select(gene_id, coding) %>% 
  unique() %>% 
  filter(!duplicated(coding)) %>% 
  rename(Gene_ID = gene_id)

fish_predictor <- predict_stability("fish")

fish_genes$predicted_stability <- fish_predictor(fish_genes$coding)

# get dizer mutant data ---------------------------------------------------

dicer <- "/Volumes/projects/smedina/projectos/181029-Paper/data/rna-time-course/dicer-data/cufflinks_fpkm_all.csv" %>% 
  read_csv()

pathways <- "../../results/19-02-24-OverlapPathwaysFig3/results_data/regulatory_pathways_matrix.csv" %>% 
  read_csv() %>% 
  filter(is_maternal) %>% 
  select(Gene_ID, PLS1, MiR430)

#mir430codes <- read_csv("../../data/19-01-17-EMBO2016DATA/datasets/Groups/mir430_genes_targets.csv") %>% 
#  rename(Gene_ID = gene_id) %>% 
#  select(Gene_ID, MicroArray)

datos <- dicer %>% 
  mutate(
    log2_WT_Dicer = log2(s_WT_6h / (s_MZdicer_6h + .0001))
  ) %>% 
  select(GeneID, log2_WT_Dicer) %>% 
  rename(Gene_ID = GeneID) %>% 
  inner_join(pathways) %>% 
  inner_join(fish_genes)

#datos <- full_join(datos, mir430codes)

datos %>% 
  ggplot(aes(x=log2_WT_Dicer)) +
  geom_density(aes(color=MiR430 > 0))


datos %>% 
  filter(MiR430 > 0, log2_WT_Dicer < 2.5) %>% 
  ggplot(aes(x = PLS1, y = log2_WT_Dicer)) +
  geom_point(shape = 16, alpha = .7, color = "grey") +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color="purple") +
  theme_tufte() +
  scale_color_viridis_d() +
  geom_vline(xintercept = median(datos$PLS1), size=1/5, linetype=3) +
  geom_hline(yintercept = 0, size=1/5, linetype=3) +
  geom_rangeframe(color = "black") +
  coord_flip(ylim = c(-4, .8)) +
  labs(
    y = "log2 WT / Dicer",
    title = "neutral genes are more sensitive to miR-430",
    subtitle = "miR-430 targets p-value = 6.54e-12",
    x = "Codon optimality level"
  )

ggsave("fig.pdf", height = 2.5, width = 4)


lm(log2_WT_Dicer ~ PLS1 + I(PLS1^2), data = datos) %>% 
  summary()

# color code by the number of sites ---------------------------------------

datos %>% 
  filter(MiR430 > 0, log2_WT_Dicer < 2.5) %>% 
  ggplot(aes(x = predicted_stability, group=MiR430, color=MiR430)) +
  geom_density()


datos %>% 
  filter(MiR430 > 0, log2_WT_Dicer < 2.5) %>% 
  ggplot(aes(x = predicted_stability, group=MiR430, color=MiR430)) +
  geom_density()

