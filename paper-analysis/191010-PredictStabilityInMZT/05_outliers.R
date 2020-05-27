library(tidyverse)

datos <- read_csv("results-data/mzt_predictionsResidualsLog2Fc.csv")

datos <- 
  datos %>% 
  group_by(specie) %>% 
  mutate(
    outlier = abs(resid) > IQR(resid)*1.5
  ) %>% 
  ungroup()

datos %>% 
  ggplot(aes(x = predicted, y = log2FC, color = outlier)) +
  geom_point() +
  facet_grid(~specie)



# get the orthologs, which genes are outliers and orthologs? --------------

ortolog <- read_csv("../../results/19-08-28-Fig04EvoAnalysis/results_data/orthologs_fishMaternal_to_Xen.csv") %>% 
  filter(!is.na(ortolog_xen))


tmp <- datos %>% 
  select(gene_id, specie, predicted_stability, outlier, log2FC)

colnames(tmp) <- colnames(tmp) %>% paste0("_fish")

tmp2 <- rename(ortolog, gene_id_fish = ortolog_fish) %>% 
  inner_join(tmp)

colnames(tmp) <- str_replace(colnames(tmp), "_fish", "_xenopus")

tmp2 <- inner_join(rename(tmp2, gene_id_xenopus = ortolog_xen), tmp)

tmp2 %>% 
  mutate(is_outlier_in_both = outlier_fish & outlier_xenopus) %>% 
  ggplot(aes(predicted_stability_fish, log2FC_fish)) +
  geom_point(color="grey", alpha=.99, shape=16) +
  geom_point(data = function(x) filter(x, is_outlier_in_both), color="black") +
  ggrepel::geom_text_repel(data = function(x) filter(x, is_outlier_in_both), aes(label=gene_id_fish), size=3) +
  labs(
    x = "predicted decay fish"
  ) +
  theme_minimal()
ggsave("tmp.pdf")
