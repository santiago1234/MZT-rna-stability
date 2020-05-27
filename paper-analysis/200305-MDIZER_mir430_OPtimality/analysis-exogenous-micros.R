library(tidyverse)
library(ggthemes)
library(scales)
library(ggforce)



quantifications <- read_csv("../../data/19-05-08-MammalianMircroRNAs/rna-seq-quants-mammlian-micro-rnas.csv")

micro_rnas <- read_csv("../../data/19-05-08-MammalianMircroRNAs/mircro_rnas_tidy.csv")
pls <- read_csv("../../results/19-08-19-OverlapFinal/results_data/pls_species.csv")


mapping <- function(x) {
  if (x > 1) return(">= 2 sites")
  if (x == 1) return("1 site")
  else return("0 sites")
}

# get the data

mir155 <- quantifications %>% 
  filter(
    sample_name %in% c("mrna_mir1_32hr", "mrna_mir1_12hr", "mrna_mock_32hr", "mrna_mock_12hr")
  ) %>% 
  left_join(filter(micro_rnas, `Representative miRNA` == "hsa-miR-155-5p")) %>% 
  replace_na(list(`Conserved sites total` = 0)) %>% 
  inner_join(pls) %>% 
  mutate(`Conserved sites total` = map_chr(`Conserved sites total`, mapping),
         `Conserved sites total` = factor(`Conserved sites total`, levels = c(
           "0 sites", "1 site", ">= 2 sites")),
         sample_name = str_replace(sample_name, "mrna_", ""),
         targets = `Conserved sites total` != "0 sites"
  )


## some duplicates
duplicados <- 
  mir155 %>% 
  count(gene_id, sample_name) %>% 
  filter(n > 1) %>% 
  pull(gene_id)

mir155_tidy <- mir155 %>% 
  filter(!gene_id %in% duplicados) %>% 
  select(gene_id, expression_level, PLS1, PLS2, sample_name, targets, `Conserved sites total`) %>% 
  spread(key = sample_name, value = expression_level) %>% 
  mutate(
    log2_FC = log2(mock_12hr / mir1_12hr)
  )

mir155_tidy %>% 
  ggplot(aes(x=PLS1, y=log2_FC, color=targets)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
  facet_wrap(~`Conserved sites total`) +
  coord_flip()

mir155_tidy %>% 
  filter(!is.na(log2_FC), !is.infinite(log2_FC)) %>% 
  filter(`Conserved sites total` == ">= 2 sites") -> x

lm(log2_FC ~ PLS1 + I(PLS1^2), data = x) %>% 
  summary
