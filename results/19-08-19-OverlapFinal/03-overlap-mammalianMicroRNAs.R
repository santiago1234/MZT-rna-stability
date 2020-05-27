library(tidyverse)
library(ggthemes)
library(scales)
library(ggforce)
theme_set(theme_tufte(base_family = "Helvetica"))


quantifications <- read_csv("../../data/19-05-08-MammalianMircroRNAs/rna-seq-quants-mammlian-micro-rnas.csv")

micro_rnas <- read_csv("../../data/19-05-08-MammalianMircroRNAs/mircro_rnas_tidy.csv")
pls <- read_csv("results_data/pls_species.csv")


mapping <- function(x) {
  if (x > 1) return(">= 2 sites")
  if (x == 1) return("1 site")
  else return("0 sites")
}

# get the data

mir155 <- quantifications %>% 
  filter(
    sample_name %in% c("mrna_mir155_32hr", "mrna_mir155_12hr", "mrna_mock_32hr", "mrna_mock_12hr")
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


mir155 <- mir155 %>% 
  filter(!str_detect(sample_name, "12hr"), expression_level > 0)

# pvalue optimality effects in mrna level
lm(log(expression_level) ~ PLS2, data = mir155) %>% summary
# PLS2        0.224319 <2e-16 ***
mir155 <- 
  mir155 %>% 
  select(gene_id, expression_level, PLS1, PLS2, targets, sample_name) %>% 
  group_by(sample_name, targets) %>% 
  mutate(optimality = ntile(-PLS2, 4)) %>% 
  ungroup() %>% 
  mutate(micro = "miR-155")


mir1 <- quantifications %>% 
  filter(
    sample_name %in% c("mrna_mir1_12hr", "mrna_mir1_32hr", "mrna_mock_32hr", "mrna_mock_12hr")
  ) %>% 
  left_join(filter(micro_rnas, `Representative miRNA` != "hsa-miR-155-5p")) %>% 
  replace_na(list(`Conserved sites total` = 0)) %>% 
  inner_join(pls) %>% 
  mutate(`Conserved sites total` = map_chr(`Conserved sites total`, mapping),
         `Conserved sites total` = factor(`Conserved sites total`, levels = c(
           "0 sites", "1 site", ">= 2 sites")),
         sample_name = str_replace(sample_name, "mrna_", ""),
         targets = `Conserved sites total` != "0 sites"
  )


mir1 <- mir1 %>% 
  filter(!str_detect(sample_name, "12hr"), expression_level > 0)

mir1 <- 
  mir1 %>% 
  select(gene_id, expression_level, PLS1, PLS2, targets, sample_name) %>% 
  group_by(sample_name, targets) %>% 
  mutate(optimality = ntile(-PLS2, 4)) %>% 
  ungroup() %>% 
  mutate(micro = "miR-1")

# do not plot the mocj results
dtaplot <- bind_rows(mir1, mir155) %>%
  filter(!str_detect(sample_name, "mock"))

dtaplot %>% 
  filter(optimality %in% c(1, 4)) %>% 
  ggplot(aes(x=targets, y=expression_level, color=as.character(optimality))) +
  geom_rangeframe(sides = "l", color="black", alpha=2/3) +
  geom_sina(size=1/5, shape=16, alpha=.99) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_discrete(labels = c("no targets", "targets")) +
  facet_grid(.~micro) +
  scale_color_manual(values = c("red", "blue")) +
  labs(
    title = "mammalian microRNAs",
    y = "RNA level\nlog10 RPKM",
    x = NULL
  ) +
  theme(
    legend.position = "none"
  )

ggsave("figures/mammalian_micros.pdf", height = 2, width = 4)

dtaplot %>% 
  group_by(sample_name, micro, targets) %>% 
  nest() %>% 
  mutate(
    fit = map(data, ~lm(log(expression_level) ~ PLS2, data = .)),
    tf = map(fit, broom::tidy)
  ) %>% 
  unnest(tf) %>% 
  filter(term != "(Intercept)")


# minimalist plot ---------------------------------------------------------

## get summary stats

summary_expression <- 
  dtaplot %>% 
  mutate(optimality = str_c("q", optimality)) %>%
  group_by(micro, targets, optimality) %>% 
  summarise(
    mediana = median(expression_level),
    n = n()
  )


dtaplot %>% 
  mutate(optimality = str_c("q", optimality)) %>%
  ggplot(aes(x=optimality, y=expression_level, color=optimality)) +
  geom_sina(shape=16, alpha=.99, size=1/10) +
  geom_rangeframe(sides = "l", color="black", size=1/5) +
  geom_errorbar(
    data = summary_expression,
    aes(ymin=mediana, ymax=mediana, y=mediana, x=optimality),
    color="black",
    size=1/5
  ) +
  geom_text(
    data = summary_expression,
    aes(x = optimality, y = 10**3.5, label = paste0("n=", n)),
    color = "grey",
    size=1
  ) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values = c("#ca0020", "#f4a582", "#92c5de", "#0571b0")) +
  facet_grid(~micro + targets) +
  labs(
    title = "mammalian microRNAs",
    y = "RNA level\nlog10 RPKM",
    x = "optimality level"
  ) +
  theme(legend.position = "none", axis.ticks = element_line(size = 1/5))
ggsave("figures/mammalian_micros_minimal.pdf", height = 2, width = 4)
