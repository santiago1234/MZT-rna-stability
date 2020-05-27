library(tidyverse)
library(gam)
library(tidybayes)

theme_set(tidybayes::theme_tidybayes())

dicer <- "/Volumes/projects/smedina/projectos/181029-Paper/data/rna-time-course/dicer-data/cufflinks_fpkm_all.csv" %>% 
  read_csv()


pathways <- "../../results/19-02-24-OverlapPathwaysFig3/results_data/regulatory_pathways_matrix.csv" %>% 
  read_csv() %>% 
  filter(is_maternal) %>% 
  select(Gene_ID, MiR430)

optimality <- read_csv("../../results/19-08-19-OverlapFinal/results_data/pls_species.csv") %>% 
  rename(Gene_ID = gene_id) %>% 
  select(-specie)

pathways <- inner_join(pathways, optimality)

datos <- dicer %>% 
  mutate(
    log2_WT_Dicer = log2(s_WT_6h / (s_MZdicer_6h + .0001))
  ) %>% 
  select(GeneID, log2_WT_Dicer) %>% 
  rename(Gene_ID = GeneID) %>% 
  inner_join(pathways)
  

# get genes with micro RNA

mir_datos <- 
  datos %>% 
  filter(MiR430 > 0, log2_WT_Dicer < 2.5)

fit <- gam(log2_WT_Dicer ~  s(PLS1, 2) + s(PLS2, 5) + MiR430, data = mir_datos)

# fold change does not depend on the codin sequence
fit0 <- gam(log2_WT_Dicer ~  MiR430, data = mir_datos)

# signficance
summary(fit)
anova(fit0, fit, test="F")

par(mfrow=c(1, 3))
plot(fit, se=TRUE,col="blue")


mir_datos %>% 
  ggplot(aes(x = PLS1, y = log2_WT_Dicer)) +
  geom_point() +
  geom_line(aes(y = prediction), color = "red")


# predictions on new data -------------------------------------------------

x_new <- tibble( 
  PLS1 = seq(from = -5.5 ,to = 7, length.out = 30),
  PLS2 = seq(from = -5.5, to = 7, length.out = 30),
  id = 1:30,
  MiR430 = 1
)

x_new$prediction <- predict(fit, newdata = x_new)

x_new %>% 
  ggplot(aes(x = PLS1, y = prediction)) +
  geom_point(
    data = mir_datos,
    aes(x = PLS1, y = log2_WT_Dicer),
    color = "forestgreen"
  ) +
  geom_hline(yintercept = 0, color = "grey20", linetype = 2) +
  geom_line() +
  coord_cartesian(ylim = c(-3.5, 1))

# bootstrap analysis ------------------------------------------------------

boostrap <- function() {
  d_b <- sample_frac(mir_datos, replace = T)
  fit <- gam(log2_WT_Dicer ~  s(PLS1, 2) + s(PLS2, 5) + MiR430, data = d_b)
  
  x_new$prediction <- predict(fit, newdata = x_new)
  
  x_new
  
}

set.seed(42)
results <- tibble(
  bid = as.character(1:100)
) %>% 
  mutate(b = map(bid, function(i) boostrap())) %>% 
  unnest(b)


results %>% 
  ggplot(aes(x = PLS1, y = prediction)) +
  geom_point(
    data = sample_frac(mir_datos, size = .7),
    aes(x = PLS1, y = log2_WT_Dicer),
    color = "grey",
    shape = 16,
    alpha = .7
  ) +
  stat_lineribbon() +
  geom_hline(yintercept = 0, color = "grey20", linetype = 2) +
  scale_y_continuous(limits = c(-3, 1.5)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_brewer("Confidence\n", palette = 2) +
  labs(
    y = "log2 (WT / Dicer)",
    x = "Codon optimality level",
    subtitle = 'p = 0.0031',
    title = 'Maternal genes with miR430 seed'
  )
ggsave('boostrap-res.pdf', height = 3, width = 3.5)


# boxplot -----------------------------------------------------------------

category <- function(x) {
  if (between(x, -Inf, -3)) return('non-optimal')
  if (between(x, -.5, -.3)) return('neutral')
  if (between(x, 3, Inf)) return('optimal')
  return('out')
  
}

mir_datos_boxplot <- 
  mir_datos %>% 
  mutate(
    opt = map_chr(PLS1, category)
  ) %>% 
  filter(opt != 'out') %>% 
  mutate(opt = factor(opt, levels = c('non-optimal', 'neutral', 'optimal')))

# significance using an annova test
fb1 <- lm(log2_WT_Dicer ~ opt, mir_datos_boxplot)
fb0 <- lm(log2_WT_Dicer ~ 1, mir_datos_boxplot)
anova(fb0, fb1)


mir_datos_boxplot %>% 
  ggplot(aes(x = opt, y = log2_WT_Dicer)) +
  geom_boxplot(aes(fill = opt), outlier.shape = NA, alpha = .9, size = .3) +
  coord_cartesian(ylim = c(-3.5, 1)) +
  scale_fill_manual(values = c('blue', 'grey', 'red')) +
  geom_hline(yintercept = 0, size=1/5, linetype=3) +
  labs(
    y = "log2 WT / Dicer",
    x = NULL,
    title = "genes with miR-430 seed",
    subtitle = 'p = 0.044'
  ) +
  theme(legend.position = 'none', axis.text.x = element_text(angle = 20, hjust = 1))
ggsave('boxplot.pdf', height = 3, width = 2.5)

