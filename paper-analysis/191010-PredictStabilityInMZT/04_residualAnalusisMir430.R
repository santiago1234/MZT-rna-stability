library(tidyverse)
library(ggthemes)
library(ggforce)
library(gridExtra)

theme_set(theme_tufte(base_family = "Helvetica"))

mztdatos <- read_csv("results-data/mzt_predictionsResidualsLog2Fc.csv")


# classify the 3utr according to the kmer type ----------------------------

classify_utr <- function(utr_seq) {
  if (str_detect(utr_seq, "AGCACTTA")) return("8-mer")
  if (str_detect(utr_seq, "(AGCACTT|GCACTTA)")) return("7-mer")
  if (str_detect(utr_seq, "(AGCACT|GCACTT)")) return("6-mer")
  return("no miR seed")
  
}

mztdatos %>% 
  filter(!is.na(`3utr`)) %>% 
  mutate(
    miR430 = map_chr(`3utr`, classify_utr),
    miR430 = factor(miR430, levels = c("no miR seed", "6-mer", "7-mer", "8-mer"))
  ) -> mztdatos


## count the number of sites in the 3' UTR

mztdatos <- 
  mztdatos %>% 
  mutate(
    n_sites = str_count(`3utr`, "GCACTT"),
    n_sites = map_chr(n_sites, ~if_else(. > 1, ">1", as.character(.))),
    n_sites = factor(n_sites, levels = c("0", "1", ">1"))
  )


# master plot -------------------------------------------------------------
# I drop the variables number of sites and seed type so I can plot a facet grid
extended_mztdatos <- 
  mztdatos %>% 
  select(gene_id, log2FC, resid, miR430, n_sites, specie) %>% 
  gather(key = regulation_type, value = var, -gene_id, -log2FC, -resid, -specie) %>% 
  mutate(
    var = factor(var, levels = c("no miR seed", "6-mer", "7-mer", "8-mer", "0", "1", ">1"))
  )


# plot by site strength ---------------------------------------------------

site_strength <- extended_mztdatos %>% 
  filter(regulation_type == "miR430")


site_stats <- 
  site_strength %>% 
  group_by(specie, var) %>% 
  summarise(mediana = median(log2FC), n=n())

site_strength %>% 
  ggplot(
    aes(x=var, y=resid, color=var)
  ) +
  geom_sina(shape=16, size=1/10, alpha=.99) +
  facet_grid(~specie) +
  geom_errorbar(data = site_stats, aes(y=mediana, x=var, ymin=mediana, ymax=mediana), color="black", size=1/5) +
  geom_text(data = site_stats, aes(x=var, y=5, label=paste0("n=", n)), color="grey", size=1.5) +
  scale_color_manual(values = c("grey", "#009E73", "#009E73","#009E73")) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1.5),
    text = element_text(size=6),
    axis.ticks = element_line(size=1/5),
    axis.line = element_line(size = 1/4)
    ) +
  labs(
    x = "miR-430 k-mer seed",
    y = "residuals (observed - predicted)"
  )
ggsave("figures/residual-mir-res-kmerType.pdf", height = 2, width = 3)
aov(resid ~ var, data=site_strength) %>% summary()

# plot by number of sites -------------------------------------------------

nsite_strength <- extended_mztdatos %>% 
  filter(regulation_type == "n_sites")

nsite_stats <- 
  nsite_strength %>% 
  group_by(specie, var) %>% 
  summarise(mediana = median(log2FC), n=n())

nsite_strength %>% 
  ggplot(
    aes(x=var, y=resid, color=var)
  ) +
  geom_sina(shape=16, size=1/10, alpha=.99) +
  facet_grid(~specie) +
  geom_errorbar(data = nsite_stats, aes(y=mediana, x=var, ymin=mediana, ymax=mediana), color="black", size=1/5) +
  geom_text(data = nsite_stats, aes(x=var, y=5, label=paste0("n=", n)), color="grey", size=1.5) +
  scale_color_manual(values = c("grey", "#009E73", "#009E73")) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1.5),
    text = element_text(size=6),
    axis.ticks = element_line(size=1/5),
    axis.line = element_line(size = 1/4)
  ) +
  labs(
    x = "miR-430 k-mer seed",
    y = "residuals (observed - predicted)"
  )
ggsave("figures/residual-mir-nSites.pdf", height = 2, width = 2.5)



# ********************** --------------------------------------------------
# m6a plot ----------------------------------------------------------------
# ********************** --------------------------------------------------

pathways <- read_csv("../../results/19-02-24-OverlapPathwaysFig3/results_data/regulatory_pathways_matrix.csv")

# add the genes containing m6a

mztdatos <- 
  pathways %>% 
  select(Gene_ID, m6A) %>% 
  rename(gene_id = Gene_ID) %>% 
  inner_join(mztdatos)


resid_median <- mztdatos  %>% 
  group_by(m6A) %>% 
  summarise(mediana_r = median(resid), n=n())

mztdatos %>% 
  ggplot(aes(x=m6A, y=resid, color=m6A)) +
  geom_sina(shape=16, alpha=.99, size=1/20) +
  geom_errorbar(data = resid_median, aes(y=mediana_r, x=m6A, ymin=mediana_r, ymax=mediana_r), color="black", size=1/5) +
  geom_text(data=resid_median, aes(x=m6A, y=4, label=paste0("n=", n)), color="grey", size=1.5) +
  labs(
    x = NULL,
    y = "residuals (observed - predicted)"
  ) +
  scale_color_manual(values = c("grey", "#E69F00")) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1.5),
    text = element_text(size=6),
    axis.ticks = element_line(size=1/5),
    axis.line = element_line(size = 1/4)
  )

ggsave("figures/residual-m6a.pdf", height = 2, width = 1.5)  



