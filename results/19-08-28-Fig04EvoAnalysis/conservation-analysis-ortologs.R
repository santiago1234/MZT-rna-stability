# protocol
# 1. take ortholog genes
# 2 take the top unstable and then take the intersection
# 3 Divide in 4 grps
# plot optimality level in each of the groups

library(tidyverse)
library(ggthemes)
library(gridExtra)

theme_set(theme_tufte(base_family = "Helvetica"))

# load the orthologs ------------------------------------------------------
# and select the top unstable
ortologs <- read_csv("results_data/orthologs_fishMaternal_to_Xen.csv") %>% 
  filter(!is.na(ortolog_xen))
stability_fc <- read_csv("results_data/fold_change_late_vs_early.csv")
optimality <- read_csv("results_data/optimality_information.csv") %>% 
  select(gene_id, PLS1, optimality_ratio, miR430) %>% 
  filter(!is.na(miR430))

convined_datos <- 
  inner_join(ortologs, rename(stability_fc, ortolog_fish = gene_id)) %>% 
  select(-specie) %>% 
  rename(log2FC_fish = log2FC) %>% 
  inner_join(rename(stability_fc, ortolog_xen = gene_id)) %>% 
  rename(log2FC_xen = log2FC) %>% 
  select(-specie)

# add the optimality to the orthologs table

rename(
  optimality,
  PLS1_fish = PLS1,
  optimality_ratio_fish = optimality_ratio,
  ortolog_fish = gene_id,
  miR430_fish = miR430
) %>% 
  inner_join(convined_datos, by = c("ortolog_fish")) -> convined_datos


rename(
  optimality,
  PLS1_xen = PLS1,
  optimality_ratio_xen = optimality_ratio,
  ortolog_xen = gene_id,
  miR430_xen = miR430
) %>% 
  inner_join(convined_datos, by = c("ortolog_xen")) -> convined_datos


# rank by stability -------------------------------------------------------

convined_datos <- 
  convined_datos %>% 
  mutate(
    stability_rank_fish = ntile(log2FC_fish, 3),
    stability_rank_xen = ntile(log2FC_xen, 3),
  )


# plot the optimality and the log2FC --------------------------------------

p_opt <- convined_datos %>% 
  ggplot(aes(optimality_ratio_fish, optimality_ratio_xen)) +
  geom_point(size=1/3, alpha=.99, shape=16) +
  labs(
    x = "optimality log2 ratio\nfish",
    y = "optimality log2 ratio\nXenopus"
  ) +
  geom_rangeframe(size=1/5) +
  ggpubr::stat_cor(size=2)

p_fc <- convined_datos %>% 
  ggplot(aes(log2FC_fish, log2FC_xen)) +
  geom_point(size=1/3, alpha=.99, shape=16) +
  labs(
    x = "mRNA stability log2 TPM\nfish",
    y = "mRNA stability log2 TPM\nXenopus"
  ) +
  geom_rangeframe(size=1/5) +
  ggpubr::stat_cor(size=2)
  

pdf("figures/05-orthologs-foldchange-optimality.pdf", height = 3.5, width = 3)
grid.arrange(p_opt, p_fc)
dev.off()
# make the 4 grps bassed on the presence of miR-430 -----------------------

which_grp <- function(mir_fish, mir_xen) {
  # assume mir_fish and mir_xen are integers given the number of miR-430 seeds
  # in the 3' UTR
  if (mir_fish > 0 & mir_xen > 0) return("target in both")
  if (mir_fish > 0 & mir_xen == 0) return("only fish")
  if (mir_fish == 0 & mir_xen == 0) return("no target")
  return("only xenopus")
}

convined_datos <- 
  convined_datos %>% 
  mutate(
    group = map2_chr(.x = miR430_fish, .y = miR430_xen, which_grp)
  ) 

# make the plot  ----------------------------------------------------------

datos_plot <- convined_datos %>% 
  filter(stability_rank_xen == stability_rank_fish, stability_rank_xen == 1)

datos_plot %>% 
  ggplot(aes(x=reorder(group, optimality_ratio_fish, mean), y=optimality_ratio_fish, fill=group)) +
  geom_boxplot(outlier.shape = NA, size=1/3) +
  geom_rangeframe(size=1/4) +
  scale_fill_viridis_d(option = "A") +
  labs(
    x = "ortholog gene miR-430 target type",
    y = "codon optimality ratio\nbazzini et al 2016",
    title = "unstable genes in both specie"
  ) +
  theme_tufte(base_family = "Helvetica") +
  theme(legend.position = "none")

ggsave("figures/05-analysis-orthologs-conserved-miR430-type.pdf", height = 3.5, width = 3.5)
