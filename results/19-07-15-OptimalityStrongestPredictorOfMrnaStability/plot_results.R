library(tidyverse)
library(ggthemes)
library(ggridges)
library(gridExtra)

theme_set(theme_tufte(base_family = "Helvetica"))
# plot regulatory pathway -------------------------------------------------

dset <- read_csv("results_data/pathways_mRNAstability.csv")
dset <- dset[!apply(dset, 1, function(x) any(is.na(x))), ,]

# example mouse -----------------------------------------------------------

get_genes_in_pathways <- function(grpspec) {
  
  m6a <- filter(grpspec, m6A != 0) %>% 
    mutate(pathway = "m6A")
  
  optimals <- arrange(grpspec, PLS_1) %>% 
    slice(1:500) %>% 
    mutate(pathway = "optimal")
  
  nonoptimals <- arrange(grpspec, -PLS_1) %>% 
    slice(1:500) %>% 
    mutate(pathway = "non-optimal")
  
  mir <- filter(grpspec, microRNAsites != 0) %>% 
    mutate(pathway = "microRNA")
  
  all <- grpspec %>% 
    mutate(pathway = "all genes")
  
  bind_rows(m6a, optimals, nonoptimals, mir, all)
  
}

pathways <- 
  dset %>% 
  group_by(specie) %>% 
  nest() %>% 
  mutate(pat = map(data, get_genes_in_pathways)) %>% 
  unnest(pat)

p1 <- pathways %>% 
  mutate(pathway = factor(pathway, levels = c("non-optimal", "m6A", "microRNA", "all genes","optimal"))) %>% 
  ggplot(aes(y=stability, x = pathway, fill = pathway)) +
  geom_boxplot(size=1/5, alpha = 3/4, outlier.shape = NA) +
  facet_grid(~specie) +
  scale_fill_manual(values = c("blue", "goldenrod2", "forestgreen" , "grey", "red")) +
  coord_flip(ylim = c(-3, 3)) +
  theme(legend.position = "none") +
  labs(
    y = "mRNA stability (scaled)",
    x = NULL
  )

# get number of genes in each pathway -------------------------------------


mir <- pathways %>% 
  group_by(specie) %>% 
  summarise(fraction = mean(microRNAsites > 0)) %>% 
  mutate(model = "microRNA")

m6a <- pathways %>% 
  group_by(specie) %>% 
  summarise(fraction = mean(m6A > 0)) %>% 
  mutate(model = "m6A")

optimality <- pathways %>% 
  select(specie) %>% 
  unique() %>% 
  mutate(fraction = 1) %>% 
  mutate(model = "optimality")

mdlw <- read_csv("results_data/mdl_weights.csv")

mdlw <- 
  bind_rows(mir, m6a, optimality) %>% 
  inner_join(mdlw)

p2 <- mdlw %>% 
  ggplot(aes(x=fraction, y=weights, color=model)) +
  geom_point(alpha=1/2, shape=16) +
  ggrepel::geom_text_repel(aes(label=model, color=model)) +
  facet_grid(~specie) +
  scale_color_manual(values = c("goldenrod2", "forestgreen", "red")) +
  geom_rangeframe(color = "black", alpha=4/5) +
  theme(legend.position = "none") +
  labs(y = "model weights\n(Bayesian bootstrap)", x="fraction of coding genes with potential regulation")

pdf("model_comparison.pdf", width = 6, height = 3)
grid.arrange(p1, p2)
dev.off()
