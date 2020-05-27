library(tidyverse)
library(ggthemes)
library(ggExtra)
library(gridExtra)
library(ggforce)

theme_set(theme_tufte(base_family = "Helvetica"))
# load the data -----------------------------------------------------------

fold <- read_csv("../../data/19-02-05-FoldChangeData/data/log2FC_earlyVSlate_tidytimecourse.csv") %>% 
  filter(sample_condition == "wt_ribo", time == 6, is_maternal) %>% 
  select(-time, -sample_condition, -is_maternal)

# pathways
pats <- read_csv("../19-02-24-OverlapPathwaysFig3/results_data/regulatory_pathways_matrix.csv") %>% 
  select(Gene_ID, m6A, MiR430) %>% 
  mutate(MiR430 = MiR430)

# predicted stability
predictions <- read_csv("results_data/stability_predictions_fish.csv") %>% 
  rename(Gene_ID = gene_id)


prportion_optimal <- function(seq) {
  # compute the proportion of optimal codons embo1026
  # this is the definition of the embo paper
  codonr::count_codons(seq) %>% 
    full_join(codonr::optimality_code_embo2016) %>% 
    replace_na(list(n = 0)) %>% 
    group_by(optimality) %>% 
    summarise(total = sum(n)) %>% 
    spread(key = optimality, value = total) %>% 
    transmute(percent = log2(optimal / (`non-optimal` + 1))) %>% 
    pull(percent)
}


predictions <- 
  predictions %>% 
  mutate(prop_opt = map_dbl(coding, prportion_optimal))


# tidy the data -----------------------------------------------------------

datum <- inner_join(fold, predictions) %>% 
  inner_join(pats) %>% 
  select(-coding)

## put the predictions in the same scale with a linear model

fit <- lm(log2FC ~ predicted_stability, data = datum)

datum <- 
  datum %>% 
  mutate(
    predicted = predict(fit),
    residual = resid(fit),
    optimality_p = ntile(prop_opt, n = 4)
  )



# define gene grps: optimal, non-optimal, miR-430, neutral ----------------

mirg <- tibble(Gene_ID = filter(datum, MiR430 > 0)$Gene_ID, grp = "miR-430")

optimal <- 
  filter(datum, !Gene_ID %in% mirg$Gene_ID, optimality_p == 4) %>% 
  select(Gene_ID) %>% 
  mutate(grp = "optimal")

nonoptimal <- 
  filter(datum, !Gene_ID %in% mirg$Gene_ID, optimality_p == 1) %>% 
  select(Gene_ID) %>% 
  mutate(grp = "non-optimal")

grps <- bind_rows(mirg, optimal, nonoptimal)

datum <- 
  full_join(grps, datum) %>% 
  replace_na(list(grp = "neutral"))


# plot scatter figure 1 ---------------------------------------------------



p <- datum %>% 
  sample_frac(1) %>% 
  filter(grp != "neutral") %>% 
  ggplot(aes(x=predicted, y=log2FC, color = grp)) +
  geom_point(shape = 16, alpha = .99, size = 1) +
  scale_color_manual(values = c("#4daf4a", "#377eb8", "#e41a1c")) +
  ggpubr::stat_cor(color = "black") +
  labs(
    x = "predicted mRNA stability",
    y = "log2 fold change\n(6hrs \ 3hrs)"
  )


pdf("preds_mzt.pdf", width = 6, height = 4)
ggMarginal(p, groupColour = TRUE)
dev.off()

# residual plot -----------------------------------------------------------

rm6a <- datum %>% 
  filter(m6A) %>% 
  mutate(grp = "m6A")

rm1 <- datum %>% 
  filter(MiR430 == 1) %>% 
  mutate(grp = "1 miR-430 site")

rm2 <- datum %>% 
  filter(MiR430 == 2) %>% 
  mutate(grp = "2 miR-430 sites")

rm3 <- datum %>% 
  filter(MiR430 > 2) %>% 
  mutate(grp = ">2 miR-430 sites")

rall <- datum %>% 
  mutate(grp = "all genes")

bind_rows(rm1, rm2, rm3, rm6a, rall) %>% 
  mutate(
    grp = factor(grp, levels = c("all genes", "m6A", ">2 miR-430 sites", "2 miR-430 sites", "1 miR-430 site"))
  ) %>% 
  ggplot(aes(x=grp, y=residual,color = grp)) +
  geom_sina(shape=16, size=1/3, alpha=.99) +
  scale_color_manual(values = c("black", "goldenrod3", "forestgreen", "forestgreen", "forestgreen")) +
  geom_hline(yintercept = c(-2.5, 0, 2.5), color = "white", size=1/2) +
  scale_y_continuous(breaks = c(-2.5, 0, 2.5)) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90),
    axis.ticks.x = element_blank()
  ) +
  labs(
    x = NULL,
    y = "observed (mzt) - predicted \n(residuals)"
  )

ggsave("residual_plots.pdf", width = 3, height = 4)

bind_rows(rm1, rm2, rm3, rm6a, rall) %>% 
  mutate(
    grp = factor(grp, levels = c("all genes", "m6A", ">2 miR-430 sites", "2 miR-430 sites", "1 miR-430 site"))
  ) %>% 
  ggplot(aes(x=grp, y=residual,fill = grp)) +
  #geom_sina(shape=16, size=1/3, alpha=.99) +
  geom_boxplot(outlier.shape = NA, size=1/3, alpha=.85) +
  scale_fill_manual(values = c("black", "goldenrod3", "forestgreen", "forestgreen", "forestgreen")) +
  scale_y_continuous(breaks = c(-2.5, 0, 2.5)) +
  geom_hline(yintercept = 0, color="white") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.ticks.x = element_blank()
  ) +
  labs(
    x = NULL,
    y = "observed (mzt) - predicted \n(residuals)"
  )
ggsave("residual_plots2.pdf", width = 3, height = 4)


############################################################################
############################################################################
# minimalist plot ---------------------------------------------------------
############################################################################
############################################################################

datum %>% 
  group_by(grp) %>% 
  mutate(lfc = median(log2FC)) %>% 
  ungroup() %>% 
  sample_frac(1) %>% 
  filter(grp != "neutral") %>% 
  mutate(grp = factor(grp, levels = c("optimal", "non-optimal", "miR-430"))) %>% 
  ggplot(aes(x=predicted, y=log2FC, color = grp)) +
  geom_point(alpha = .79, size = 1, shape=16) +
  scale_color_manual(values = c("black", "grey", "steelblue")) +
  coord_cartesian(ylim = c(-6, 3)) +
  theme(legend.position = "none")


ggsave("test_minaml.pdf", width = 5, height = 2.5)
 

## residual plot


rm6a <- datum %>% 
  filter(m6A) %>% 
  mutate(grp = "m6A", meta = "m6A")

rm1 <- datum %>% 
  filter(MiR430 == 1) %>% 
  mutate(grp = "1 site", meta = "miR-430 targets")

rm2 <- datum %>% 
  filter(MiR430 == 2) %>% 
  mutate(grp = "2 sites", meta = "miR-430 targets")

rm3 <- datum %>% 
  filter(MiR430 > 2) %>% 
  mutate(grp = ">2 sites" , meta = "miR-430 targets")
  

rall <- datum %>% 
  mutate(grp = "all genes", meta = "all")


res_p_data <- 
  bind_rows(rm1, rm2, rm3, rm6a, rall) %>% 
  mutate(
    grp = factor(grp, levels = c("all genes", "m6A", ">2 sites", "2 sites", "1 site"))
  )

res_p_data_summary <- res_p_data %>% 
  group_by(meta, grp) %>% 
  summarise(medianval = median(residual))
  
res_p_data %>% 
  ggplot(aes(x=grp, y=residual)) +
  geom_sina(shape=16, size=1/3, alpha=.9) +
  geom_rangeframe(size=1/5) +
  geom_point(
    data = res_p_data_summary, aes(y = medianval),
    color="red",
    size = 2,
    alpha = .99,
    shape = 16
  ) +
  facet_grid(~meta, scales = "free_x", space = "free_x")
ggsave("residual_minimal.pdf", width = 3, height = 3)
