library(tidyverse)
library(ggthemes)

theme_set(theme_tufte(base_family = "Helvetica"))

## parameters for analysis
mrna_seq <- "wt_ribo"

timecourse <- read_csv("../../data/19-01-09-RNAseqProfilesFish/rna-seq-profiles/RNAseq_tidytimecourse.csv") %>% 
  rename(gene_id = Gene_ID)

timecourse <- 
  timecourse %>% 
  filter(sample_condition == mrna_seq, is_maternal) 

observed_decay <- 
  timecourse %>% 
  group_by(gene_id) %>% 
  nest() %>% 
  mutate(
    fit = map(data, ~lm(log(TPM + 0.001) ~ time, data = .)),
    tidyfit = map(fit, broom::tidy)
  )


observed_decay <- 
  observed_decay %>% 
  select(-data, -fit) %>% 
  unnest(cols = tidyfit) %>% 
  ungroup() %>% 
  select(-std.error, -statistic, -p.value) %>% 
  spread(key = term, value = estimate) %>% 
  rename(A0 = `(Intercept)`, b_observed = time)


predictions_mzt <- read_csv("../191010-PredictStabilityInMZT/results-data/mzt_data_predictionsAndObservedStability.csv") %>% 
  select(gene_id, predicted_stability) %>% 
  filter(str_detect(gene_id, "ENSDAR.*")) %>% 
  rename(b_predicted = predicted_stability)


results <- inner_join(observed_decay, predictions_mzt)

## put the coeficient in the same scale
fit_to_scale <- lm(b_observed ~ b_predicted, data = results)

results$b_predicted <- predict(fit_to_scale)

results %>% 
  ggplot(aes(x=b_predicted, y=b_observed)) +
  geom_point() +
  ggpubr::stat_cor()



# make the plot -----------------------------------------------------------



rna_level_predictions <- 
  timecourse %>% 
  select(time) %>% 
  unique() %>% 
  crossing(results) %>% 
  select(-b_observed) %>% 
  mutate(predicted_RNA = exp(A0 + b_predicted * time)) %>% 
  inner_join(select(timecourse, gene_id, TPM, time))


rna_level_predictions %>% 
  ggplot(aes(x=predicted_RNA, y=TPM)) +
  geom_point(shape=".") +
  scale_y_log10() +
  scale_x_log10() +
  geom_abline() +
  ggpubr::stat_cor() +
  facet_grid(~time)

cor.test(rna_level_predictions$TPM, rna_level_predictions$predicted_RNA)

# herarchical clustering for sorting genes --------------------------------

expression <- rna_level_predictions %>% 
  select(gene_id, time, TPM) %>% 
  mutate(logTPM = log(TPM+0.001)) %>% 
  select(-TPM) %>% 
  spread(key=time, value=logTPM) %>% 
  as.data.frame()


row.names(expression) <- expression$gene_id
expression$gene_id <- NULL
expression <- as.matrix(expression)
expression <- scale(expression)

cluster_data <- hclust(d = dist(x = expression))

# make the rows to show as clustered --------------------------------------

rna_level_predictions <- rna_level_predictions %>% 
  mutate(gene_id = factor(gene_id, levels = cluster_data$labels[cluster_data$order]))



rna_level_predictions %>% 
  ggplot(aes(x=time, y=gene_id, fill=log(TPM+0.001) - mean(log(TPM+0.001)))) +
  geom_tile() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  scale_fill_viridis_c()

# I will log transform the RNA level
expression_plot_data <- 
  rna_level_predictions %>% 
  select(-A0, -b_predicted) %>% 
  gather(key = data_From, value = rna_level, -time, -gene_id) %>% 
  mutate(
    rna_level = log(rna_level),
    time = str_c(time, " h")
  )


# get mir430 gene information ---------------------------------------------

expression_plot_data <- 
  read_csv("~/rna_stability/paper-analysis/191010-PredictStabilityInMZT/results-data/mzt_data_predictionsAndObservedStability.csv")  %>% 
  select(gene_id, `3utr`) %>% 
  mutate(mir430 = str_detect(`3utr`, "GCACTT")) %>% 
  select(-`3utr`) %>% 
  filter(gene_id %in% expression_plot_data$gene_id) %>% 
  mutate(gene_id = factor(gene_id, levels = levels(expression_plot_data$gene_id))) %>% 
  inner_join(expression_plot_data)


# make the plot -----------------------------------------------------------

expression_plot_data %>% 
  mutate(data_From = factor(data_From, levels = c("TPM", "predicted_RNA"))) %>% 
  ggplot(aes(x=time, y=gene_id, fill=rna_level)) +
  geom_tile() +
  facet_grid(mir430~data_From, scales = "free_y", space="free_y") +
  scale_fill_viridis_c(option = "A", limits=c(1.5, 5), oob=scales::squish) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 80, hjust = 1)
  )

ggsave("heatmap.pdf", height = 8, width = 4)


expression_plot_data %>% 
  filter(mir430, time %in% c("4 h", "5 h", "6 h", "7 h")) %>%
  ggplot(aes(x=time, rna_level, fill=data_From)) +
  geom_boxplot(outlier.shape = NA, size=1/3) +
  coord_cartesian(ylim = c(-1, 7)) +
  geom_rangeframe() +
  scale_fill_colorblind() +
  labs(
    x = "time post feritlization (hrs)",
    y = "log RNA level (TPM)",
    title = "genes with miR-430 seed in 3' UTR"
  )
ggsave("boxplot.pdf", height = 2, width = 4)


