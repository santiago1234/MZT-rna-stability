library(tidyverse)
library(ggthemes)

theme_set(theme_tufte())

time_course <- read_csv("../../data/19-06-12-XenopusData/time_course_xenopus.csv")

time_course <-
  time_course %>%
  rename(Gene_ID = ensembl_gene_id) %>%
  filter(
    time < 13,
    sample_condition == "wt_polya" # PARAMETER
  ) %>% 
  rename(gene_id = Gene_ID) %>% 
  select(-jgi_id, -external_gene_name, -sample_condition)


top_expressed_genes <- 
  filter(time_course, time == 2) %>% 
  arrange(-expression_quantification) %>% 
  slice(1:6000) %>% 
  pull(gene_id)


time_course <- 
  filter(time_course, gene_id %in% top_expressed_genes)


# estimate the obsertved decay rate] --------------------------------------

observed_decay <- 
  time_course %>% 
  group_by(gene_id) %>% 
  nest() %>% 
  mutate(
    fit = map(data, ~lm(log(expression_quantification + 0.001) ~ time, data = .)),
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
  filter(str_detect(gene_id, "ENSXETG.*")) %>% 
  rename(b_predicted = predicted_stability)


results <- inner_join(observed_decay, predictions_mzt)

## put the coeficient in the same scale
fit_to_scale <- lm(b_observed ~ b_predicted, data = results)

results$b_predicted <- predict(fit_to_scale)

results %>% 
  ggplot(aes(x=b_predicted, y=b_observed)) +
  geom_point() +
  ggpubr::stat_cor()



rna_level_predictions <- 
  time_course %>% 
  select(time) %>% 
  unique() %>% 
  crossing(results) %>% 
  select(-b_observed) %>% 
  mutate(predicted_RNA = exp(A0 + b_predicted * time)) %>% 
  inner_join(select(time_course, gene_id, expression_quantification, time))


rna_level_predictions %>% 
  ggplot(aes(x=predicted_RNA, y=expression_quantification)) +
  geom_point(shape=".") +
  scale_y_log10() +
  scale_x_log10() +
  geom_abline() +
  ggpubr::stat_cor() +
  facet_grid(~time)

# herarchical clustering for sorting genes --------------------------------

expression <- rna_level_predictions %>% 
  select(gene_id, time, expression_quantification) %>% 
  mutate(expression_quantification = log(expression_quantification+0.001)) %>% 
  spread(key=time, value=expression_quantification) %>% 
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
  ggplot(aes(x=time, y=gene_id, fill=log(expression_quantification+0.001) - mean(log(expression_quantification+0.001)))) +
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
    #time = str_c(time, " h")
  )

expression_plot_data %>% 
  ggplot(aes(x=time, y=gene_id, fill=rna_level)) +
  geom_tile() +
  facet_grid(~data_From, scales = "free_y", space="free_y") +
  scale_fill_viridis_c(option = "A", limits=c(8, 18), oob=scales::squish) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 80, hjust = 1)
  )

