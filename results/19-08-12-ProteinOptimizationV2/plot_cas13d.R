library(tidyverse)
library(ggrepel)

res <- read_csv("results_data/deoptimizationCas13d_allregion_best_log.csv")

res %>% 
  ggplot(aes(x=iteration, y=fitness)) +
  geom_point() +
  geom_line(aes(y=ci_l)) +
  geom_line(aes(y=ci_u)) 



# load predictions for fish to show histrogram ----------------------------

fishpreds <- 
  read_csv("../19-08-06-UncertantyPredictionsModel/results_data/prediction_interavals_test_data.csv") %>% 
  filter(specie == "fish", datatype == "aamanitin polya")

## make histogram


# selections --------------------------------------------------------------

cas13d <- read_csv("results_data/cas13dS.csv") %>% 
  mutate(pos_y = c(10, 20, 30, 40))

cas13d %>% 
  ggplot(aes(x=prediction, y=name)) +
  geom_point() +
  coord_cartesian(xlim = c(-1, 1))


fishpreds %>% 
  ggplot(aes(x=median_prediction)) +
  geom_histogram(bins=100, color="grey", fill='grey') +
  geom_point(
    data = cas13d,
    aes(x=prediction, y=pos_y)
  ) +
  ggrepel::geom_text_repel(data = cas13d, aes(x=prediction, y=pos_y, label=name)) +
  ggthemes::theme_tufte()
ggsave("figures/cas13d_results.pdf")



# reporter sequences ------------------------------------------------------
# compare previous optimizatin --------------------------------------------

oldopt <- read_csv("results_data/old_optimization_with_new_predictions_synReps.csv") %>% 
  select(iteration, fitness, optimization, median_pred:ci_u)

oldopt %>% 
  ggplot(aes(x=fitness, y=median_pred, color=iteration)) +
  geom_point() +
  facet_grid(~optimization)

oldopt %>% 
  ggplot(aes(x=iteration, y=fitness, group=optimization, color=optimization)) +
  geom_point() +
  geom_point(aes(y=median_pred),shape=1) +
  geom_errorbar(aes(ymin=ci_l, ymax=ci_u), alpha=2/5)

ggsave("figures/optimization_path_old_and_new_synreps.pdf")
# reporters ---------------------------------------------------------------

reporters <- read_csv("results_data/reporters_predicitons.csv") %>% 
  select(-predicted_stability_human293t)
## load q's experimenta results

qs_res <- readxl::read_excel("../19-07-29-ExperimentalResultsOptimizationGfpMajo/results_data/Q_silent_Reporters.xlsx") %>% 
  rename(id_seq = X__1)

reporters <- inner_join(qs_res, reporters, by=c("id_seq"))

reporters %>% 
  ggplot(aes(x=median_pred, y=`mcherry/GFP`)) +
  geom_point() +
  geom_errorbarh(aes(xmin=ci_l, xmax=ci_u), alpha=1/2) +
  geom_text_repel(aes(label=id_seq)) +
  scale_y_log10()

ggsave("figures/synonimous_reporters.pdf", width = 3, height = 2)
