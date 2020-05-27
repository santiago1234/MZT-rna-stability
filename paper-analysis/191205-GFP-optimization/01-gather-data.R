library(tidyverse)
library(optimalcodonR)

# experimental results majo -----------------------------------------------

experiment <- readxl::read_excel("../../results/19-07-29-ExperimentalResultsOptimizationGfpMajo/results_data/SecondQuantificationFinalGFPs.xlsx", col_names = c("sequence", "intensity"))

results <- 
  experiment %>% 
  filter(intensity != "err") %>%
  mutate(
    seq_id = str_replace(sequence, "GFP", ""),
    seq_id = str_c(seq_id, "GFP"),
    intensity = as.numeric(intensity)
  )

write_csv(results, "results-data/gfps-experiments-intensities.csv")

# sequences ---------------------------------------------------------------
# now generate the prediction with the new algorthm (not overfitting)

predictor <- predict_stability(specie = 'human')
gfps <- read_csv("../../results/19-08-12-ProteinOptimizationV2//results_data/gfps_predictions_uncertanty.csv") %>% 
  select(seq_id, seqs) %>% 
  mutate(
    prediction = map_dbl(seqs, predictor),
    half_life = map_dbl(prediction, unscale_decay_to_mouse)
  )

write_csv(gfps, "results-data/gfp_sequences.csv")
