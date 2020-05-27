library(boot)
library(tidyverse)

# get a confidence interval  ----------------------------------------------

dta <- read_csv("results_data/prediction_interavals_test_data.csv")

rsq <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample 
  fit <- lm(formula, data=d)
  return(summary(fit)$r.square)
} 
# bootstrapping with 1000 replications 
results <- boot(data=dta, statistic=rsq, 
                R=10000, formula=y_observed~median_prediction)

# view results
results 
plot(results)

# get 95% confidence interval 
boot.ci(results, type="bca")

# the result is 95% ( 0.3515,  0.3910


)
# do the estimation for the model with the 3’ UTR regulation --------------

dta2 <- read_csv("../19-07-18-PredictiveModelWithM6AandMicroRNAs/results_data/val_gbm.csv")

results2 <- boot(data=dta2, statistic=rsq, 
                R=10000, formula=observed~predicted)
boot.ci(results2, type="bca")

# the result is 
# 95%   ( 0.3601,  0.3998 ) 