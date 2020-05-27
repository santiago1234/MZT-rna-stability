library(tidyverse)
library(caret)


get_resample_performance_mdl <- function(path_to_mdl) {
  mdl_name <- path_to_mdl %>% 
    basename() %>% 
    str_replace(".rds", "")
  
  message("loading ", mdl_name, " ...")
  mdl <- read_rds(path_to_mdl)
  
  mdl$resample %>% 
    as_tibble() %>% 
    mutate(model = mdl_name)
  
}


trained_mdls <- list.files("results/trained-models", full.names = T)

resamples_res <- trained_mdls %>% 
  map_df(get_resample_performance_mdl)



resamples_res %>% 
  ggplot(aes(x=Rsquared, y=reorder(model, Rsquared, median), group=Resample)) +
  geom_line() +
  geom_point()


