library(tidyverse)
library(Biostrings)

# load the data -----------------------------------------------------------

datos <- read_csv("../191010-PredictStabilityInMZT/results-data/mzt_predictionsResidualsLog2Fc.csv") %>% 
  filter(!is.na(`3utr`))

generate_data_for_specie <- function(specie_) {
  # this function ranks the gene ids by 3 variables:
  # log2FC, residual and predicted stability
  # and saves the data including the 3' UTRs to run sylamer
  
  datos_specie <- 
    datos %>% 
    dplyr::filter(specie == specie_)
  
  utrs <- DNAStringSet(datos_specie$`3utr`)
  names(utrs) <- datos_specie$gene_id
  writeXStringSet(utrs, str_c("results_data/utrs_", specie_, ".fa"))
  ## rankings
  datos_specie %>% 
    arrange(log2FC) %>% 
    pull(gene_id) %>% 
    write_lines(
      file.path("results_data", str_c(specie_, "_rank_by_observedlog2FC.txt"))
    )
  
  datos_specie %>% 
    arrange(predicted_stability) %>% 
    pull(gene_id) %>% 
    write_lines(
      file.path("results_data", str_c(specie_, "_rank_by_predictedStability.txt"))
    )
  
  datos_specie %>% 
    arrange(resid) %>% 
    pull(gene_id) %>% 
    write_lines(
      file.path("results_data", str_c(specie_, "_rank_by_residual.txt"))
    )
  
}

map(c("fish", "xenopus"), generate_data_for_specie)