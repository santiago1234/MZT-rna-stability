library(tidyverse)


# input -------------------------------------------------------------------
# save the snakemake object into a log file
write_rds(snakemake, "logs/01-rank-genes.rds")

# input parameters
predictions_by_time <- snakemake@input$timepoints %>% 
  read_csv()
time_point <- snakemake@wildcards$time %>% as.double() # we need a numeric type
animal <- snakemake@wildcards$specie # the species
outfile <- snakemake@output[[1]]

# subset the species and time point ---------------------------------------

pred_time <- predictions_by_time %>% 
  filter(
    time == time_point,
    specie == animal
  )

# if the table is empty invalid params
if (nrow(pred_time) == 0) {
  stop("Invalid time point or species given")
}

# now rank by the residuals -----------------------------------------------
pred_time %>% 
  arrange(residual) %>% 
  pull(gene_id) %>% 
  write_lines(outfile)

