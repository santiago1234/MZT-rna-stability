library(tidyverse)

# script parameters
meta <- snakemake@input$metadata
infiles <- snakemake@input$ari
outfile <- snakemake@output[[1]]
  
  

meta <- read_csv(meta)


loadf <- function(infile) {
  samname <- basename(infile) %>% 
    str_replace(pattern = "\\.1_optimality_mir\\.csv", replacement = "")
  
  read_csv(infile) %>% 
    mutate(sample_name = samname)
}

alldat <- infiles %>% 
  map_df(loadf) %>% 
  inner_join(meta) %>% 
  select(sample_name, time, condition, everything())

write_csv(alldat, outfile)