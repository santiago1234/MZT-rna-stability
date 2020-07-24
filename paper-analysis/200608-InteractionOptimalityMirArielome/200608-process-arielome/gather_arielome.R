library(tidyverse)


# script parameters
meta <- snakemake@input$metadata
meta_xen <- snakemake@input$metadata_xen

# combine fish and xenopus samples
infiles <- c(
  snakemake@input$ari,
  snakemake@input$ari_xen
)
outfile <- snakemake@output[[1]]
  
  

meta <- read_csv(meta) %>% mutate(species = "fish")
meta_xen <- read_csv(meta_xen) %>% mutate(species = "xenopus")

meta <- bind_rows(meta, meta_xen)

loadf <- function(infile) {
  samname <- basename(infile) %>% 
    str_replace(pattern = "\\.?1?_optimality_mir\\.csv", replacement = "")
  
  read_csv(infile) %>% 
    mutate(sample_name = samname)
}

alldat <- infiles %>% 
  map_df(loadf) %>% 
  inner_join(meta) %>% 
  select(sample_name, time, condition, everything())

write_csv(alldat, outfile)