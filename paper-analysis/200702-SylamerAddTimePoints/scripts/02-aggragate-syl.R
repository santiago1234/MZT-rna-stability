library(tidyverse)

# save log file
write_rds(snakemake, "logs/02-aggregate-log.rds")
fish_samples <- snakemake@input$fish
xen_samples <- snakemake@input$xenopus
all_samples <- c(fish_samples, xen_samples)
outfile <- snakemake@output[[1]]

# function for processing data --------------------------------------------

get_sample_info <- function(sample_name) {
  sample_name <- sample_name %>% 
    basename()
  timepoint <- str_extract(sample_name, pattern = "\\d\\d?hrs")
  animal <- str_extract(sample_name, pattern = "(fish|xenopus)")
  c(time = timepoint, specie = animal)
}

load_syl <- function(syl_file) {
  syld <- read_tsv(syl_file)
  meta <- get_sample_info(syl_file)
  
  syld %>% 
    mutate(`0` = 0) %>% # Cei suggested this
    pivot_longer(-upper, names_to = "rank", values_to = "log10pval") %>% 
    mutate(rank = as.integer(rank)) %>% 
    dplyr::rename(kmer = upper) %>% 
    mutate(
      species = meta["specie"],
      time = meta["time"]
    )
}

# runall ------------------------------------------------------------------

all_samples %>% 
  map_df(load_syl) %>% 
  write_csv(outfile)
