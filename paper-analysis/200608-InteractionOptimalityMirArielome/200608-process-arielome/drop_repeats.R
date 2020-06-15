library(tidyverse)

# In this script i take advantage that if a sequence shares the same gene id, then
# the sequences (cds & utr) are the same.


# input parameter ---------------------------------------------------------

infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]
minrepetition <- snakemake@params$min_repeat

# load data ---------------------------------------------------------------

d <- read_csv(infile)

# run code ----------------------------------------------------------------


d_counts <- d %>% 
  count(seq_id)

to_keep <- d_counts %>% 
  filter(n <= minrepetition) %>% 
  pull(seq_id)


# save sequences ----------------------------------------------------------

d %>% 
  filter(seq_id %in% to_keep) %>% 
  write_csv(outfile)
