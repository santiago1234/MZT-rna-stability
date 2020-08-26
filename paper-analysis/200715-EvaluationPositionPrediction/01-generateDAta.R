library(tidyverse)
library(optimalcodonR)

split_seq_by_position <- function(aseq) {
  # this function split a sequence in three equal size regions
  # 5' region, medium regin, and 3' region
  # returns a tibble with the regions
  stopifnot(nchar(aseq) %% 3 == 0)
  

  ncodons <- nchar(aseq) / 3
  region_size <- ncodons %/% 3
  
  first <- str_sub(aseq, start = 1, end = region_size * 3)
  second <- str_sub(aseq, start = region_size * 3 + 1, end = region_size * 3 * 2)
  third <- str_sub(aseq, start = region_size * 3 * 2 + 1, end = -1)
  
  tibble(
    first = first,
    second = second,
    third = third
  )
  
}

# get the data for fish ---------------------------------------------------

testing_f <- testing %>% 
  filter(datatype == "aamanitin polya") %>% 
  select(gene_id, coding, decay_rate) %>% 
  mutate(
    allocation = "testing"
  )

training_f <- training %>% 
  filter(datatype == "aamanitin polya") %>% 
  select(gene_id, coding, decay_rate) %>% 
  mutate(
    allocation = "training"
  )

d <- bind_rows(testing_f, training_f)

# add positional info -----------------------------------------------------

d <- d %>% 
  mutate(
    position = map(coding, split_seq_by_position)
  ) %>% 
  unnest(position)

write_csv(d, "data/position_seqs.csv")
