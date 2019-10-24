library(tidyverse)


ranks <- list.files(path = "results_data", pattern = "*rank_by*", full.names = T)

names_rankes <- ranks %>% 
  basename() %>% 
  str_remove("rank_by_") %>% 
  str_remove(".txt")

names(ranks) <- names_rankes

# define sylamer vars

sylamer_bin <- "./../../bin/sylamer-12-342/sylamer"

run_sylamer <- function(rank_file, id_name, k = 6) {
  
  out_file_name <- file.path("results_data/", paste0("sylamerResult_", id_name, "_k", k, ".tsv"))
  
  utrs_path <- if_else(str_detect(rank_file, "fish"), "results_data/utrs_fish.fa", "results_data/utrs_xenopus.fa")
  
  sylamer_cmd <- paste(
    sylamer_bin,
    "-k ",
    k,
    " -m 4 -grow-times 20 -v 0 -fasta",
    utrs_path, 
    "-universe",
    rank_file,
    "-o", 
    out_file_name
  )
  
  #system(sylamer_cmd)
  sylamer_cmd
  
}


to_run_sylmer <- 
  tibble(ranks = ranks, id_name = names(ranks)) %>% 
  crossing(k = 6:7)
  
to_run_sylmer <- 
  to_run_sylmer %>% 
  mutate(
    run_syl = purrr::pmap_chr(list(ranks, id_name, k), run_sylamer)
  )

# execute the run with system
map(to_run_sylmer$run_syl, .f = system)
