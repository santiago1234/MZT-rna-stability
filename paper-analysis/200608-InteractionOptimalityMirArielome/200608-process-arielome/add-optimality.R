library(tidyverse)
library(optimalcodonR)
library(glmnet)

# get the optimaliy by using the predictive model
# we use the predicted stability as an indicator of
# codon optimality

infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]

arielome <- read_csv(infile)
fish_predictor <- predict_stability('fish') # fom optimcalcodonR


# first predict the stability of the sequences ----------------------------
# get unique sequences to make the computation faster

arielomeseqs <- arielome %>% 
  select(cds) %>% 
  unique()

arielomeseqs$optimality <- fish_predictor(arielomeseqs$cds)

arielome <- inner_join(arielome, arielomeseqs)


# add an indicator variable for the miR-430 -------------------------------

mir430kmer <- "GCACTT"
utr <- select(arielome, utr3) %>% 
  unique()

utr$mir430mer6 <- str_detect(utr$utr3, mir430kmer)

arielome <- inner_join(arielome, utr) %>% 
  select(-cds, -utr3)

write_csv(arielome, outfile)

