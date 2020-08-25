library(tidyverse)
library(Biostrings)


# params ------------------------------------------------------------------
infile <- snakemake@input$seqs
specie_ <- snakemake@wildcards$species
seqtype <- snakemake@wildcards$seqtype
outfile <- snakemake@output[[1]]

# load the data -----------------------------------------------------------

seqs <- read_csv(infile) %>% 
  filter(!is.na(`3utr`), specie == specie_) %>% 
  select(gene_id, `3utr`, coding, specie)


# make a window of -+200bp around the stop codon --------------------------

seqs <- seqs %>% 
  mutate(
    stopcodon = str_c(
      str_sub(coding, -200), # the last 200 bp of the coding
      str_sub(`3utr`, start = 1, end = 50) # the first 50 bp
    )
  )


# retreive the seqs for given seqtype -------------------------------------

secuencias <- seqs[, seqtype] %>% 
  pull(1)
secuencias <- DNAStringSet(secuencias)
names(secuencias) <- seqs$gene_id

writeXStringSet(secuencias, outfile)

