library(tidyverse)
library(Biostrings)


# load the sequences and save them as fasta -------------------------------
# I manually created the file adding/removing codon to make the neutralGFP match
seqs <- readBStringSet("seqs-trimmed-for-counting-differences.csv") %>% 
  DNAStringSet()

seqs
s1 <- as.character(seqs)[3]
s2 <- as.character(seqs)[2]


get_nucs <- function(sequence) {
  seq(from=1, to=str_length(sequence), by = 1) %>% 
    map_chr(~str_sub(string = sequence, start = ., end = .))
}

get_codons <- function(sequence) {
  seq(from=1, to=str_length(sequence)-3+1, by = 3) %>% 
    map_chr(~str_sub(string = sequence, start = ., end = . + 2))
}


seq_distance <- function(s1, s2, iter_function) {
  sum(iter_function(s1) != iter_function(s2))
  
}

codon_distance
seq_distance(s1, s2, get_codons)



# create table with function ----------------------------------------------

x_1 <- tibble(
  id_1 = names(seqs),
  seq_1 = as.character(seqs)
)

x_2 <- tibble(
  id_2 = names(seqs),
  seq_2 = as.character(seqs)
)

sequencias <- crossing(x_1, x_2)

## compute the distance

distancias <- 
  sequencias %>% 
  mutate(
    dist_nucs = map2_dbl(seq_1, seq_2, function(x, y) seq_distance(x, y, get_nucs)),
    dist_codons = map2_dbl(seq_1, seq_2, function(x, y) seq_distance(x, y, get_codons))
  )

# grafico -----------------------------------------------------------------

distancias %>% 
  mutate(
    id_1 = factor(id_1, levels = c("supremaGFP", "neutralGFP", "infimaGFP", "eGFP")),
    labeler = paste0(dist_codons, "/",dist_nucs)
  ) %>% 
  ggplot(aes(x=id_1, y=id_2, fill = dist_codons)) +
  scale_x_discrete(expand = c(0, 0)) +
  geom_tile() +
  geom_text(aes(label=labeler), color="grey90") +
  scale_y_discrete(expand = c(0, 0)) +
  ggthemes::theme_tufte() +
  theme(legend.position = "none")

ggsave("distance-matrix.pdf", height = 3, width = 4)


