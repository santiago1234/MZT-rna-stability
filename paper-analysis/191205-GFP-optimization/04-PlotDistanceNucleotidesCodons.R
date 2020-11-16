library(tidyverse)
library(Biostrings)
library(optimalcodonR) # tienes que instalar esto: https://github.com/santiago1234/optimalcodonR

# load the sequences and save them as fasta -------------------------------
# I manually created the file adding/removing codon to make the neutralGFP match
seqs <- readBStringSet("seqs-trimmed-for-counting-differences.csv") %>% 
  DNAStringSet()


# create table with function ----------------------------------------------
# la idea es tener una tabla con dos columnas: (nombre de secuencia) y (secuencia)
x_1 <- tibble(
  id_1 = names(seqs),
  seq_1 = as.character(seqs)
)

x_2 <- tibble(
  id_2 = names(seqs),
  seq_2 = as.character(seqs)
)

# aqui tengo una tabla con las combinaciones pareadas
sequencias <- crossing(x_1, x_2)

# uso dos funciones del paquete optimalcodonR
# checa este ejemplo
s1 <- sequencias$seq_1[1]
s2 <- sequencias$seq_2[2]

# cuantos nucleotidos son diferentes?
nucleotide_distance(s1, s2)

# cuantos codones son diferentes?
codon_distance(s1, s2)

## compute the distance
## unas funciones para correrlo en la tabla

distancias <- 
  sequencias %>% 
  mutate(
    dist_nucs = map2_dbl(seq_1, seq_2, function(x, y) nucleotide_distance(x, y)),
    dist_codons = map2_dbl(seq_1, seq_2, function(x, y) codon_distance(x, y))
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
  theme_minimal() +
  theme(legend.position = "none")

ggsave("distance-matrix.pdf", height = 3, width = 4)


