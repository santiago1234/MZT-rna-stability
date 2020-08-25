library(tidyverse)

infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]

sylamer <- read_csv("results-data/syll_aggregated.csv") %>% 
  mutate(
    time = str_extract(time, "\\d\\d?") %>% as.numeric()
  )

m6a_motif <- "[GA]GAC[CT]"

m6a_like <- sylamer %>% 
  filter(str_detect(kmer, pattern = m6a_motif) | kmer == "GCACTT") %>% 
  mutate(
    seqtype = factor(seqtype, levels = c("coding", "stopcodon", "3utr"))
  )

m6a_like %>% 
  filter(species == "fish") %>% 
  ggplot(aes(x = rank, y = log10pval, group = kmer, color = kmer == "GCACTT")) +
  geom_line() +
  scale_color_manual(values = c("grey", "steelblue")) +
  labs(
    subtitle = "k-mers with m6a [GA]GAC[CT] motif\n miR-430 kmer shown in blue for comparison"
  ) +
  facet_grid(.~seqtype, scales = "free_x", space = "free_x") +
  theme(legend.position = 'none')

ggsave(outfile, height = 2, width = 6)


