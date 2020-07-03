library(tidyverse)

write_rds(snakemake, "logs/03-make-plot.rds")
infile <- snakemake@input[[1]]
outplot <- snakemake@output[[1]]

sylamer <- read_csv("results-data/syll_aggregated.csv") %>% 
  mutate(
    time = str_extract(time, "\\d\\d?") %>% as.numeric()
  ) %>% 
  filter(time <= 9, time > 3)

elementos <- bind_rows(
  tibble(kmer=c("GGACTT", "TAGGAC"), motif = "m6a"),
  tibble(kmer=c("TCTATC", "CTATCT", "TATCTA"), motif = "m5c"),
  tibble(kmer=c("GCACTT"), motif = "miR-430")
) %>% 
  inner_join(sylamer) %>% 
  mutate(time = paste0(time, " hrs"))

set.seed(11)
sample_mers <- unique(sylamer$kmer) %>% 
  sample(500)

# make the plot supplemental plot -----------------------------------------

sylamer %>% 
  filter(
    kmer %in% sample_mers,
    !kmer %in% elementos$kmer
  ) %>% 
  mutate(time = paste0(time, " hrs")) %>% 
  ggplot(aes(x = rank, y = log10pval, group = kmer)) +
  geom_line(position = "jitter", color = "grey", size = 1/10) +
  geom_hline(yintercept = c(-1, 1) * log10(0.001), linetype=3, size=1/3) +
  geom_line(
    data = elementos,
    aes(x=rank, y=log10pval, group=kmer, color=motif, size=kmer),
    size=.4
  ) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 2000, 3500), labels = c("0k", "2k", "3.5k")) +
  scale_y_continuous(breaks = c(-2, 0, 5, 10)) +
  scale_color_manual(values = c("#7F00F9", "#E69F00", "#009E73")) +
  facet_grid(time~species, scales = "free_x") +
  theme_bw(base_family = 'Helvetica') +
  theme(panel.grid = element_blank())

ggsave(outplot, height = 6, width = 4)  
