library(tidyverse)
library(ggthemes)
library(gridExtra)

theme_set(theme_tufte(base_family = 'Helvetica'))


sylamer_results <- list.files("../191012-Fig02SylamerAnalysis/results_data", pattern = "sylamerResult*", full.names = T)

load_syla_file <- function(syl_file) {
  id_name <- basename(syl_file) %>% 
    str_remove("sylamerResult_") %>% 
    str_remove(".tsv") %>% 
    str_remove("_k\\d")
  
  read_tsv(syl_file) %>% 
    mutate(utrs_sorted_from = id_name) %>% 
    gather(key = rank, value = log10pavl, -upper, -utrs_sorted_from) %>% 
    mutate(rank = as.integer(rank)) %>% 
    dplyr::rename(kmer = upper) %>% 
    mutate(
      pvalue = 10**-(abs(log10pavl)),
      rank = as.integer(rank),
      k = str_length(kmer)
    ) %>% 
    separate(utrs_sorted_from, into = c("specie", "ranked_by"), sep = "_")
  
}


sylamer <- map_df(sylamer_results, load_syla_file) %>% 
  filter(k == 6)

syla <- sylamer
# normalize rank position

syla <- 
  syla %>% 
  group_by(specie) %>% 
  mutate(rank = rank / max(rank)) %>% 
  ungroup()


# get the peak with more signal

syla <- 
  syla %>% 
  group_by(specie, kmer, ranked_by) %>% 
  filter(abs(log10pavl) == max(abs(log10pavl))) %>% 
  group_map(~head(., 1), keep = T) %>% # get the first row of each group, in some case
  # i get more than 1 row
  reduce(bind_rows)

# get the mean position of the peek in both specie

syla <- syla %>% 
  group_by(kmer, ranked_by) %>% 
  mutate(relative_peak_position = mean(rank))

# tidy the data for visualization

syla <- 
  syla %>% 
  select(-k, -rank, -pvalue) %>% 
  spread(key = specie, value = log10pavl)

# visualization -----------------------------------------------------------


syla %>% 
  filter(ranked_by == "residual") %>% 
  ggplot(aes(x=xenopus, y=fish, color=relative_peak_position)) +
  geom_point(size=1/2) +
  scale_color_viridis_c(option = "B") +
  coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme(legend.position = "bottom")
  
ggsave("peak-with-high-signal.pdf", height = 3, width = 4)  
# visualize m6a like and m5c like, mir430 like ----------------------------


sylamer <- 
  sylamer %>% 
  filter(k == 6)

mr430_like <- sylamer %>% 
  filter(str_detect(kmer, "GCACT"))
  

m6a_like <- sylamer %>% 
  filter(str_detect(kmer, "GGACT"))

m5c_like <- sylamer %>% 
  filter(str_detect(kmer, "CATCT"))


do_plot <- function(df_like, title) {
  df_like %>% 
    ggplot(aes(x=rank, y=log10pavl, group=kmer, color=kmer)) +
    geom_line(data  = sample_n(sylamer, 1e5), aes(), color="grey", size=1/10) +
    geom_line(size=1/2) +
    geom_hline(yintercept = c(-1, 1) * log10(.01), linetype=3) +
    scale_color_colorblind() +
    facet_wrap(specie~ranked_by,scales = "free_x") +
    coord_cartesian(ylim = c(-8, 8)) +
    labs(
      title = title,
      color = title
    )
}

p1 <- do_plot(mr430_like, "miR-430 (GCACT)")
p2 <- do_plot(m6a_like, "m6A (GGACT)")
p3 <- do_plot(m5c_like, "m5C (CATCT)")

pdf("destabalizing.pdf", height = 12, width = 7)
grid.arrange(p1,p2,p3, ncol=1)
dev.off()
