library(tidyverse)
library(ggthemes)

theme_set(theme_tufte(base_family = "Helvetica"))

# gather the data ---------------------------------------------------------
# xenopus data
expression_xen <- read_csv("../../data/19-06-12-XenopusData/time_course_xenopus.csv")
pathways_xen <- read_csv("../../data/19-06-12-XenopusData/regulatory_pathways_xenopus.csv") %>% 
  select(-PLS1, -PLS2)

dataset_xen <- expression_xen %>% 
  rename(Gene_ID = ensembl_gene_id) %>% 
  # estimate the log2 fold change
  filter(
    time %in% c(2, 9)
  ) %>% 
  spread(key = time, value = expression_quantification) %>% 
  mutate(
    log2FC = log2(`9` / (`2` + 0.00001))
  ) %>% 
  left_join(pathways_xen) %>% 
  filter(
    !is.na(MiR430),
    in_top_4k_expressed_at_0hrs,
    sample_condition == "wt_ribo"
  ) %>% 
  select(Gene_ID, log2FC, MiR430) %>% 
  mutate(
    target = MiR430 > 0,
    pathway = "miR430"
  ) %>% 
  select(-MiR430)


# fish data

fish_fc <- read_csv("../../data/19-02-05-FoldChangeData/data/log2FC_earlyVSlate_tidytimecourse.csv") %>% 
  filter(
    sample_condition == "wt_ribo",
    is_maternal,
    time == 6
  ) %>% 
  select(-is_maternal, -sample_condition, -time)

# fish pathways

dataset_fish <- 
  read_csv("../19-02-24-OverlapPathwaysFig3/results_data/regulatory_pathways_matrix.csv") %>% 
  filter(is_maternal) %>% 
  select(Gene_ID, m6A, MiR430) %>% 
  mutate(miR430 = MiR430 > 0) %>% 
  select(-MiR430) %>% 
  gather(key = pathway, value = target, -Gene_ID) %>% 
  inner_join(fish_fc)


datasets <- bind_rows(dataset_fish, dataset_xen) %>% 
  rename(gene_id = Gene_ID)

# load optimality (PLS) ---------------------------------------------------

optimality <- read_csv("results_data/pls_species.csv")
data <- inner_join(datasets, optimality)

# run the analysis --------------------------------------------------------

# make optimality quantiles
data <- 
  data %>% 
  mutate(pathway = factor(pathway, levels = c("miR430", "m6A"))) %>% 
  group_by(specie, pathway) %>% 
  mutate(optimality = ntile(-PLS1, 4)) %>% 
  ungroup()

data %>% 
  filter(optimality %in% c(1, 4)) %>% # plot the extreme quantiles
  ggplot(aes(x=target, y=log2FC, color=as.character(optimality))) +
  geom_tufteboxplot() +
  geom_rangeframe(sides="l", color="black", alpha=2/3) +
  scale_x_discrete(labels = c("targets", "no targets")) +
  scale_color_manual(values = c("red", "blue")) +
  facet_grid(pathway~specie) +
  labs(
    title = "combinatorial code MZT",
    y = "log2 fold change\n late stage (fish = 6 hrs, xenopus = 9 hrs) / 2 hrs"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
    )

ggsave("figures/mzt_overlap.pdf", width = 4, height = 4)


data %>% 
  group_by(specie, pathway, target) %>% 
  nest() %>% 
  mutate(
    fit = map(data, ~lm(log2FC ~ PLS1, data = .)),
    tf = map(fit, broom::tidy)
  ) %>% 
  unnest(tf) %>% 
  filter(term == "PLS1")
