library(tidyverse)
library(ggthemes)
library(ggforce)
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



data <- data %>% 
  mutate(optimality = str_c("q", optimality))

# mir-430 sites -----------------------------------------------------------

mir430_data <- data %>% 
  filter(pathway == "miR430")

### get summary stats for number of points and median values

summary_mir430 <- mir430_data %>% 
  group_by(specie, target, optimality) %>% 
  summarise(
    mediana = median(log2FC),
    n = n()
  )

mir430_data %>% 
  ggplot(aes(x=optimality, y=log2FC, color=optimality)) +
  geom_hline(yintercept = 0, size = .1, linetype = 2) +
  geom_sina(shape=16, alpha=.99, size=1/4) +
  geom_rangeframe(color="black", size=1/5, sides = "l") +
  geom_errorbar(
    data = summary_mir430,
    aes(ymin=mediana, ymax=mediana, y=mediana, x=optimality),
    color="black",
    size=1/5
  ) +
  geom_text(
    data = summary_mir430,
    aes(x = optimality, y = 4.5, label = paste0("n=", n)),
    color = "grey",
    size=2
  ) +
  scale_y_continuous(breaks = c(-4, 0, 4)) +
  scale_color_manual(values = c("#ca0020", "#f4a582", "#92c5de", "#0571b0")) +
  facet_grid(specie ~ target) +
  coord_cartesian(ylim = c(-5, 5)) +
  labs(
    title = "combinatorial code optimality and miR-430, MZT",
    y = "log2 fold change\n late stage (fish = 6 hrs, xenopus = 9 hrs) / 2 hrs",
    x = "codon optimality level"
  ) +
  theme(legend.position = "none", axis.ticks = element_line(size=1/5))

ggsave("figures/mzt_mir430_overlap.pdf", height = 3.5, width = 4)


# overlap m6a regulation --------------------------------------------------
m6a_data <- data %>% 
  filter(pathway == "m6A")

### get summary stats for number of points and median values

summary_m6a <- m6a_data %>% 
  group_by(specie, target, optimality) %>% 
  summarise(
    mediana = median(log2FC),
    n = n()
  )

m6a_data %>% 
  ggplot(aes(x=optimality, y=log2FC, color=optimality)) +
  geom_hline(yintercept = 0, size = .1, linetype = 2) +
  geom_sina(shape=16, alpha=.99, size=1/4) +
  geom_rangeframe(color="black", size=1/5, sides = "l") +
  geom_errorbar(
    data = summary_m6a,
    aes(ymin=mediana, ymax=mediana, y=mediana, x=optimality),
    color="black",
    size=1/5
  ) +
  geom_text(
    data = summary_m6a,
    aes(x = optimality, y = 4.5, label = paste0("n=", n)),
    color = "grey",
    size=2
  ) +
  scale_y_continuous(breaks = c(-4, 0, 4)) +
  scale_color_manual(values = c("#ca0020", "#f4a582", "#92c5de", "#0571b0")) +
  facet_grid(. ~ target) +
  coord_cartesian(ylim = c(-5, 5)) +
  labs(
    title = "combinatorial code optimality and m6A, MZT",
    y = "log2 fold change (6hrs / 2hrs)",
    x = "codon optimality level"
  ) +
  theme(legend.position = "none", axis.ticks = element_line(size=1/5))

ggsave("figures/mzt_m6a_overlap.pdf", height = 2.5, width = 4)



# linear model p-value ----------------------------------------------------

data %>% 
  group_by(specie, pathway, target) %>% 
  nest() %>% 
  mutate(
    fit = map(data, ~lm(log2FC ~ PLS1, data = .)),
    tf = map(fit, broom::tidy)
  ) %>% 
  unnest(tf) %>% 
  filter(term == "PLS1")
