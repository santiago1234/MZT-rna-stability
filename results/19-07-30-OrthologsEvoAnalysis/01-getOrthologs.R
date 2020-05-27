library(tidyverse)
library(ggthemes)
library(gridExtra)

theme_set(theme_tufte(base_family = 'Helvetica'))

# get the orthologs -------------------------------------------------------

orto <- read_csv("results_data/fish_xenopus_ortologs.csv") %>% 
  select(`Gene stable ID`, `Xenopus gene stable ID`) %>% 
  distinct() %>% 
  rename(Gene_ID = `Gene stable ID`, xen_ort = `Xenopus gene stable ID`)

## load the maternal genes only

fish_maternal <- read_csv("../19-02-24-OverlapPathwaysFig3/results_data/regulatory_pathways_matrix.csv") %>% 
  filter(is_maternal) %>% 
  select(Gene_ID) %>% 
  distinct()

orto <- inner_join(fish_maternal, orto)

# keep just 1 to 1 mapping

to_keep <- 
  orto %>% 
  group_by(Gene_ID) %>% 
  summarise(n=n()) %>% 
  filter(n == 1) %>% 
  pull(Gene_ID)

orto <- filter(orto, Gene_ID %in% to_keep)

# load the log2 fold change data ------------------------------------------

xen_rna <- read_csv("../../data/19-06-12-XenopusData/time_course_xenopus.csv")

xen_rna <- 
  xen_rna %>% 
  rename(Gene_ID = ensembl_gene_id) %>% 
  filter(
    time %in% c(2, 9),
    sample_condition == "wt_polya"
  ) %>% 
  spread(key = time, value = expression_quantification) %>% 
  mutate(
    log2FC = log2(`9` / (`2` + 0.00001))
  ) %>% 
  select(Gene_ID, log2FC) %>% 
  rename(xen_ort = Gene_ID, xen_fc = log2FC)

# fish log2 FC ------------------------------------------------------------

fish_fc <- read_csv("../../data/19-02-05-FoldChangeData/data/log2FC_earlyVSlate_tidytimecourse.csv") %>% 
  filter(
    sample_condition == "wt_polya",
    time == 6,
    is_maternal
  ) %>% 
  select(Gene_ID, log2FC) %>% 
  rename(fish_fc = log2FC)


FC_data <- 
  left_join(orto, xen_rna) %>% 
  left_join(fish_fc) %>% 
  filter(xen_fc < 10)


FC_data %>% 
  ggplot(aes(x=fish_fc, y=xen_fc)) +
  geom_point(size=1/3) +
  geom_smooth() +
  ggpubr::stat_cor()


ggsave("cor_fc.pdf", height = 3, width = 3)


# get the coding for those genes ------------------------------------------

orf_fish <- read_csv("../../data/19-01-17-Get-ORFS-UTRS-codon-composition/sequence-data/fish_seq_data_cds_3utr.csv") %>% 
  rename(Gene_ID = ensembl_gene_id, f_orf = coding, f_utr = `3utr`)

orf_xen <- read_csv("../../data/19-06-12-XenopusData/xenopus_seq_data_cds_3utr.csv") %>% 
  rename(xen_ort = ensembl_gene_id, x_orf = coding, x_utr = `3utr`)

FC_data <- 
  left_join(FC_data, orf_fish) %>% 
  left_join(orf_xen)

# compute the optimality ration according to EMBO 2016 --------------------


opt_ratio <- function(orf) {
  # computes the optimality ratio according to bazzini EMBO 2016
  # r = log2(n optnal / n non optiaml)
  orf %>% 
    codonr::count_codons() %>% 
    right_join(codonr::optimality_code_embo2016) %>% 
    replace_na(list(n=0)) %>% 
    filter(!optimality %in% c("stop", "neutral")) %>% 
    group_by(optimality) %>%
    summarise(total = sum(n)) %>% 
    spread(key = optimality, value = total) %>% 
    mutate(opt_ratio = log2(`optimal` / `non-optimal`)) %>% 
    pull(opt_ratio)
}

FC_data <- 
  FC_data %>% 
  filter(!is.na(x_orf)) %>% 
  mutate(
    opt_fish = map_dbl(f_orf, opt_ratio),
    opt_xen = map_dbl(x_orf, opt_ratio)
  )

# plot optimality ratio ---------------------------------------------------

p1 <- FC_data %>% 
  ggplot(aes(x=opt_xen, y=xen_fc)) +
  geom_point() +
  ggpubr::stat_cor() +
  geom_rangeframe(color="black", alpha=1/4) +
  labs(x='optimality Xen', y = "log2 FC Xen")


p2 <- FC_data %>% 
  ggplot(aes(x=opt_fish, y=fish_fc)) +
  geom_point() +
  ggpubr::stat_cor() +
  geom_rangeframe(color="black", alpha=1/4) +
  labs(x='optimality fish', y = "log2 FC fish")


p3_opt <- FC_data %>% 
  ggplot(aes(opt_fish, opt_xen, color=fish_fc)) +
  geom_point() +
  scale_color_viridis_c() +
  ggpubr::stat_cor() +
  geom_rangeframe(color="black", alpha=1/4) +
  geom_smooth(method = "lm", color="black") +
  labs(x="fish Optimality", y = "xen Optimality")

p4_fc <- FC_data %>% 
  ggplot(aes(x=fish_fc, y=xen_fc)) +
  geom_point(size=1/3) +
  geom_rangeframe(color="black", alpha=1/4) +
  ggpubr::stat_cor()


pdf("result.pdf", width = 8, height = 6)
grid.arrange(p1, p2, p3_opt, p4_fc)
dev.off()
