library(tidyverse)
library(ggthemes)

theme_set(theme_tufte(base_family = "Helvetica"))

load_gfp <- function(dpath, replica_id) {
  readxl::read_excel(dpath, col_names = c("GFP", "Intensity")) %>% 
    filter(Intensity != "err") %>% 
    mutate(
      Intensity = as.numeric(Intensity),
      Replica = replica_id
    )
}


gfps <- 
  bind_rows(
  load_gfp("results_data/firstQuantificationGFPs.xlsx", "A"),
  load_gfp("results_data/SecondQuantificationFinalGFPs.xlsx", "B")
)


gfps %>% 
  ggplot(aes(x=GFP, y=Intensity)) +
  geom_point() +
  facet_grid(~Replica) +
  scale_y_log10()


# only B (A esta del asco) ------------------------------------------------

gfps <- filter(gfps, Replica == "B")

gfp_mean <- gfps %>% 
  group_by(GFP) %>% 
  summarise(m_i = mean(Intensity))

gfps %>% 
  filter(Replica == "B") %>% 
  ggplot(aes(x=GFP, y=Intensity, color=log10(Intensity))) +
  geom_jitter(width = 1/4, shape=16, alpha=4/5) +
  scale_color_viridis_c(option = "D", direction = -1) +
  scale_y_log10() +
  geom_rangeframe(color="black", alpha=3/4) +
  labs(
    y="Intensity",
    x = "synonymous GFP",
    title = "mRNA stability optimization"
  )

ggsave("gfp_res.pdf", height = 3, width = 4.5)


# plot Q results ----------------------------------------------------------

q_res <- readxl::read_excel("results_data/Q_silent_Reporters.xlsx") %>% 
  select(X__1, `mcherry/GFP`) %>% 
  rename(GFP = X__1, Intensity = `mcherry/GFP`) %>% 
  mutate(experiment = "silent Reporter\n7X")

gfps %>% 
  mutate(experiment = "GFPs\n1500X") %>% 
  bind_rows(q_res) %>% 
  ggplot(aes(x=reorder(GFP, Intensity, mean), y=Intensity, color=log10(Intensity))) +
  geom_jitter(width = 1/4, shape=16, alpha=4/5) +
  facet_wrap(~experiment, scales="free_x") +
  scale_y_log10() +
  scale_color_viridis_c(option = "C", direction = -1) +
  geom_rangeframe(color="black") +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(
    x = "Reporter Sequence",
    y = "log10 Intensity"
  )
ggsave("gfp_with_synRepsres.pdf", height = 3, width = 5.5)


# plot distance matrix GFPs -----------------------------------------------

d_gfp <- read_csv("results_data/DistanceData.csv") %>% 
  filter(X1 != "IDT_optimized_neutral_GFP") %>% 
  select(-IDT_optimized_neutral_GFP)

same_gfps <- d_gfp %>% 
  gather(key = GFP, value = d, -X1) %>% 
  filter(is.na(d) & X1 == GFP) %>% 
  replace_na(list(d=0))

d_gfp <- 
  d_gfp %>% 
  gather(key = GFP, value = d, -X1) %>% 
  filter(!is.na(d)) %>% 
  bind_rows(same_gfps)


d_gfp %>% 
  mutate(
    GFP = factor(GFP, levels = c('eGFP', "supremaGFP", "infimaGFP", "neutral_GFP")),
    X1 = factor(X1, levels = c("neutral_GFP", "infimaGFP", "supremaGFP", 'eGFP'))
  ) %>% 
  ggplot(aes(x=X1, y=GFP, fill=d, label=d)) +
  geom_tile() +
  geom_text(color="white") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Number of nucleotide differences"
  )

ggsave("n_mutations_gfps.pdf", height = 2.5, width = 3)

