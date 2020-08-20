library(tidyverse)
library(optimalcodonR)
library(ggthemes)

theme_set(theme_tufte(base_family = "Helvetica") + theme(axis.line = element_line(size = .3)))

d <- read_csv("../200608-process-arielome/data/arielome_opt_mir.csv")

# make the time var early or late

d <- 
  d %>%
  mutate(
    time = map_chr(time, ~if_else(. %in% c(2, 1), "early", "late")),
    condition == "wt" # some samples are morfolino
  )


# plot optimality ---------------------------------------------------------
# endogenous genes

predictor <- predict_stability(specie = "fish")

endogenes <- 
  training %>% 
  filter(specie == "fish") %>% 
  select(coding) %>% 
  unique()
endogenes$optimality <- predictor(endogenes$coding)


# make plot comparing distribution to endogenous genes --------------------

tmp <- d %>% 
  select(optimality) %>% 
  mutate(source = "arielome")

tmp <- endogenes %>% 
  select(optimality) %>% 
  mutate(source = "endogenous") %>% 
  bind_rows(tmp)

tmp %>% 
  ggplot(aes(x = optimality)) +
  geom_density(aes(linetype = source, color = source), size = .5) +
  scale_color_manual(values = c("steelblue", "grey30"))

ggsave("figures/01-ComparingOptimalityToEndogeounes.pdf", height = 2, width = 4)


# now look at the time, do we see mirR430 depletion? ----------------------
# draw median line

mopt <- d %>% 
  group_by(species, time) %>% 
  summarise(
    mopt = median(optimality)
  )

d %>% 
  ggplot(aes(x = optimality, color = time)) +
  geom_density() +
  geom_vline(
    data = mopt,
    aes(xintercept = mopt, color = time),
    linetype = 2, size = .3) +
  scale_color_colorblind() +
  facet_grid(. ~ species)

ggsave("figures/01-DoWeSeeOptimality.pdf", height = 2, width = 5)


# compute p-values --------------------------------------------------------


d %>% 
  group_by(species) %>% 
  nest() %>% 
  mutate(
    fit = map(data, ~lm(optimality ~ time, data = .)),
    tf = map(fit, broom::tidy)
  ) %>% 
  unnest(tf) %>% 
  select(-c(data, fit)) %>% 
  write_csv("results-data/stats-optimality.csv")


# compute the p-value for mirDepletion ------------------------------------

d %>% 
  group_by(species) %>% 
  nest() %>% 
  mutate(
    fit = map(data, ~glm(mir430mer6 ~ time , data = ., family=binomial(link='logit'))),
    tf = map(fit, broom::tidy)
  ) %>%
  unnest(tf) %>% 
  select(-c(data, fit)) %>% 
  write_csv("results-data/stats-mir.csv")

## draw plot

boostra_mir_dep <- function(index)  {
  d %>% 
    sample_frac(size = 1, replace = T) %>% 
    filter(condition == "wt") %>% 
    group_by(species, time) %>% 
    count(mir430mer6) %>% 
    mutate(
      p = n / sum(n) # this gives the probability of mir in each sample
    ) %>% 
    ungroup() %>% 
    filter(mir430mer6) %>% 
    select(-mir430mer6, -n) %>% 
    pivot_wider(names_from = time, values_from = p) %>% 
    mutate(
      log2fc = log2(late / early),
      i = index
    )
}

mir_dep_res <- 
  1:100 %>% 
  map_df(boostra_mir_dep)


mir_dep_res %>% 
  ggplot(aes(x = species, y = log2fc)) +
  geom_hline(yintercept = 0) +
  geom_boxplot(fill = "grey")
ggsave("figures/02-MirDepletion.pdf", height = 2, width = 1.5)
