library(tidyverse)
library(ggthemes)

theme_set(theme_tufte(base_family = "Helvetica") + theme(axis.line = element_line(size = .3)))

d <- read_csv("../200608-process-arielome/data/arielome_opt_mir.csv")

# make the time var early or late

d <- 
  d %>%
  mutate(
    time = map_chr(time, ~if_else(. %in% c(2, 1), "early", "late")),
    condition == "wt"
  )



# boostrap ----------------------------------------------------------------

boostrap_sample <- function(i) {
  sample_frac(d, prop = 1, replace = T) %>% 
    mutate(
      opt = ntile(optimality, n = 7)
    ) %>% 
    group_by(species, time, opt) %>% 
    summarise(mir_p = mean(mir430mer6)) %>% 
    ungroup() %>% 
    mutate(boot = i)
}

resultado <- 1:100 %>% 
  map_df(boostrap_sample)

resultado %>% 
  pivot_wider(values_from = mir_p, names_from = time) %>% 
  mutate(
    log2 = log2((late + .000001) / (early + .000001))
  ) %>% 
  ggplot(aes(x = opt, y = log2)) +
  geom_hline(yintercept = 0, size = .1, linetype = 2) +
  geom_boxplot(aes(group = opt, fill = opt), outlier.shape = NA, size = .3) +
  scale_x_continuous(breaks = 1:15) +
  facet_wrap(~species, scales = "free_y") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 3.5)

ggsave("figures/03-result-maon.pdf", height = 2, width = 5)
