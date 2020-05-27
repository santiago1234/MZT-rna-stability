library(brms)
library(tidybayes)

mir1 <- mutate(mir1, )


# plot the response distribution ------------------------------------------

mir1 %>% 
  ggplot(aes(x = sites,y=expression_level)) +
  geom_boxplot() +
   scale_y_log10() +
  facet_grid(.~sample_name)



# tidy the data for the model ---------------------------------------------

data <- 
  mir155 %>% 
  mutate(
    sites = `Conserved sites total`,
    time = str_extract(sample_name, "..hr"),
    condition = str_extract(sample_name, "(mock|mir[0-9]+)")
  ) %>% 
  select(gene_id, expression_level, sample_name, time, sites, condition) %>% 
  filter(expression_level > 0) %>% 
  mutate(logrna = log2(expression_level))


fml <- bf(logrna ~ 1 + (1 | time + sites + condition))

fit <- brm(fml, data = data, family = gaussian(), chains = 2, cores = 2)



fit %>%
  spread_draws(b_Intercept, r_sites[sites,Intercept]) %>% 
  median_qi() %>% 
  ggplot(aes(x=r_sites + b_Intercept, y=sites)) +
  geom_errorbarh(aes(xmin = r_sites.lower + b_Intercept.lower, xmax=r_sites.upper + b_Intercept.upper)) +
  geom_point()
