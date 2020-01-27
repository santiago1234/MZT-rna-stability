plot(model)

dt2 <- datum %>% 
  mutate(
    utrlen = str_length(`3utr`)
  ) %>% 
  filter(specie == "fish")

model <- brm(data = dt2, family = binomial,
          miR430 | trials(utrlen) ~ 1 + PLS1 + PLS2 + log2FC,
          prior(normal(0, 10), class = Intercept),
          seed = 10, cores = 2, chains = 2)
plot(model)

marginal_effects(model)

dt2 %>% 
  ggplot(aes(PLS1, PLS2)) +
  geom_point()

topredict <- tibble(
  PLS1 = seq(-6, 6, length.out = 50),
  PLS2 = seq(-6, 6, length.out = 50)
  
) %>% 
  crossing(
    tibble(log2FC = seq(-6, 6, length.out = 50))
  ) %>% 
  mutate(utrlen = median(dt2$utrlen))


topredict <- 
  predict(model, newdata = topredict) %>% 
  as_tibble() %>% 
  bind_cols(topredict)

smarca2 <- dt2 %>% 
  filter(gene_id == "ENSDARG00000008904") %>% 
  mutate(Estimate=10) %>% 
  mutate(name = "smarca2")

topredict %>% 
  ggplot(aes(x=log2FC, y=PLS1, color=Estimate)) +
  geom_point(shape=15, size=4) +
  geom_text(data = smarca2, aes(x = log2FC, y = PLS1, label=name), color="white") +
  geom_point(data = smarca2, aes(x = log2FC, y = PLS1), color="white") +
  scale_color_viridis_c(option = "A")


paramspreds <- 
  topredict %>% 
  select(PLS1:utrlen) %>% 
  select(-log2FC) %>% 
  distinct() %>% 
  crossing(tibble(log2FC = c(-5, 0, 5)))

res2 <- fitted(model, newdata = paramspreds) %>% 
  as_tibble() %>% 
  bind_cols(paramspreds)

res2 %>% 
  ggplot(aes(x=PLS1, y=Estimate)) +
  geom_line() +
  geom_line(aes(y=Q2.5), linetype=2) +
  geom_line(aes(y=Q97.5), linetype=2) +
  facet_grid(log2FC~.)


