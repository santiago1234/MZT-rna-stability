library(tidyverse)
library(ggthemes)

theme_set(theme_tufte())

sylreporters <- read_csv("silent_predictions.csv") %>% 
  filter(!id_seq %in% c("infima", "suprema"))


sylreporters %>% 
  mutate(p = str_extract(id_seq, pattern = "\\d\\d\\d?") %>% as.numeric()) %>% 
  ggplot(aes(x=prediction, y=`mcherry/GFP`)) +
  #geom_smooth(method = "lm", color="grey10", size=1/4) +
  geom_point(aes(color=p), shape = 16, alpha = .99) +
  scale_color_gradient2(low = "blue", high = "red", midpoint = 60) +
  ggrepel::geom_text_repel(aes(label=id_seq)) +
  scale_y_log10() +
  geom_rangeframe() +
  ggpubr::stat_cor(method = "spearman") +
  labs(
    y = "log10 mCherry/GFP",
    x = "predicted stability"
  ) +
  theme(legend.position = "none")

ggsave("sylentscatter.pdf", height = 2, width = 3)

