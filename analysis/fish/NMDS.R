library(tidyverse)
library(here)
library(ggrepel)
library(vegan)

nmds <- metaMDS(fish.list$abund, distance="bray", k=2, trymax=1000, plot = TRUE)

points <- data.frame(nmds$points) %>% 
  mutate(site_month = rownames(fish.list$abund)) %>% 
  separate_wider_delim(site_month, delim = "_", names = c("site", "month"), cols_remove = FALSE) %>% 
  mutate(site_type = ifelse(site %in% c("MA", "WA", "HO", "TL", "LL", "PR"), "jubilee", "core"))

hulls <- points %>%
  group_by(site) %>% 
  slice(chull(MDS1,MDS2))

#all
ggplot() +
  geom_point(data = points, aes(x = MDS1, y = MDS2, color = site, shape = site_type)) + 
  geom_text_repel(data = points %>% filter(!site %in% c("MA", "WA", "HO", "TL", "LL", "PR")), aes(x = MDS1, y = MDS2, color = site, label = month)) +
  geom_text_repel(data = points %>% filter(site %in% c("MA", "WA", "HO", "TL", "LL", "PR")), aes(x = MDS1, y = MDS2, color = site, label = site)) +
  geom_polygon(data = hulls, aes(x = MDS1, y = MDS2, fill = site), alpha = 0.2) +
  theme_minimal() 

ggsave("docs/figures/fish_nmds_sites.png")

#core
ggplot() +
  geom_point(data = points %>% filter(!site %in% c("MA", "WA", "HO", "TL", "LL", "PR")), aes(x = MDS1, y = MDS2, color = site)) + 
  geom_text_repel(data = points %>% filter(!site %in% c("MA", "WA", "HO", "TL", "LL", "PR")), aes(x = MDS1, y = MDS2, color = site, label = month)) +
  geom_polygon(data = hulls %>% filter(!site %in% c("MA", "WA", "HO", "TL", "LL", "PR")), aes(x = MDS1, y = MDS2, fill = site), alpha = 0.2) +
  theme_minimal() 

# ggsave("docs/figures/nmds_core_sites.png")

