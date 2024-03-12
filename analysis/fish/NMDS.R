library(tidyverse)
library(here)
library(ggrepel)
library(vegan)

load(here("data", "fish.list.Rdata")) #object created in 02_tidy_data
fish.adj.abund <- fish.list[[2]] 
#need to add a dummy species because you can't calculate BC distances on a zero (SHR_04)
fish.adj.abund <- cbind(fish.adj.abund, dummy = rep(1, times = nrow(fish.adj.abund)))
fish.adj.abund <- sqrt(fish.adj.abund) #data transformation (sqrt) to avoid weighting by very abundant taxa

K <- 3 #number of dimensions, more than 3 generally gets unweildy
nmds <- metaMDS(fish.adj.abund, distance="bray", k= K, trymax=1000, plot = FALSE)
S <- round(nmds$stress, 2) #stress decreases with increasing K, >0.2 is generally a poor fit
#though stress is just one piece of the puzzle and shouldn't be the only thing considered

points <- data.frame(nmds$points) %>% 
  mutate(site_month = rownames(fish.list$abund)) %>% 
  separate_wider_delim(site_month, delim = "_", names = c("site", "month"), cols_remove = FALSE) %>% 
  mutate(site_type = ifelse(site %in% c("MA", "WA", "HO", "TL", "LL", "PR"), "jubilee", "core"))

hulls <- points %>%
  group_by(site) %>% 
  slice(chull(MDS1,MDS2))

#all
ggplot() +
  geom_point(data = points, aes(x = MDS1, y = MDS2, color = site, shape = site_type), size = 3) + 
  geom_text_repel(data = points %>% filter(!site %in% c("MA", "WA", "HO", "TL", "LL", "PR")), aes(x = MDS1, y = MDS2, color = site, label = month)) +
  geom_text_repel(data = points %>% filter(site %in% c("MA", "WA", "HO", "TL", "LL", "PR")), aes(x = MDS1, y = MDS2, color = site, label = site)) +
  geom_polygon(data = hulls, aes(x = MDS1, y = MDS2, fill = site), alpha = 0.2) +
  annotate("text", x = Inf, y = Inf, label = paste("stress = ", S), vjust = 2, hjust = 2) +
  annotate("text", x = Inf, y = Inf, label = paste("k = ", K), vjust = 2, hjust = 2) +
  theme_minimal() 

ggsave("docs/figures/fish_nmds_sites.png")

#core
ggplot() +
  geom_point(data = points %>% filter(!site %in% c("MA", "WA", "HO", "TL", "LL", "PR")), aes(x = MDS1, y = MDS2, color = site)) + 
  geom_text_repel(data = points %>% filter(!site %in% c("MA", "WA", "HO", "TL", "LL", "PR")), aes(x = MDS1, y = MDS2, color = site, label = month)) +
  geom_polygon(data = hulls %>% filter(!site %in% c("MA", "WA", "HO", "TL", "LL", "PR")), aes(x = MDS1, y = MDS2, fill = site), alpha = 0.2) +
  theme_minimal() 

# ggsave("docs/figures/nmds_core_sites.png")

############ follow up tests #################

#SIMPER
# fish.dist <- metaMDS(fish.adj.abund, distance = "bray", k = 3)
# 
# simper(fish.adj.abund,
#        points$site,
#        permutations = 999,
#        parallel = 1)
#not sure why I'm getting an error

#ISA
library(indicspecies)
set.seed(42)

fish.ISA <- multipatt(x = fish.adj.abund,
                      cluster = points$site,
                      duleg = TRUE)

summary(fish.ISA) #hm not all that informative??
