library(tidyverse)
library(here)
library(ggrepel)
library(vegan)

load(here("data", "fish.list.june.Rdata")) #object created in 02_tidy_data
load(here("data", "RLQ1.Rdata"))

fish.adj.abund <- fish.list.june[[2]] 
#need to add a dummy species because you can't calculate BC distances on a zero (SHR_04)
fish.adj.abund <- cbind(fish.adj.abund, dummy = rep(1, times = nrow(fish.adj.abund)))
fish.adj.abund <- sqrt(fish.adj.abund) #data transformation (sqrt) to avoid weighting by very abundant taxa

K <- 2 #number of dimensions, more than 3 generally gets unweildy
nmds <- metaMDS(fish.adj.abund, distance="bray", k= K, trymax=1000, plot = FALSE)
S <- round(nmds$stress, 2) #stress decreases with increasing K, >0.2 is generally a poor fit
#though stress is just one piece of the puzzle and shouldn't be the only thing considered

points <- data.frame(nmds$points) %>% mutate(site = rownames(fish.list.june$abund))

# hulls <- points %>%
#   slice(chull(MDS1,MDS2))

#all
ggplot() +
  geom_point(data = points, aes(x = MDS1, y = MDS2, color = RLQ1), size = 3) + 
  geom_text_repel(data = points, aes(x = MDS1, y = MDS2, color = RLQ1, label = site)) +
  # annotate("text", x = Inf, y = Inf, label = paste("stress = ", S), vjust = 2, hjust = 2) +
  # annotate("text", x = Inf, y = Inf, label = paste("k = ", K), vjust = 2, hjust = 2) +
  theme_minimal() 

ggsave("docs/figures/fish_nmds_sites.png")

#all by month
hulls2 <- points %>%
  group_by(month) %>% 
  slice(chull(MDS1,MDS2))

ggplot() +
  geom_point(data = points, aes(x = MDS1, y = MDS2, color = month, shape = site_type), size = 3) + 
  geom_text_repel(data = points, aes(x = MDS1, y = MDS2, color = month, label = site)) +
  geom_polygon(data = hulls2, aes(x = MDS1, y = MDS2, fill = month), alpha = 0.2) +
  annotate("text", x = Inf, y = Inf, label = paste("stress = ", S), vjust = 2, hjust = 2) +
  annotate("text", x = Inf, y = Inf, label = paste("k = ", K), vjust = 2, hjust = 2) +
  theme_minimal() 

ggsave("docs/figures/fish_nmds_months.png")

#core
ggplot() +
  geom_point(data = points %>% filter(!site %in% c("MA", "WA", "HO", "TL", "LL", "PR")), aes(x = MDS1, y = MDS2, color = site)) + 
  geom_text_repel(data = points %>% filter(!site %in% c("MA", "WA", "HO", "TL", "LL", "PR")), aes(x = MDS1, y = MDS2, color = site, label = month)) +
  geom_polygon(data = hulls %>% filter(!site %in% c("MA", "WA", "HO", "TL", "LL", "PR")), aes(x = MDS1, y = MDS2, fill = site), alpha = 0.2) +
  theme_minimal() 

# ggsave("docs/figures/nmds_core_sites.png")

############ follow up tests #################

#PERMANOVA
site_result <- adonis2(fish.adj.abund ~ points$site, method = "bray")
month_result <- adonis2(fish.adj.abund ~ points$month, method = "bray")


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
