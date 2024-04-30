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
  mutate(sample = rownames(fish.list$abund)) %>% 
  separate_wider_delim(sample, delim = "_", names = c("year", "month", "site"), cols_remove = FALSE) 

hulls <- points %>%
  group_by(site) %>% 
  slice(chull(MDS1,MDS2))

#all
ggplot() +
  geom_point(data = points, aes(x = MDS1, y = MDS2, color = site, shape = year), size = 3) + 
  geom_text_repel(data = points, aes(x = MDS1, y = MDS2, color = site, label = month)) +
  geom_polygon(data = hulls, aes(x = MDS1, y = MDS2, fill = site), alpha = 0.2) +
  annotate("text", x = Inf, y = Inf, label = paste("stress = ", S), vjust = 2, hjust = 2) +
  annotate("text", x = Inf, y = Inf, label = paste("k = ", K), vjust = 2, hjust = 2) +
  theme_minimal() 

ggsave("docs/figures/fish_nmds_sites.png")

#all by month
hulls2 <- points %>%
  group_by(month) %>% 
  slice(chull(MDS1,MDS2))

ggplot() +
  geom_point(data = points, aes(x = MDS1, y = MDS2, color = month, shape = year), size = 3) + 
  geom_text_repel(data = points, aes(x = MDS1, y = MDS2, color = month, label = site)) +
  geom_polygon(data = hulls2, aes(x = MDS1, y = MDS2, fill = month), alpha = 0.2) +
  annotate("text", x = Inf, y = Inf, label = paste("stress = ", S), vjust = 2, hjust = 2) +
  annotate("text", x = Inf, y = Inf, label = paste("k = ", K), vjust = 2, hjust = 2) +
  theme_minimal() 

ggsave("docs/figures/fish_nmds_months.png")

############ follow up tests #################

#PERMANOVA
site_result <- adonis2(fish.adj.abund ~ points$site + points$month + points$year, method = "bray")
site_result


#SIMPER
# fish.dist <- metaMDS(fish.adj.abund, distance = "bray", k = 3)
# 
# simper(fish.adj.abund,
#        points$site,
#        permutations = 999,
#        parallel = 1)
#not sure why I'm getting an error

#ISA
# library(indicspecies)
# set.seed(42)
# 
# fish.ISA <- multipatt(x = fish.adj.abund,
#                       cluster = points$site,
#                       duleg = TRUE)
# 
# summary(fish.ISA) #hm not all that informative??
