library(tidyverse)
library(vegan)

#### same data aggregation procedure as the FD analysis ####
load(here("data", "fish.list.Rdata")) #object created in 03_create_matrices
SOS_core_sites <- factor(c("FAM", "TUR", "COR", "SHR", "DOK", "EDG"), 
                         levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG"))

fish_L <- fish.list$abund
fish_L_long <- fish_L %>% 
  rownames_to_column("sample") %>% 
  pivot_longer(!sample, names_to = "species", values_to = "avg_n")

#### alpha diversity ####
alpha_div <- fish_L_long %>% 
  group_by(sample) %>% 
  summarize(richness = specnumber(avg_n),
            shannon = diversity(avg_n, index = "shannon"),
            simpson = diversity(avg_n, index = "simpson"),
            invsimpson = diversity(avg_n, index = "invsimpson"),
            n = sum(avg_n)) %>% 
  ungroup() %>% 
  separate_wider_delim(sample, delim = "_", names = c("site", "ipa"), cols_remove = FALSE) %>% 
  mutate(site = factor(site, levels = SOS_core_sites))
  
alpha_div %>% 
  pivot_longer(cols = c(richness, shannon, invsimpson, simpson), names_to = "metric") %>% 
  mutate(metric = factor(metric, levels = c("richness", "shannon", "simpson", "invsimpson"))) %>% 
  ggplot(aes(x = site, y = value, color = site)) +
    geom_point(size = 2) +
    facet_wrap(~metric, scales = "free_y") +
    theme_classic()

#### beta diversity ####
library(betapart)

fish_bcdist <- vegdist(fish_L, method = "bray")
fish_meta <- fish_L %>% 
  rownames_to_column("sample") %>% 
  separate_wider_delim(sample, delim = "_", names = c("site", "ipa"), cols_remove = FALSE) %>% 
  select(1:3)

adonis2(fish_bcdist~fish_meta$site) #test for differences in means
bd_site <- betadisper(fish_bcdist, fish_meta$site)
permutest(bd_site) #no significant dispersion differences, so we can trust the adonis output

fish_L_pa <- decostand(fish_L, method = "pa")
beta.multi(fish_L_pa)
beta.pair(fish_L_pa)
