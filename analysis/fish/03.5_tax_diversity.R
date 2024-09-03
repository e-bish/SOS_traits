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
beta.multi(fish_L_pa, index.family="jaccard")
beta.pair(fish_L_pa, index.family="jaccard")




#### taxonomic diversity with raw counts ####

fish_counts <- net_core %>% 
  filter(!is.na(ComName)) %>% 
  mutate(year_month = paste(year, month, sep = "_")) %>% 
  group_by(site, ComName, year_month) %>% 
  summarize(spp_sum = sum(species_count)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = ComName, values_from = spp_sum, values_fill = 0) %>% 
  arrange(site, year_month) %>% 
  clean_names() %>% 
  mutate(sample = paste(site, year_month, sep = "_"), .before = site) %>% 
  select(!2:3) %>% 
  ungroup() %>% 
  column_to_rownames("sample") 

fish_counts_long <- fish_counts %>% 
  rownames_to_column("sample") %>% 
  pivot_longer(!sample, names_to = "species", values_to = "count")

alpha_div_counts <- fish_counts_long %>% 
  group_by(sample) %>% 
  summarize(richness = specnumber(count),
            shannon = diversity(count, index = "shannon"),
            simpson = diversity(count, index = "simpson"),
            invsimpson = diversity(count, index = "invsimpson"),
            n = sum(count)) %>% 
  ungroup() %>% 
  separate_wider_delim(sample, delim = "_", names = c("site", "year", "month"), cols_remove = FALSE) %>% 
  mutate(site = factor(site, levels = SOS_core_sites))

alpha_div_counts %>% 
  pivot_longer(cols = c(richness, shannon, invsimpson, simpson), names_to = "metric") %>% 
  mutate(metric = factor(metric, levels = c("richness", "shannon", "simpson", "invsimpson"))) %>% 
  ggplot(aes(x = site, y = value, color = site)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(~metric, scales = "free_y") +
  theme_classic()


