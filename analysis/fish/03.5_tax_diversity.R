library(tidyverse)
library(vegan)

## rarefy? can only rarefy count data

#load tidy fish data frame created in 02_tidy_data
load(here("data", "net_tidy.Rdata")) 

#### Create the L (abundance) matrix ####
expand_species <- net_tidy %>% 
  expand(nesting(year, month, site, ipa), ComName) %>% 
  filter(!is.na(ComName))

fish_counts <- net_tidy %>% #L is referring to the RLQ analysis
  filter(!is.na(ComName)) %>% 
  group_by(year, month, site, ipa, ComName) %>%
  summarize(spp_sum = sum(species_count)) %>% #sum across depths within a site
  ungroup() %>%
  full_join(expand_species) %>% #add back in all of the events so we can capture 0s
  mutate(ComName = replace(ComName, ComName == "Pacific Sandfish", "Pacific sandfish")) %>% #if it's capitolized it gets confused about the proper order
  mutate(spp_sum =replace_na(spp_sum, 0)) %>% 
  group_by(year, month, site, ComName) %>% 
  summarize(spp_sum_all = sum(spp_sum))  %>% 
  ungroup()

fish_counts_mat <- fish_counts %>% 
  pivot_wider(id_cols = !c(year, month, site), names_from = ComName, values_from = spp_sum_all)
#why isn't this working??

  group_by(site, year, ComName) %>% 
  summarize(spp_avg = mean(spp_sum)) %>% #average across sampling events at each site/ipa (unbalanced)
  pivot_wider(names_from = ComName, values_from = spp_avg, values_fill = 0) %>% 
  clean_names() %>% 
  ungroup() %>% 
  mutate(sample = paste(site, year, sep = "_"), .after = year) %>% 
  select(!1:2) %>% 
  column_to_rownames(var = "sample")


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
    geom_boxplot() +
    geom_point() +
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


