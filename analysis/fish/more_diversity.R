
##did we sample enough to adequately describe the community?
make_filtered_spp_mat <- function(site_ID) {
  
  spp_mat <- net_core %>% 
    filter(!is.na(ComName)) %>% 
    mutate(year_month = ym(paste(year, month, sep = "-"))) %>% 
    filter(site == site_ID) %>% 
    group_by(ComName, year_month) %>% 
    summarize(spp_sum = sum(species_count)) %>% 
    pivot_wider(names_from = ComName, values_from = spp_sum) %>% 
    arrange(year_month) %>% 
    clean_names() %>% 
    replace(is.na(.), 0) %>% 
    select(-1)
  
  return(spp_mat)
}

make_spp_curve <- function(site_ID) {
  
  site_mat <- make_filtered_spp_mat(site_ID)
  site_curve <- specaccum(site_mat, method = "rarefaction", permutations = 100) #collector method preserves the order
  return(site_curve)
  
}

curve_list <- lapply(SOS_core_sites, make_spp_curve)
curve_df <- data.frame() 

for (i in 1:length(SOS_core_sites)) {
  sites <- curve_list[[i]]$sites
  richness <- curve_list[[i]]$richness
  tmp <- data.frame(site = SOS_core_sites[i], samples = sites, richness = richness)
  
  curve_df <- rbind(curve_df, tmp)
}

ggplot() +
  # geom_point(data = curve_df, aes(x = samples, y = richness, color = factor(site, level = SOS_sites))) +
  geom_line(data = curve_df, aes(x = samples, y = richness, color = factor(site, level = SOS_sites))) +
  theme_classic() +
  labs(x = "Times sampled", y = "Number of species", color = "Site")

### rarefy 

spp_mat <- net_core %>% 
    filter(!is.na(ComName)) %>% 
    mutate(year_month = ym(paste(year, month, sep = "-"))) %>% 
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
  
rarefy(spp_mat, sample = 1) #how many species you should expect to find for a given level of effort

rarecurve(spp_mat, sample = 1,tidy = TRUE) 


### other indices
shared <- net_core %>% 
  filter(!is.na(ComName)) %>% 
  group_by(site, ComName) %>% 
  tally()
  
shared %>% 
  group_by(site) %>% 
  summarize(richness = specnumber(n),
            shannon = diversity(n, index = "shannon"),
            simpson = diversity(n, index = "simpson"),
            invsimpson = 1/simpson,
            n = sum(n)) %>% 
  pivot_longer(cols = c(richness, shannon, invsimpson, simpson), names_to = "metric") %>% 
  ggplot(aes(x = n, y = value)) +
  geom_point() +
  geom_smooth() + 
  facet_wrap(~metric, nrow=4,scales="free_y")

shared_2 <- spp_mat %>% 
  decostand(method = "pa") %>% 
  mutate(sum = rowSums(across(where(is.numeric)))) %>% 
  rownames_to_column("sample") %>% 
  select(sample, sum)

shared_2 <- net_core %>% 
  filter(!is.na(ComName)) %>% 
  group_by(site, date, ComName) %>% 
  summarize(n = sum(species_count))

alpha_div <- shared_2 %>% 
  group_by(site, date) %>% 
  summarize(richness = specnumber(n),
            shannon = diversity(n, index = "shannon"),
            simpson = diversity(n, index = "simpson"),
            invsimpson = 1/simpson,
            n = sum(n))

alpha_div %>% 
  pivot_longer(cols = c(richness, shannon, invsimpson, simpson), names_to = "metric") %>% 
  ggplot(aes(x = n, y = value, color = site)) +
  geom_point() +
  geom_smooth() + 
  facet_wrap(~metric, nrow=4,scales="free_y")

#watch this: https://www.youtube.com/watch?v=ht3AX5uZTTQ&t=439s
### format matrix to have the number species in a sample
no_spp_in_sample <- net_core %>% 
  filter(!is.na(ComName)) %>% 
  select(site, year, month, ComName, species_count) %>% 
  group_by(site, year, month, ComName) %>% 
  summarize(sample_count = sum(species_count)) %>% #sum at each site on each sample day
  ungroup() %>% 
  group_by(site, year, month) %>% 
  mutate(n_individuals = sum(sample_count)) %>% 
  ungroup() %>% 
  filter(n_individuals >=3) %>% #only keep sample days where we caught at least 3 fish
  select(-n_individuals) 

abu_mat <- no_spp_in_sample %>% 
  pivot_wider(names_from = "ComName", values_from = "sample_count", values_fill = 0) %>% 
  mutate(sample = paste(site, month, year, sep = "_"), .before = site)

abu_meta <- abu_mat %>% 
  select(sample, site, year, month) %>% 
  mutate(season = ifelse(month %in% c("May", "Jun", "Jul"), "peak", "shoulder"))

abu_tbl <- abu_mat %>% 
  select(-c(site, year, month)) %>% 
  column_to_rownames("sample")

### rarefy
(raremax <- min(rowSums(abu_tbl)))
abu_dist <- avgdist(abu_tbl, sample = raremax)

set.seed(0416)
abu_nmds <- metaMDS(abu_dist) %>% 
  scores() %>% 
  as_tibble(rownames = "sample")

fish_nmds<- inner_join(abu_meta, abu_nmds, by = "sample")

centroids <- fish_nmds %>% 
  group_by(site) %>% 
  summarize(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2))

ggplot(fish_nmds, aes(x = NMDS1, y = NMDS2, color = site)) +
  geom_point() + 
  stat_ellipse(show.legend = FALSE) +
  theme_classic() +
  geom_point(data = centroids, size = 5, color = "black", shape = 21, 
             aes(fill = site), show.legend = FALSE)

adonis2(abu_dist~abu_meta$site)

season_centroids <- fish_nmds %>% 
  group_by(season) %>% 
  summarize(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2))

ggplot(fish_nmds, aes(x = NMDS1, y = NMDS2, color = season)) +
  geom_point() + 
  stat_ellipse(show.legend = FALSE) +
  theme_classic() +
  geom_point(data = season_centroids, size = 5, color = "black", shape = 21, 
             aes(fill = season), show.legend = FALSE)

adonis2(abu_dist~abu_meta$season, strata = abu_meta$site)

adonis2(abu_dist~abu_meta$site*abu_meta$season)

bd_site <- betadisper(abu_dist, abu_meta$site)
bd_season <- betadisper(abu_dist, abu_meta$season)

permutest(bd_site)
permutest(bd_season)
#no statistically significant disperson, so we can be confident that's not creating a false positive in the adonis test

### rarefy abundance matrix 

n_iter <- 100
alpha_list <- list()

for (i in 1:n_iter) {
  
  set.seed(i)
  
  abu_rare <- rrarefy(abu_tbl, sample = 3) %>% #tells you how many unique species you'd expect to see if you caught 3 fish
    as.data.frame() %>% 
    rownames_to_column("sample")
    
  abu_rare_tbl <- inner_join(abu_meta, abu_rare) %>% 
    pivot_longer(cols = 6:47, names_to = "ComName", values_to = "rare_n")
  
  ### calculate alpha diversity indices
  alpha_div <- abu_rare_tbl %>% 
    group_by(site) %>% 
    summarize(richness = specnumber(rare_n),
              shannon = diversity(rare_n, index = "shannon"),
              simpson = diversity(rare_n, index = "simpson"),
              invsimpson = 1/simpson,
              n = sum(rare_n))
  
  alpha_list[[i]] <- alpha_div
}

alpha_results <- alpha_list %>% 
  map_df(as_tibble) 
  
alpha_results %>% 
  select(!n) %>% 
  pivot_longer(!site, names_to = "metric") %>% 
  ggplot(aes(x = site, y = value)) +
    geom_violin()+
    theme_classic() +
    facet_wrap(~metric, scales = "free_y")
  

### calculate beta diversity indices

#do randomized rarefaction approach for calculating alpha diversity indices because you don't know the right number of times to sample to capture the species