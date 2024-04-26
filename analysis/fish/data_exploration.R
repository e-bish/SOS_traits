library(tidyverse)
library(here)
library(janitor)
library(vegan)
# library(lme4)
# library(sjPlot)

load_data <- function() {
  #Load data
  #load field data (raw data downloaded/formatted in "import_data.R")
  net_2018.19 <- here::here("data","net_18.19.csv") %>% 
    read_csv()
  
  net_2021 <- here::here("data","net_2021.csv") %>% 
    read_csv()
  
  net_2022 <- here::here("data","net_2022.csv") %>% 
    read_csv()
  
  net_tidy <- bind_rows(net_2018.19, net_2021, net_2022) %>% 
    filter(!ipa == "Armored_2") %>% #remove second armored site from Titlow in 2021
    mutate(ipa = replace(ipa, site == "TUR" & ipa == "Restored", "Natural")) %>% #no restoration at Turn Island
    mutate(month = if_else(site == "MA", "06", month)) %>% # we did a July 1st survey at Maylor that we want to count as a June survey
    # mutate(site_ipa = paste(site, ipa, sep = "_")) %>% 
    mutate(date = make_date(year, month, day)) %>% 
    mutate(species_count = as.numeric(species_count)) %>%
    mutate(species_count = ifelse(is.na(org_type), 0, species_count)) %>% 
    filter(org_type == "Fish") %>% #### this is where you lose zero counts, if those are needed later
    rename(ComName = "species") %>% 
    mutate(ComName = case_when(tax_group == "Salmon" ~ paste(ComName, "salmon"), 
                               ComName == "Sand Sole" ~ "Pacific sand sole",
                               ComName == "Herring" ~ "Pacific herring", 
                               ComName == "Snake Prickleback" ~ "Pacific snake prickleback",
                               ComName == "Speckled Sand Dab" ~ "Speckled sanddab",
                               ComName == "Poacher" ~ "UnID Poacher",
                               ComName == "Striped Perch" ~ "Striped Seaperch",
                               ComName == "Staghorn Sculpin" ~ "Pacific staghorn sculpin",
                               ComName == "Stickleback" ~ "Three-spined stickleback",
                               ComName == "Tom Cod" ~ "Pacific tomcod",
                               ComName == "Sand Lance" ~ "Pacific sand lance",
                               ComName == "Pipefish" ~ "Bay pipefish",
                               ComName == "Irish Lord" ~ "UnID Irish Lord",
                               ComName == "Midshipman" ~ "Plainfin midshipman",
                               ComName == "White Spotted Greenling" ~ "Whitespotted greenling",
                               ComName == "Silver Spotted Sculpin" ~ "Silverspotted sculpin",
                               ComName == "Tube Snout" ~ "Tube-snout",
                               TRUE ~ ComName)) %>% 
    mutate(ComName = case_when(ComName == "Steelhead salmon" ~ "Steelhead trout",
                               ComName == "Cutthroat salmon" ~ "Cutthroat trout",
                               TRUE ~ ComName)) 
  return(net_tidy)
}

net_tidy <- load_data()


#in how many sites does each species occur?
#what is the abundance of each species when it occurs?
#is the mean abundance correlated with the number of sites where it occurs?
#is the total abundance of fish correlated with the number of species in each site?

#species accumulation curves for each site
make_filtered_spp_mat <- function(site_ID) {

  spp_mat <- net_tidy %>% 
    mutate(year_month = ym(paste(year, month, sep = "-"))) %>% 
    filter(site == site_ID) %>% 
    group_by(ComName, year_month) %>% 
    summarize(spp_sum = sum(species_count)) %>% 
    pivot_wider(names_from = ComName, values_from = spp_sum) %>% 
    select(-contains("UnID")) %>% 
    arrange(year_month) %>% 
    clean_names() %>% 
    replace(is.na(.), 0) %>% 
    select(-1)
  
  return(spp_mat)
}

make_spp_curve <- function(site_ID) {

  site_mat <- make_filtered_spp_mat(site_ID)
  site_curve <- specaccum(site_mat, method = "collector", permutations = 100) #collector method preserves the order
  return(site_curve)

}

SOS_sites <- c("FAM", "TUR", "COR", "MA", "WA", "HO", "SHR", "DOK", "LL", "TL", "PR", "EDG")

curve_list <- lapply(SOS_sites, make_spp_curve)
curve_df <- data.frame()

for (i in 1:length(SOS_sites)) {
  sites <- curve_list[[i]]$sites
  richness <- curve_list[[i]]$richness
  tmp <- data.frame(site = SOS_sites[i], samples = sites, richness = richness)
  
  curve_df <- rbind(curve_df, tmp)
}
  
ggplot() +
  geom_point(data = curve_df, aes(x = samples, y = richness, color = factor(site, level = SOS_sites))) +
  geom_line(data = curve_df, aes(x = samples, y = richness, color = factor(site, level = SOS_sites))) +
  theme_classic() +
  labs(x = "Times sampled", y = "Number of species", color = "Site")

### rarefaction curve
# spp_mat <- make_spp_mat()
# S <- specnumber(spp_mat)
# (raremax <- min(rowSums(spp_mat)))
# Srare <- rarefy(spp_mat, raremax)
# plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
# rarecurve(spp_mat, step = 20, sample = raremax, col = "blue", cex = 0.6)

#for each species, how many sites were they encountered?
times_encountered <- net_tidy %>% 
  group_by(ComName, tax_group, site) %>% 
  summarize(freq_obs = n())
  
ggplot(times_encountered, aes(x = tax_group, y = freq_obs, fill = factor(site, levels = SOS_sites))) + 
  geom_bar(position = "stack", stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "", y = "Frequency observed", fill = "Site") 

ggplot(times_encountered, aes(x = factor(site, levels = SOS_sites), y = freq_obs, fill = tax_group)) + 
  geom_bar(position = "stack", stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "", y = "Frequency observed", fill = "tax group") 

sites_encountered <- times_encountered %>% 
  mutate(pa = ifelse(freq_obs > 0, 1, 0)) %>% 
  group_by(ComName) %>% 
  summarize(encounters_site = sum(pa)) %>% 
  arrange(desc(encounters_site))

ggplot(sites_encountered, aes(x = reorder(ComName, -encounters_site), y = encounters_site)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs(x = "", y = "No. of sites encountered")

#number of species by site
spp_by_site <- net_tidy %>% 
  select(-contains("UnID")) %>% 
  group_by(site) %>% 
  summarize(n_spp = n_distinct(ComName)) %>% 
  ungroup()

#beta diversity by site
make_spp_mat <- function(pa) {
  
  spp_mat <- net_tidy %>% 
    mutate(site_month = paste(site, month, sep = "_")) %>% 
    group_by(site_month, ComName) %>% 
    summarize(spp_sum = sum(species_count)) %>% 
    ungroup() %>% 
    pivot_wider(names_from = ComName, values_from = spp_sum) %>% 
    select(-contains("UnID")) %>% 
    clean_names() %>% 
    replace(is.na(.), 0) %>% 
    select(-1) 
  
  spp_pa_mat <- decostand(spp_mat, method = "pa") #change from abundance to presence absence
  
  if (pa == T) {
    return(spp_pa_mat)
  } else {
    return(spp_mat)
  }
}

spp_mat <- make_spp_mat(pa = T)
rownames(spp_mat) <- SOS_sites

spp_dist_mat <- betadiver(spp_mat, method = "sor") #sorensen dissimilarity index that accounts for the fact that you inherently see more spp with larger sample sizes
#but actually need a different index for PA data??!
plot(spp_dist_mat)
str(spp_dist_mat)
metaMDS(spp_dist_mat)


spp_b_disp <- betadisper(spp_dist_mat, group = SOS_sites)
spp_b_disp

adonis2(spp_dist_mat~net_tidy$site)

anova(spp_b_disp)
plot(spp_b_disp)
boxplot(spp_b_disp)

permutest(spp_b_disp, pairwise = TRUE, permutations = 99)


#species frequency by site
spp_freq <- net_tidy %>% 
  complete(nesting(year, month, site), station, ipa) %>% 
  filter(!((site == "MA" & year == 2021 & ipa == "Armored") |
             (site == "SHR" & year == 2021 & month == "08" & ipa == "Armored")) ) %>% # no surveys on these days
  #group_by() something?? CPUE??


#kruskal wallace non parametric test
#spearman rank correlation (urbanization to spp richness)

## species richness
spp_richness <- net_tidy %>% 
  mutate(species_count = ifelse(is.na(species_count), 0, species_count)) %>% 
  group_by(year, month, site, ipa, station) %>% 
  summarize(richness = n_distinct(ComName, na.rm = FALSE)) %>% 
  ungroup() %>% 
  complete(nesting(year, month, site, ipa), station) %>% #use nesting to only include combinations already present in the data
  replace(is.na(.),0) %>% 
  mutate(ipa = factor(ipa, levels = c("Natural", "Armored", "Restored")), 
         station = factor(station, labels = c("shallow", "mid", "deep"))) %>% 
  mutate(season = ifelse(month %in% c("04", "05", "06"), "spring", "summer"))

ggplot(spp_richness) + 
  geom_violin(aes(x = station, y = richness))

ggplot(spp_richness) + 
  geom_violin(aes(x = ipa, y = richness))

ggplot(spp_richness) + 
  geom_violin(aes(x = site, y = richness))

ggplot(spp_richness) + 
  geom_violin(aes(x = season, y = richness))



