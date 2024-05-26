library(tidyverse)
library(here)
library(janitor)
library(vegan)

net_2018.19 <- here::here("data","net_18.19.csv") %>% 
  read_csv()

net_2021 <- here::here("data","net_2021.csv") %>% 
  read_csv()

net_2022 <- here::here("data","net_2022.csv") %>% 
  read_csv()

SOS_sites <- c("FAM", "TUR", "COR", "MA", "WA", "HO", "SHR", "DOK", "LL", "TL", "PR", "EDG")
SOS_core_sites <- factor(c("FAM", "TUR", "COR", "SHR", "DOK", "EDG"), levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG"))

##### tidy the imported data

net_tidy <- bind_rows(net_2018.19, net_2021, net_2022) %>% 
  filter(!ipa == "Armored_2") %>% #remove second armored site from Titlow in 2021
  mutate(ipa = replace(ipa, site == "TUR" & ipa == "Restored", "Natural")) %>% #no restoration at Turn Island
  mutate(month = if_else(site == "MA", "06", month)) %>% # we did a July 1st survey at Maylor that we want to count as a June survey
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
                             TRUE ~ ComName)) %>% 
  filter(!grepl("UnID", ComName)) %>% #remove unidentified species
  select(!c(day, org_type, tax_group, mean_length_mm)) %>% 
  mutate(month = recode(month, `04` = "Apr", `05` = "May", `06` = "Jun", `07` = "Jul", `08` = "Aug", `09` = "Sept")) %>% 
  mutate(year = factor(year), month = factor(month, levels = c("Apr", "May", "Jun", "Jul", "Aug", "Sept")))

#### create alternative abundance matrices
#retain 4 years of data at the core sites
Core_4years <- net_tidy %>% 
  filter(site %in% SOS_core_sites) %>% #filter to core sites
  group_by(year, month, site, ComName) %>% 
  summarize(species_count = sum(species_count)) %>% #sum across shorelines and stations, still need to find a way to account for sampling effort at times when we missed a station
  pivot_wider(names_from = ComName, values_from = species_count, values_fill = 0) %>% 
  clean_names() 

Core_4years.mat <- sqrt(Core_4years[4:45]) #square root transformation for count data

#this one would be paired with the 2022 bird data  
Core_1year <- net_tidy %>% 
  filter(site %in% SOS_core_sites & year == 2022) %>% 
  group_by(month, site, ComName) %>% 
  summarize(species_count = sum(species_count)) %>% #sum across shorelines and stations, need to account for sampling effort?
  pivot_wider(names_from = ComName, values_from = species_count, values_fill = 0) %>% 
  clean_names()

Core_1year.mat <- sqrt(Core_1year[3:32]) #square root transformation for count data

#this one would be paired with the two years of bird data (though 2022 is more complete)
All_2years <- net_tidy %>% 
  filter(year == 2021 | year == 2022 & month == "Jun") %>% #all sites but only june in the last two years when we sampled 12
  group_by(year, site, ComName) %>% 
  summarize(species_count = sum(species_count)) %>%  #need to account for sampling effort?
  pivot_wider(names_from = ComName, values_from = species_count, values_fill = 0) %>%
  clean_names()

All_2years.mat <- sqrt(All_2years[3:35]) #square root transformation for count data

#### NMDS 4 years core
Core_4years.mds <- metaMDS(Core_4years.mat) 
Core_4years.scores <- scores(Core_4years.mds)
Core_4years.site.scores <- as.data.frame(Core_4years.scores$sites)

ggplot(data = Core_4years.site.scores) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Core_4years$site), level = 0.5) +
  geom_point(aes(x = NMDS1, y = NMDS2, color = Core_4years$site), size = 4) + 
  theme_classic()

ggplot(data = Core_4years.site.scores) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Core_4years$year), level = 0.5) +
  geom_point(aes(x = NMDS1, y = NMDS2, color = Core_4years$year), size = 4) + 
  theme_classic()

ggplot(data = Core_4years.site.scores) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Core_4years$month), level = 0.5) +
  geom_point(aes(x = NMDS1, y = NMDS2, color = Core_4years$month), size = 4) + 
  theme_classic()

#### PERMANOVA 4 years core
Core_4years.perm <- adonis2(Core_4years.mat ~ Core_4years$site + Core_4years$year + Core_4years$month, method = "bray", perm = 999)
print(Core_4years.perm) #significant differences between sites, years, and months
#but this isn't correct because it isn't restricting permutations??
#don't need to restrict between depths and stations because we summed across them
#what would we restrict?
# between years? months?

Core_4years.aov <- aov(Core_4years ~ site + month + Error(year), data = Core_4years)


#### NMDS 1 year core
Core_1year.mds <- metaMDS(Core_1year.mat) 
Core_1year.scores <- scores(Core_1year.mds)
Core_1year.site.scores <- as.data.frame(Core_1year.scores$sites)

ggplot(data = Core_1year.site.scores) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Core_1year$site), level = 0.5) +
  geom_point(aes(x = NMDS1, y = NMDS2, color = Core_1year$site), size = 4) + 
  theme_classic()

ggplot(data = Core_1year.site.scores) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Core_1year$month), level = 0.5) +
  geom_point(aes(x = NMDS1, y = NMDS2, color = Core_1year$month), size = 4) + 
  theme_classic()

#### PERMANOVA 1 year core
Core_1year.perm <- adonis2(Core_1year.mat ~ Core_1year$site + Core_1year$month, method = "bray", perm = 999)
print(Core_1year.perm) #significant differences between sites, years, and months

#### NMDS all sites June 2 years
All_2years.mds <- metaMDS(All_2years.mat)
All_2years.scores <- scores(All_2years.mds)
All_2years.site.scores <- as.data.frame(All_2years.scores$sites)

u_group <- All_2years %>% 
  mutate(armor_group = case_when(site %in% c("FAM", "TUR", "WA", "EDG") ~ "low_urbanization",
                                 site %in% c("COR", "DOK", "PR", "LL") ~ "med_urbanization",
                                 site %in% c("HO", "SHR", "MA", "TL") ~ "high_urbanization")) %>% #one idea for post-hoc groupings
  ungroup() %>% 
  mutate(armor_group = factor(armor_group, levels = c("low_urbanization", "med_urbanization", "high_urbanization"))) %>% 
  pull(armor_group)

ggplot(data = All_2years.site.scores) +
  #stat_ellipse(aes(x = NMDS1, y = NMDS2, color = All_2years$site), level = 0.5) + #too few points
  geom_point(aes(x = NMDS1, y = NMDS2, color = u_group), size = 4) + 
  theme_classic() 






