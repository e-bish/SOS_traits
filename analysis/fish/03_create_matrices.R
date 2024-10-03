library(tidyverse)
library(here)
library(rfishbase)
library(janitor)
library(GGally)

#load tidy fish data frame created in 02_tidy_data
load(here("data", "net_tidy.Rdata")) 


## sum across the year ##
# sampling_events <- net_tidy %>% 
#   select(year, month, day, site, ipa, station) %>% 
#   distinct() %>% 
#   filter(!(site == "COR" & year == "2019" & month == "Apr" & day == "30")) %>% #these COR seem like data entry errors because they were only sampled at one station and we already had complete sampling at COR in april and may. No fish recorded in either entry
#   filter(!(site == "COR" & year == "2019" & month == "May" & day == "01")) %>% 
#   filter(!(site == "TUR" & year == "2018" & month == "Jul" & day == "11")) %>%  #this is an incomplete sampling event. Complete sampling at TUR occured on 7/12/18
#   filter(!(site == "TUR" & year == "2018" & month == "Sept" & day == "11")) %>% #another month with repeat sampling; keeping only the second september TUR sampling event
#   group_by(year, site, ipa) %>% 
#   summarize(no_net_sets = n()) %>% 
#   ungroup()
# 
# sampling_events_ipa <- sampling_events %>% 
#   group_by(site, ipa) %>% 
#   summarize(total_net_sets = sum(no_net_sets)) %>% 
#   ungroup()
# 
# sampling_events_yr <- sampling_events %>% 
#   group_by(year, site) %>% 
#   summarize(total_net_sets = sum(no_net_sets)) %>% 
#   ungroup()

## average across the year ##

sampling_events <- net_tidy %>% 
  select(year, month, day, site, ipa, station) %>% 
  distinct() %>% 
  filter(!(site == "COR" & year == "2019" & month == "Apr" & day == "30")) %>% #these COR seem like data entry errors because they were only sampled at one station and we already had complete sampling at COR in april and may. No fish recorded in either entry
  filter(!(site == "COR" & year == "2019" & month == "May" & day == "01")) %>% 
  filter(!(site == "TUR" & year == "2018" & month == "Jul" & day == "11")) %>%  #this is an incomplete sampling event. Complete sampling at TUR occured on 7/12/18
  filter(!(site == "TUR" & year == "2018" & month == "Sept" & day == "11")) %>% #another month with repeat sampling; keeping only the second september TUR sampling event
  group_by(year, month, site, ipa) %>% 
  summarize(no_net_sets = n()) %>% 
  ungroup()

sampling_events_ipa <- sampling_events %>%
  group_by(year, month, site, ipa) %>%
  summarize(total_net_sets = sum(no_net_sets)) %>%
  ungroup()
# 
# sampling_events_yr <- sampling_events %>% 
#   group_by(year, site) %>% 
#   summarize(total_net_sets = sum(no_net_sets)) %>% 
#   ungroup()

sampling_events_mo.yr <- sampling_events %>%
  group_by(year, month, site) %>%
  summarize(total_net_sets = sum(no_net_sets)) %>%
  ungroup()

#### Create the L (abundance) matrix ####

## sum across the year ##
expand_species_ipa <- net_tidy %>%
  expand(nesting(site, ipa), ComName) %>%
  filter(!is.na(ComName))
# 
# expand_species_yr <- net_tidy %>% 
#   expand(nesting(year, site), ComName) %>% 
#   filter(!is.na(ComName))
# 
# expand_species_both <- net_tidy %>% 
#   expand(nesting(year, site, ipa), ComName) %>% 
#   filter(!is.na(ComName))
# 
# fish_L.ipa <- net_tidy %>% #L is referring to the RLQ analysis
#   filter(!is.na(ComName)) %>% 
#   group_by(site, ipa, ComName) %>%
#   summarize(spp_sum = sum(species_count)) %>% #sum across net sets within a shoreline by year
#   ungroup() %>%
#   full_join(sampling_events_ipa) %>% 
#   mutate(catch_per_set = spp_sum/total_net_sets) %>% 
#   full_join(expand_species_ipa) %>% #add back in all of the events so we can capture 0s
#   mutate(ComName = replace(ComName, ComName == "Pacific Sandfish", "Pacific sandfish")) %>% #if it's capitolized it gets confused about the proper order
#   arrange(site, ipa, ComName) %>% 
#   mutate(catch_per_set =replace_na(catch_per_set, 0)) %>% 
#   select(!c(spp_sum, total_net_sets)) %>% 
#   pivot_wider(names_from = ComName, values_from = catch_per_set, values_fill = 0) %>% 
#   clean_names() %>% 
#   ungroup() %>% 
#   mutate(sample = paste(site, ipa, sep = "_"), .after = ipa) %>% 
#   select(!1:2) %>% 
#   column_to_rownames(var = "sample")
# 
# fish_L.year <- net_tidy %>% #L is referring to the RLQ analysis
#   filter(!is.na(ComName)) %>% 
#   group_by(year, site, ComName) %>% 
#   summarize(spp_sum = sum(species_count)) %>% #sum across sampling events 
#   ungroup() %>%
#   full_join(sampling_events_yr) %>% 
#   mutate(catch_per_set = spp_sum/total_net_sets) %>% 
#   full_join(expand_species_yr) %>% #add back in all of the events so we can capture 0s
#   mutate(ComName = replace(ComName, ComName == "Pacific Sandfish", "Pacific sandfish")) %>% #if it's capitolized it gets confused about the proper order
#   arrange(year, site, ComName) %>% 
#   mutate(catch_per_set  =replace_na(catch_per_set , 0)) %>% 
#   select(!c(spp_sum, total_net_sets)) %>% 
#   pivot_wider(names_from = ComName, values_from = catch_per_set, values_fill = 0) %>% 
#   clean_names() %>% 
#   ungroup() %>% 
#   mutate(sample = paste(site, year, sep = "_"), .after = site) %>% 
#   select(!1:2) %>% 
#   column_to_rownames(var = "sample")
# 
# fish_L.both <- net_tidy %>% #L is referring to the RLQ analysis
#   filter(!is.na(ComName)) %>% 
#   group_by(year, site, ipa, ComName) %>% 
#   summarize(spp_sum = sum(species_count)) %>% #sum across sampling events 
#   ungroup() %>%
#   full_join(sampling_events) %>% 
#   mutate(catch_per_set = spp_sum/no_net_sets) %>% 
#   full_join(expand_species_both) %>% #add back in all of the events so we can capture 0s
#   mutate(ComName = replace(ComName, ComName == "Pacific Sandfish", "Pacific sandfish")) %>% #if it's capitolized it gets confused about the proper order
#   arrange(year, site, ipa, ComName) %>% 
#   mutate(catch_per_set  = replace_na(catch_per_set , 0)) %>% 
#   select(!c(spp_sum, no_net_sets)) %>% 
#   pivot_wider(names_from = ComName, values_from = catch_per_set, values_fill = 0) %>% 
#   clean_names() %>% 
#   ungroup() %>% 
#   mutate(sample = paste(site, year, ipa, sep = "_"), .after = ipa) %>% 
#   select(!1:3) %>% 
#   column_to_rownames(var = "sample")
# 
# check_matrix <- fish_L.both %>% 
#   decostand(method = "pa") %>% 
#   mutate(no_spp = rowSums(.)) %>% 
#   filter(no_spp < 3)
# 
# fish_L.both <- fish_L.both %>% 
#   filter(!row.names(.) %in% row.names(check_matrix))

## average across the year ##
# expand_species_ipa <- net_tidy %>% 
#   expand(nesting(site, ipa), ComName) %>% 
#   filter(!is.na(ComName))
# 
expand_species_mo.yr <- net_tidy %>%
  expand(nesting(year, month, site), ComName) %>%
  filter(!is.na(ComName))

# expand_species_both <- net_tidy %>% 
#   expand(nesting(year, site, ipa), ComName) %>% 
#   filter(!is.na(ComName))

# doesn't need to be updated for averaging?
fish_L.ipa <- net_tidy %>% #L is referring to the RLQ analysis
  filter(!is.na(ComName)) %>%
  group_by(site, ipa, ComName) %>%
  summarize(spp_sum = sum(species_count)) %>% #sum across net sets within a shoreline
  ungroup() %>%
  full_join(sampling_events_ipa) %>%
  mutate(catch_per_set = spp_sum/total_net_sets) %>%
  full_join(expand_species_ipa) %>% #add back in all of the events so we can capture 0s
  mutate(ComName = replace(ComName, ComName == "Pacific Sandfish", "Pacific sandfish")) %>% #if it's capitolized it gets confused about the proper order
  arrange(site, ipa, ComName) %>%
  mutate(catch_per_set =replace_na(catch_per_set, 0)) %>%
  select(!c(spp_sum, total_net_sets)) %>%
  pivot_wider(names_from = ComName, values_from = catch_per_set, values_fill = 0) %>%
  clean_names() %>%
  ungroup() %>%
  mutate(sample = paste(site, ipa, sep = "_"), .after = ipa) %>%
  select(!1:2) %>%
  column_to_rownames(var = "sample")

fish_L.year <- net_tidy %>% #L is referring to the RLQ analysis
  filter(!is.na(ComName)) %>% 
  group_by(year, month, site, ComName) %>% 
  summarize(spp_sum = sum(species_count)) %>% #sum across sampling events 
  ungroup() %>%
  full_join(sampling_events_mo.yr) %>% 
  mutate(catch_per_set = spp_sum/total_net_sets) %>% 
  filter(!is.na(ComName)) %>% #these are accounted for in the next step
  full_join(expand_species_mo.yr) %>% #add back in all of the events so we can capture 0s
  mutate(ComName = replace(ComName, ComName == "Pacific Sandfish", "Pacific sandfish")) %>% #if it's capitolized it gets confused about the proper order
  arrange(year, month, site, ComName) %>% 
  mutate(catch_per_set  = replace_na(catch_per_set , 0)) %>% 
  group_by(year, site, ComName) %>% 
  summarize(avg_cps = mean(catch_per_set)) %>% #average across months within a year
  pivot_wider(names_from = ComName, values_from = avg_cps, values_fill = 0) %>% 
  clean_names() %>% 
  ungroup() %>% 
  mutate(sample = paste(site, year, sep = "_"), .after = site) %>% 
  select(!1:2) %>% 
  column_to_rownames(var = "sample")

# haven't updated this for the averaging yet
# fish_L.both <- net_tidy %>% #L is referring to the RLQ analysis
#   filter(!is.na(ComName)) %>% 
#   group_by(year, site, ipa, ComName) %>% 
#   summarize(spp_sum = sum(species_count)) %>% #sum across sampling events 
#   ungroup() %>%
#   full_join(sampling_events) %>% 
#   mutate(catch_per_set = spp_sum/no_net_sets) %>% 
#   full_join(expand_species_both) %>% #add back in all of the events so we can capture 0s
#   mutate(ComName = replace(ComName, ComName == "Pacific Sandfish", "Pacific sandfish")) %>% #if it's capitolized it gets confused about the proper order
#   arrange(year, site, ipa, ComName) %>% 
#   mutate(catch_per_set  = replace_na(catch_per_set , 0)) %>% 
#   select(!c(spp_sum, no_net_sets)) %>% 
#   pivot_wider(names_from = ComName, values_from = catch_per_set, values_fill = 0) %>% 
#   clean_names() %>% 
#   ungroup() %>% 
#   mutate(sample = paste(site, year, ipa, sep = "_"), .after = ipa) %>% 
#   select(!1:3) %>% 
#   column_to_rownames(var = "sample")
# 
# check_matrix <- fish_L.both %>% 
#   decostand(method = "pa") %>% 
#   mutate(no_spp = rowSums(.)) %>% 
#   filter(no_spp < 3)
# 
# fish_L.both <- fish_L.both %>% 
#   filter(!row.names(.) %in% row.names(check_matrix))

#### Create the Q (trait) matrix ####

#extract the names of the species present in the lampara net dataset
spp_names <- net_tidy %>% 
  distinct(ComName) %>% 
  mutate(Species = NA) 

#link these common names to their scientific names in fishbase
sci_names <- vector(mode = 'list', length = length(spp_names))

for (i in 1:nrow(spp_names)) {
  sci_names[[i]] <- rfishbase::common_to_sci(spp_names[i,])
  spp_names[i,2] <- ifelse(nrow(sci_names[[i]]) == 1, sci_names[[i]][[1]], NA)
}

#specify which entry to use for species that have multiple common names
spp_names <- spp_names %>% 
  mutate(Species = case_when(ComName == "Pacific herring" ~ "Clupea pallasii", 
                             ComName == "Pacific Cod" ~ "Gadus macrocephalus",
                             ComName == "Northern Anchovy" ~ "Engraulis mordax",
                             ComName == "Striped Seaperch" ~ "Embiotoca lateralis",
                             ComName == "Rock Sole" ~ "Lepidopsetta bilineata",
                             ComName == "Steelhead trout" ~ "Oncorhynchus mykiss",
                             ComName == "Cutthroat trout" ~ "Oncorhynchus clarkii",
                             ComName == "Tube-snout" ~ "Aulorhynchus flavidus",
                             TRUE ~ Species)) %>% 
  arrange(ComName) %>% 
  filter(!str_detect(ComName, 'UnID'))

#calculate the mean fork length of each species from the subsamples taken in the field
fork_length <- net_tidy %>% 
  group_by(ComName) %>% 
  summarize(mean_length_mm = mean(mean_length_mm)) %>% 
  inner_join(spp_names) %>% 
  mutate(mean_length_mm = ifelse(ComName == "Tidepool Sculpin", 89.0, mean_length_mm)) #no length in our df so taking the max length from fishbase

#extract other trait data from rfishbase
milieu <- species(spp_names$Species) %>% 
  mutate(DemersPelag = ifelse(Genus == "Oncorhynchus", "epipelagic", DemersPelag)) %>% #juvenile salmon are generally near the surface, 
  #as reviewed in Quinn 2018 and observed by Munsch et al. 2017, and for steelhead cite Daly et al. 2014
  select(Species, BodyShapeI, DemersPelag, AnaCat) %>% 
  mutate(DemersPelag = ifelse(DemersPelag == "pelagic-neritic", "pelagic", DemersPelag)) %>% #simplify this category because there is only one pelagic and two pelagic-neritic species
  mutate(migrations = ifelse(is.na(AnaCat), "non-migratory", AnaCat)) %>% #presumed non migratory if no information is available
  select(!AnaCat)

feeding_guild1 <- fooditems(spp_names$Species) %>% 
  select(Species, FoodI, FoodII, PredatorStage) %>% 
  group_by(Species, FoodI) %>% 
  summarize(count = n()) %>% 
  group_by(Species) %>% 
  mutate(per =  100 *count/sum(count)) %>% 
  filter(per >= 60) %>% #classify main food source if more than 60% of diet is in one category
  inner_join(spp_names) %>% 
  arrange(ComName) %>% 
  mutate(feeding_guild = case_when(FoodI == "zoobenthos" ~ "Zoobenthivorous",
                                   FoodI == "zooplankton" ~ "Planktivorous", 
                                   FoodI == "nekton" ~ "Piscivorous",
                                   TRUE ~ FoodI)) %>% 
  ungroup() %>% 
  select(Species, feeding_guild) %>% 
  mutate(feeding_guild = ifelse(Species == "Microgadus proximus", "Zoobenthivorous",feeding_guild)) %>% #juvenile diet
  mutate(feeding_guild = ifelse(Species == "Ophiodon elongatus", "Planktivorous",feeding_guild)) #juvenile lingcod eat copepods and other small cruscaceans - Fishbase citing Pacific Fishes of Canada

#if food items described don't include a dominant category, classify as omnivorous
feeding_guild2 <- fooditems(spp_names$Species) %>% 
  select(Species, FoodI, FoodII, PredatorStage) %>% 
  group_by(Species, FoodI) %>% 
  summarize(count = n()) %>% 
  group_by(Species) %>% 
  mutate(per = 100 *count/sum(count)) %>% 
  mutate(category = ifelse(per >= 60, FoodI, NA)) %>% 
  mutate(onlyNA = all(is.na(category))) %>% 
  filter(onlyNA == TRUE) %>% 
  inner_join(spp_names) %>% 
  select(!c(category, onlyNA)) %>% 
  mutate(feeding_guild = "Omnivorous") %>% 
  ungroup() %>% 
  select(Species, feeding_guild) %>% 
  distinct() %>% 
  mutate(feeding_guild = ifelse(Species == "Engraulis mordax", "Planktivorous",feeding_guild)) #phyto and zoo plankton


feeding_guild <- rbind(feeding_guild1, feeding_guild2) %>% 
  add_row(Species = "Blepsias cirrhosus", feeding_guild = "Zoobenthivorous") %>% #fishbase "Diet"
  add_row(Species = "Liparis florae", feeding_guild = "Zoobenthivorous")  %>% #don't have a great source for this one
  mutate(feeding_guild = str_to_lower(feeding_guild))

#chinook are omnivorous - duffy et al. 2011
#chinook and chum in eelgrass mostly eat epifaunal invertebrates - Kennedy et al. 2018
#coho and chinook eat some fish, otherwise they (and other species) eat a combination of zoobenthos and zooplankton - Beamish et al. 2003 (A history of reserach)

fish_traits <- full_join(fork_length, milieu) %>% 
  select(3,1,2,4,5,6) %>% 
  left_join(feeding_guild) %>% 
  mutate(ComName = replace(ComName, ComName == "Pacific Sandfish", "Pacific sandfish")) %>%
  arrange(ComName)

# write_csv(fish_traits, "data/fish_traits.csv")

fish_Q <- fish_traits %>% 
  select(-c(Species, ComName)) %>% 
  mutate_if(is.character, as.factor) %>% 
  clean_names() %>% 
  as.data.frame()

rownames(fish_Q) <- colnames(fish_L.year)

confirm_proper_names <- fish_Q %>% 
  mutate(Species = fish_traits$ComName, .before = mean_length_mm)

### Conduct data quality and integrity check using steps outlined in Palacio et al. 2022 ####
##Step 1: Plot the community data matrix to assess the prevalence of zeros 

fish_L.year %>% 
  rownames_to_column(var = "sample") %>% 
  pivot_longer(cols = 2:43, names_to = "species") %>% 
  ggplot() +
  geom_tile(aes(x = species, y = sample, fill = value))
#lots of zeros!

## Step 2: Check species sampling coverage (e.g. rarefaction)
#did this in the original data exploration phase

## Step 3: Plot the distribution of continuous traits to check for outliers and continuous traits to check the balance of levels 
par(mfrow=c(1,2))
hist(fish_Q$mean_length_mm, main="Histogram", xlab="Mean Length (mm)")
boxplot(fish_Q$mean_length_mm, main="Boxplot", ylab="Mean Length (mm)")
#trait is generally well distributed with few outliers

# #check the coefficient of variation
# CV <- function(x) { 100 * sd(x) / mean(x) }
# CV(fish_Q$mean_length_mm)
# # <50 is small, so we don't necessarily need to transform the length data
#log transform the length data
fish_Q.t <- fish_Q %>% mutate_if(is.numeric, log)
# CV(fish_Q.t$mean_length_mm)

par(mfrow=c(1,1))

barplot(table(fish_Q$body_shape_i, useNA="ifany"))
#body shape is well distributed with no missing data

barplot(table(fish_Q$demers_pelag, useNA="ifany"))
table(fish_Q$demers_pelag, useNA="ifany")
#combined pelagic and neritic groups above to better distribute

barplot(table(fish_Q$migrations, useNA="ifany"))
table(fish_Q$migrations, useNA="ifany")
#some categories could be grouped

barplot(table(fish_Q$feeding_guild, useNA="ifany"))
table(fish_Q$feeding_guild, useNA="ifany")
#distributed well enough

## Step 4: Evaluate multicollinearity among continuous traits and associations with categorical traits 
ggpairs(fish_Q.t)
#no evidence of multicollinearity

## Step 5: Identify missing trait data 
# missing data were added above

#### save final L and Q matrices #### 
fish.list <- list("trait" = fish_Q, 
                  "trait.t" = fish_Q.t,
                  # "abund.ipa" = fish_L.ipa,
                  "abund" = fish_L.year
                  # "abund.both" = fish_L.both
                  ) 

save(fish.list, file = here("data", "fish.list.Rdata"))
