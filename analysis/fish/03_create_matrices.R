library(tidyverse)
library(here)
library(rfishbase)
library(janitor)
library(GGally)

#load tidy fish data frame created in 02_tidy_data
load(here("data", "net_tidy.Rdata")) 

#### Create the L (abundance) matrix ####
expand_species <- net_tidy %>% 
  expand(nesting(year, month, site, ipa), ComName) %>% 
  filter(!is.na(ComName))

fish_L <- net_tidy %>% #L is referring to the RLQ analysis
  filter(!is.na(ComName)) %>% 
  group_by(year, month, site, ipa, ComName) %>%
  summarize(spp_sum = sum(species_count)) %>% #sum across depths within a site
  ungroup() %>%
  full_join(expand_species) %>% #add back in all of the events so we can capture 0s
  mutate(ComName = replace(ComName, ComName == "Pacific Sandfish", "Pacific sandfish")) %>% #if it's capitolized it gets confused about the proper order
  arrange(year, month, site, ipa) %>% 
  mutate(spp_sum =replace_na(spp_sum, 0)) %>% 
  group_by(site, ipa, ComName) %>% 
  summarize(spp_avg = mean(spp_sum)) %>% #average across sampling events at each site/ipa (unbalanced)
  pivot_wider(names_from = ComName, values_from = spp_avg, values_fill = 0) %>% 
  clean_names() %>% 
  ungroup() %>% 
  mutate(sample = paste(site, ipa, sep = "_"), .after = ipa) %>% 
  select(!1:2) %>% 
  column_to_rownames(var = "sample")

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

rownames(fish_Q) <- colnames(fish_L)

confirm_proper_names <- fish_Q %>% 
  mutate(Species = fish_traits$ComName, .before = mean_length_mm)

### Conduct data quality and integrity check using steps outlined in Palacio et al. 2022 ####
##Step 1: Plot the community data matrix to assess the prevalence of zeros 

fish_L %>% 
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
                  "abund" = fish_L) 

save(fish.list, file = here("data", "fish.list.Rdata"))
