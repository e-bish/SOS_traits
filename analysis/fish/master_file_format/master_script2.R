# Load libraries
# tidy data
library(tidyverse)
library(here)
library(janitor)
library(vegan)

#create matrices
library(tidyverse)
library(here)
library(rfishbase)
library(janitor)
library(GGally)

#FD analysis
library(tidyverse)
library(here)
library(ggrepel)
library(FD)
library(ggordiplots)
library(PNWColors)
library(patchwork)
library(vegan)

#null model
library(here)
library(tidyverse)
library(picante)
library(FD)

#bootstrapping
library(here)
library(tidyverse)
library(janitor)
library(rsample)
library(mFD)
library(FD)
library(ggridges)
library(patchwork)

set.seed(1993)

#### Import Data in 01_import_data####
#### Tidy Data

SOS_sites <- factor(c("FAM", "TUR", "COR", "MA", "WA", "HO", "SHR", "DOK", "LL", "TL", "PR", "EDG"), 
                    levels = c("FAM", "TUR", "COR", "MA", "WA", "HO", "SHR", "DOK", "LL", "TL", "PR", "EDG"))
SOS_core_sites <- factor(c("FAM", "TUR", "COR", "SHR", "DOK", "EDG"), 
                         levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG"))
SOS_jubilee_sites <- factor(c("MA", "WA", "HO", "LL", "TL", "PR"),
                            levels = c("MA", "WA", "HO", "LL", "TL", "PR"))

load_data <- function() {
  #Load data
  #load field data (raw data downloaded/formatted in "01_import_data.R")
  net_2018.19 <- here::here("data","net_18.19.csv") %>% 
    read_csv()
  
  net_2021 <- here::here("data","net_2021.csv") %>% 
    read_csv()
  
  net_2022 <- here::here("data","net_2022.csv") %>% 
    read_csv()
  
  net_tidy <- bind_rows(net_2018.19, net_2021, net_2022) %>% 
    mutate(month = if_else(site == "MA", "06", month)) %>% # we did a July 1st survey at Maylor that we want to count as a June survey
    mutate(ipa = replace(ipa, site == "TUR" & ipa == "Restored", "Natural2")) %>% #no restoration at Turn Island
    mutate(date = make_date(year, month, day)) %>% 
    mutate(month = recode(month, `04` = "Apr", `05` = "May", `06` = "Jun", `07` = "Jul", `08` = "Aug", `09` = "Sept")) %>% 
    mutate(year = factor(year), 
           month = factor(month, levels = c("Apr", "May", "Jun", "Jul", "Aug", "Sept")),
           site = factor(site, levels = c("FAM", "TUR", "COR", "MA", "WA", "HO", "SHR", "DOK", "LL", "TL", "PR", "EDG"))) %>% 
    mutate(species_count = as.numeric(species_count)) %>%
    mutate(species_count = ifelse(is.na(org_type), 0, species_count)) 
  
  #check that we're keeping only the data we want
  sampling_events <- net_tidy %>% 
    select(year, month, day, site, ipa, station) %>% 
    distinct() %>% 
    filter(!(site == "COR" & year == "2019" & month == "Apr" & day == "30")) %>% #these COR seem like data entry errors because they were only sampled at one station and we already had complete sampling at COR in april and may. No fish recorded in either entry
    filter(!(site == "COR" & year == "2019" & month == "May" & day == "01")) %>% 
    filter(!(site == "TUR" & year == "2018" & month == "Jul" & day == "11")) %>%  #this is an incomplete sampling event. Complete sampling at TUR occured on 7/12/18
    filter(!(site == "TUR" & year == "2018" & month == "Sept" & day == "11")) #another month with repeat sampling; keeping only the second september TUR sampling event
  
  net_tidy2 <- net_tidy %>% 
    semi_join(sampling_events, join_by("year", "month", "day", "site")) %>%  #keep only the legit sampling events (days)
    filter(org_type == "Fish") %>% #### this is where you lose zero counts, if those are needed later
    full_join(sampling_events) %>% #add back in sampling events with zeros so they're represented in the species accumulation chart
    arrange(year, month, day, site, ipa, station) %>% 
    mutate(species_count=replace_na(species_count, 0)) %>% 
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
    filter(!grepl("UnID", ComName)) #remove unidentified species
  
  return(net_tidy2)
}

net_tidy <- load_data()
net_core <- net_tidy %>% 
  filter(site %in% SOS_core_sites) %>% 
  mutate(site = factor(site, levels = SOS_core_sites))

# final tidy data frame for analysis 
net_tidy <- net_tidy %>% 
  filter(site %in% SOS_core_sites) #use only core sites because we didn't sample jubilee sites enough to capture the community

# save(net_tidy, file = here("data", "net_tidy.Rdata"))

#### Create Matrices ####
#load tidy fish data 
# load(here("data", "net_tidy.Rdata")) 

## create the abundance matrix 
# extract unique sampling events to quantify sample effort (slightly unbalanced between depths and shorelines)
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

#expand species by sampling event
expand_species <- net_tidy %>%
  expand(nesting(year, month, site, ipa), ComName) %>%
  filter(!is.na(ComName))

#create a version of the abundance matrix that has years and ipas
fish_L <- net_tidy %>% #L is referring to the RLQ analysis
  filter(!is.na(ComName)) %>% 
  group_by(year, month, site, ipa, ComName) %>% 
  summarize(spp_sum = sum(species_count)) %>% #sum across sampling events 
  ungroup() %>%
  full_join(sampling_events) %>% 
  mutate(catch_per_set = spp_sum/no_net_sets) %>% 
  filter(!is.na(ComName)) %>% #these are accounted for in the next step
  full_join(expand_species) %>% #add back in all of the events so we can capture 0s
  mutate(ComName = replace(ComName, ComName == "Pacific Sandfish", "Pacific sandfish")) %>% #if it's capitolized it gets confused about the proper order
  arrange(year, month, site, ipa, ComName) %>% 
  mutate(catch_per_set  = replace_na(catch_per_set , 0)) %>% 
  group_by(year, site, ipa, ComName) %>% 
  summarize(avg_cps = mean(catch_per_set)) %>% #average across months within a year for each ipa
  pivot_wider(names_from = ComName, values_from = avg_cps, values_fill = 0) %>% 
  clean_names() %>% 
  ungroup() %>% 
  mutate(sample = paste(site, ipa, year, sep = "_"), .after = ipa) %>% 
  select(!1:3) %>% 
  column_to_rownames(var = "sample")

## create the trait matrix 

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

#log transform continuous variables
fish_Q.t <- fish_Q %>% mutate_if(is.numeric, log)

#save final matrices
fish.list <- list("trait" = fish_Q, 
                  "trait.t" = fish_Q.t,
                  "abund" = fish_L) 

# save(fish.list, file = here("data", "fish.list.Rdata"))

#### FD Analysis ####
# load(here("data", "fish.list.Rdata")) #object created in 03_create_matrices

#could try checking this with the mFD package, but the lit supports retaining 4
n_axes_to_keep <- 4

fishFD <- dbFD(x = fish.list$trait.t, 
               a = fish.list$abund,
               ord = "podani",
               corr = "cailliez", 
               m = n_axes_to_keep,
               calc.FDiv = TRUE, 
               print.pco = FALSE)

FD_values <- cbind(fishFD$nbsp, fishFD$FRic, fishFD$FEve, fishFD$FDiv, fishFD$FDis) #extract indices
colnames(FD_values) <- c("Species_Richness", "FRic", "FEve", "FDiv", "FDis")

FD_results <- FD_values %>% 
  as_tibble(rownames = "sample") %>% 
  replace(is.na(.),0) %>% #some rows don't have enough species to calculate indices
  separate_wider_delim(sample, delim = "_", names = c("site", "ipa", "year"), cols_remove = TRUE) %>% 
  mutate(site = factor(site, levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG")),
         year = factor(year, levels = c("2018", "2019", "2021", "2022")),
         region = ifelse(site %in% c("FAM", "TUR", "COR"), "North", "South"), 
         eelgrass = ifelse(site %in% c("TUR", "COR", "SHR"), "yes", "no"), .after = site,
         ipa = ifelse(ipa == "Natural2", "Natural", ipa)) #combine the two natural sites at TUR

# save(FD_results, file = "data/FD_results.Rdata")  

#plot the results
FD_results %>% 
  pivot_longer(!c(site, ipa, region, eelgrass, year), names_to = "metric", values_to = "value") %>% 
  ggplot(aes(x = site, y = value, fill = eelgrass)) +
  geom_boxplot() +
  geom_point(show.legend = FALSE) +
  theme_classic() +
  facet_wrap(~metric, scales = "free_y")

FD_results %>% 
  pivot_longer(!c(site, ipa, region, eelgrass, year), names_to = "metric", values_to = "value") %>% 
  ggplot(aes(x = ipa, y = value)) +
  geom_boxplot() +
  geom_point(show.legend = FALSE) +
  theme_classic() +
  facet_wrap(~metric, scales = "free_y")

#check to see if sites are different while controlling for year
adonis2(FD_results[,c("Species_Richness","FDis", "FEve", "FRic", "FDiv")] ~ site, 
        strata = FD_results$year, data = FD_results, method = "euc")
#yes

#check to see if years are the same
adonis2(FD_results[,c("Species_Richness","FDis", "FEve", "FRic", "FDiv")] ~ year, 
        strata = FD_results$site, data = FD_results, method = "euc")
#very slightly different

#### compute null model ####
#load tidy fish data frame created in 02_tidy_data
# load(here("data", "net_tidy.Rdata"))
# load(here("data", "fish.list.Rdata")) #object created in 03_create_matrices

null_L <- list(fish.list$abund)
n_iter <- 1000 #observed data + 999 permutations

# #frequency is C4 null model in Gotzenberger et al. that works well for detection of environmental filtering
# #independentswap seems to be the method that other FD papers have used (e.g., Zhang) and is recommended by Swenson

for (i in 2:n_iter) {
  set.seed(i)
  null_L[[i]] <- randomizeMatrix(fish.list$abund, null.model ="independentswap", iterations = 1000)
}

# # calculate FD indices with the FD package
FD_null_results <- list()
for (i in 1:n_iter){
  
  null_fishFD <- dbFD(x = fish.list$trait.t, #must be a df where character columns are factors or a distance matrix
                      a = null_L[[i]],
                      corr = "cailliez",
                      m = n_axes_to_keep,
                      calc.FDiv = TRUE,
                      print.pco = FALSE)
  
  null_FD_values <- cbind(null_fishFD$nbsp, null_fishFD$FRic, null_fishFD$FEve, null_fishFD$FDiv,
                          null_fishFD$FDis) #extract indices
  
  FD_null_df <- null_FD_values %>%
    as_tibble(rownames = "sample")
  
  FD_null_results <- rbind(FD_null_results, FD_null_df)
}

colnames(FD_null_results)[2:6] <- c("Species_Richness", "FRic", "FEve", "FDiv", "FDis")

FD_null_results <- FD_null_results %>%
  separate_wider_delim(sample, delim = "_", names = c("site", "year"), cols_remove =  TRUE) %>%
  mutate(site = factor(site, levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG")),
         year = factor(year, levels = c("2018", "2019", "2021", "2022"))) %>%
  mutate(region = ifelse(site %in% c("FAM", "TUR", "COR"), "North", "South"), .after = site)

# save(FD_null_results, file = "data/FD_null_results.Rda")

## load null dataset 
# load("data/FD_null_results.Rda")  

FD_null_summary <- FD_null_results %>% 
  group_by(site) %>% 
  summarize(across(where(is.numeric), list(mean = mean, sd = sd)))

#### calculate SES samples ####
# load("data/FD_results.Rda")  

FD_means <- FD_results %>% 
  replace(is.na(.), 0) %>% #replace na with zero for EDG armored where there was only two species caught at one ipa so cant calculate FRic, FEve, FDiv
  group_by(site) %>% 
  summarize(across(where(is.numeric), mean))

SES_tab <- as.data.frame(FD_means$site)
SES_tab[,2] <- FD_means$Species_Richness 
SES_tab[,3] <- (FD_means$FRic - FD_null_summary$FRic_mean) / FD_null_summary$FRic_sd
SES_tab[,4] <- (FD_means$FEve - FD_null_summary$FEve_mean) / FD_null_summary$FEve_sd
SES_tab[,5] <- (FD_means$FDiv - FD_null_summary$FDiv_mean) / FD_null_summary$FDiv_sd
SES_tab[,6] <- (FD_means$FDis - FD_null_summary$FDis_mean) / FD_null_summary$FDis_sd

names(SES_tab) <- names(FD_means)

# calculate p values for SES
#cant get this to work in a function!!! had to create a different function for each metric
pull_FRic_ntiles <- function(site_ID) {
  #commented out code left here for checking that the proper order has been extracted
  lower <- FD_null_results %>% 
    filter(site == site_ID) %>% 
    # rownames_to_column(var = "index") %>% 
    arrange(FRic) %>% 
    # rownames_to_column(var = "arranged") %>% 
    slice(25) %>%
    select(FRic)
  # select(index, arranged, FRic)
  
  upper <- FD_null_results %>%
    filter(site == site_ID) %>%
    # rownames_to_column(var = "index") %>% 
    arrange(FRic) %>% 
    # rownames_to_column(var = "arranged") %>% 
    ungroup() %>% 
    slice(975) %>%
    select(FRic)
  # select(index, arranged, FRic)
  
  names <- c("site", 
             # "index_lower", 
             # "arranged_lower",
             "FRic_lower",
             # "index_upper", 
             # "arranged_upper", 
             "FRic_upper")
  
  df <- data.frame(site_ID)
  df <- cbind(df, lower, upper)
  names(df) <- names
  
  return(df)
}
pull_FEve_ntiles <- function(site_ID) {
  
  lower <- FD_null_results %>% 
    filter(site == site_ID) %>% 
    arrange(FEve) %>% 
    slice(25) %>%
    select(FEve)
  
  upper <- FD_null_results %>%
    filter(site == site_ID) %>%
    arrange(FEve) %>% 
    ungroup() %>% 
    slice(975) %>%
    select(FEve)
  
  names <- c("site", "FEve_lower","FEve_upper")
  
  df <- data.frame(site_ID)
  df <- cbind(df, lower, upper)
  names(df) <- names
  
  return(df)
}
pull_FDiv_ntiles <- function(site_ID) {
  
  lower <- FD_null_results %>% 
    filter(site == site_ID) %>% 
    arrange(FDiv) %>% 
    slice(25) %>%
    select( FDiv)
  
  upper <- FD_null_results %>%
    filter(site == site_ID) %>%
    arrange(FDiv) %>% 
    ungroup() %>% 
    slice(975) %>%
    select(FDiv)
  
  names <- c("site", "FDiv_lower","FDiv_upper")
  
  df <- data.frame(site_ID)
  df <- cbind(df, lower, upper)
  names(df) <- names
  
  return(df)
}
pull_FDis_ntiles <- function(site_ID) {
  
  lower <- FD_null_results %>% 
    filter(site == site_ID) %>% 
    arrange(FDis) %>% 
    slice(25) %>%
    select(FDis)
  
  upper <- FD_null_results %>%
    filter(site == site_ID) %>%
    arrange(FDis) %>% 
    ungroup() %>% 
    slice(975) %>%
    select(FDis)
  
  names <- c("site",  "FDis_lower","FDis_upper")
  
  df <- data.frame(site_ID)
  df <- cbind(df, lower, upper)
  names(df) <- names
  
  return(df)
}

#extract the 95% confidence intervals for each metric
FRic_CI <- lapply(SOS_core_sites, pull_FRic_ntiles)
FEve_CI <- lapply(SOS_core_sites, pull_FEve_ntiles)
FDiv_CI <- lapply(SOS_core_sites, pull_FDiv_ntiles)
FDis_CI <- lapply(SOS_core_sites, pull_FDis_ntiles)

#combine the lists using mapply
combine_CI <- t(mapply(c, FRic_CI, FEve_CI, FDiv_CI, FDis_CI)) %>% as.data.frame()

#turn the confidence intervals into long format
CI_uppers_lowers <- combine_CI %>% 
  select(!starts_with("site")) %>% 
  mutate(site = SOS_core_sites) %>% 
  pivot_longer(!site, names_to = "metric") %>% 
  separate_wider_delim(metric, delim = "_", names = c("metric", "range"), cols_remove = TRUE) %>% 
  mutate(value = unlist(value))

#turn the mean values into long format
FD_means_long <- FD_means %>% 
  pivot_longer(!c(site, Species_Richness), names_to = "metric", values_to = "value") %>% 
  select(!Species_Richness) %>% 
  mutate(range = "mean", .before = value)

#combine the intervals with the means into one df and check significance
df_for_p_vals <- bind_rows(CI_uppers_lowers,FD_means_long) %>% 
  pivot_wider(names_from = range, values_from = value) %>% 
  mutate(significant = ifelse(mean > lower & mean < upper, "no", "yes")) 

#format significance checks into a table that's the same format as the SES table
p_vals_tbl <- df_for_p_vals %>% 
  select(site, metric, significant) %>% 
  pivot_wider(names_from = metric, values_from = significant)

#### Bootstrapping ####

#create abundance matrix
mat_to_boot <- net_tidy %>% 
  filter(!is.na(ComName)) %>% 
  group_by(year, month, site, ComName) %>% 
  summarize(spp_sum = sum(species_count)) %>% #sum across samples on each day
  ungroup() %>% 
  full_join(expand_species_mo.yr) %>% 
  mutate(site = factor(site, levels = SOS_core_sites)) %>% 
  mutate(spp_sum = replace_na(spp_sum, 0)) %>% 
  mutate(ComName = replace(ComName, ComName == "Pacific Sandfish", "Pacific sandfish")) %>% #if it's capitolized it gets confused about the proper order
  arrange(ComName) 

#bootstrap & keep equal proportions for each month to maintain seasonal variability
#using the strata argument here instead of the group option, because the group option results in 
#some assessment datasets with zero rows (meaning that the analysis dataset is the exact same as the raw data)
#the strata argument works well when categorical variables are unbalanced but you want similar proportions
boot_obj <- bootstraps(mat_to_boot, strata = month, times = 999)

#store the bootstrap dataframes in a list
boot_list <- list()
for (i in seq_along(boot_obj$splits)) {
  resample <- boot_obj$splits[[i]]
  boot_list[[i]] <- analysis(resample)
  
}

#format the bootstrapped data frame into properly formatted matrices
format_fish_L <- function(df) {
  
  fish_L <- df %>% 
    full_join(sampling_events_mo.yr) %>% 
    mutate(catch_per_set = spp_sum/total_net_sets) %>% 
    group_by(year, site, ComName) %>%  
    summarize(avg_cps = mean(catch_per_set)) %>% #average across months within a year
    ungroup() %>% 
    arrange(year, site, ComName) %>% 
    pivot_wider(names_from = ComName, values_from = avg_cps, values_fill = 0, names_sort = TRUE) %>% #you have to use names sort or else species missing in the first matrix will move to the end 
    clean_names() %>% 
    mutate(sample = paste(site, year, sep = "_"), .after = site) %>% 
    select(!1:2) %>% 
    column_to_rownames(var = "sample") 
  
  return(fish_L)
}

boot_L <- lapply(boot_list, format_fish_L)
# save(boot_L, file = here("data", "boot_L.Rda"))

#remove columns for species that aren't represented in some assemblages
remove_missing_spp <- function(df) {
  filtered_df <- df %>% 
    select_if(colSums(.) != 0)
  
  return(filtered_df)
}

boot_L_filtered <- lapply(boot_L, remove_missing_spp)

#create a list of trait matrices based on species represented in the bootstrapped L matrices
fish_Q_list <- list() 

for (i in 1:length(boot_L)) {
  
  fish_Q_list[[i]] <- fish.list$trait.t %>%
    rownames_to_column(var = "species") %>%
    filter(species %in% colnames(boot_L_filtered[[i]])) %>%
    column_to_rownames(var = "species")
}

# calculate alpha diversity with the FD package
gowdist.list <- lapply(fish_Q_list, gowdis, ord = "podani")
FD_boot_output <- list()

for (i in 1:length(boot_L)){
  
  fishFD.boot <- dbFD(x = gowdist.list[[i]], #must be a distance object or df where character columns are factors
                 a = data.matrix(boot_L_filtered[[i]]),
                 corr = "cailliez", 
                 m = n_axes_to_keep,
                 calc.FDiv = TRUE,
                 print.pco = TRUE)
  
  FD_values.boot <- cbind(fishFD.boot$nbsp, fishFD.boot$FRic, fishFD.boot$FEve, fishFD.boot$FDiv, fishFD.boot$FDis) #extract indices
  
  FD_boot_results_df <- FD_values.boot %>%
    as_tibble(rownames = "site")
  
  FD_boot_output <- rbind(FD_boot_output, FD_boot_results_df)
}

colnames(FD_boot_output)[2:6] <- c("Species_Richness", "FRic", "FEve", "FDiv", "FDis")

FD_boot_results <- FD_boot_output %>%
  separate_wider_delim(site, delim = "_", names = c("site", "year"), cols_remove =  TRUE) %>%
  replace(is.na(.), 0) %>% 
  mutate(site = factor(site, levels = SOS_core_sites), 
         region = ifelse(site %in% c("FAM", "TUR", "COR"), "North", "South"), 
         veg = ifelse (site %in% c("TUR", "COR", "SHR"), "present", "absent"), .after = year) 

# save(FD_boot_results, file = "data/FD_boot_results.Rda")

index_names <- c("Species Richness", "F. Richness", "F. Eveness", "F. Divergence", "F. Dispersion")

plot_site_index <- function (index){
  ggplot(data = FD_boot_results, aes(x = .data[[index]], 
                                y = factor(site, levels = rev(SOS_core_sites)), 
                                fill = region, color = region)) +
    geom_density_ridges(alpha = 0.9) + 
    theme_classic() +
    theme(axis.title.y = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
}

index_plots <- lapply(names(FD_boot_results[6:9]), plot_site_index)

for (i in 1:length(index_plots)) {
  index_plots[[i]] <- index_plots[[i]] + xlab(index_names[i+1])
}

#FEve and FDiv are scaled between 0,1
index_plots$"F. Evenness" <- index_plots$"F. Evenness" + xlim(0,1)
index_plots$"F. Divergence" <- index_plots$"F. Divergence" + xlim(0,1)

index_plots[[5]] <- ggplot(data = FD_boot_results, 
                           aes(x = Species_Richness, 
                               y = factor(site, levels = rev(SOS_core_sites)), 
                               color = region)) +
  geom_density_ridges2(aes(fill = region), 
                       stat = "binline", 
                       binwidth = 1, 
                       scale = 0.95, 
                       alpha = 0.9) +
  theme_classic() +
  xlab("Species Richness") +
  theme(axis.title.y = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

index_plots[[5]]  + index_plots[[1]]  + guide_area() + index_plots[[2]] + index_plots[[3]] + index_plots[[4]] + 
  plot_layout(ncol = 3, guides = "collect")

# ggsave("docs/figures/fish_FDbootpatch.png")

FD_boot_results %>% 
  filter(site == "COR") %>% 
  ggplot(aes(x = Species_Richness)) + 
  geom_histogram()
#why does this look weird?