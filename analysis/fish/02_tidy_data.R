library(tidyverse)
library(here)
library(janitor)
library(rfishbase)

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
  mutate(site_ipa = paste(site, ipa, sep = "_")) %>% 
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

################################################################################
#prepare field abundance data for analysis
fish_N <- net_tidy %>% 
  select(date, site_ipa, ComName, species_count) %>% 
  arrange(ComName)

spp_names <- fish_N %>% 
  distinct(ComName) %>% 
  mutate(Species = NA)

sci_names <- vector(mode = 'list', length = length(spp_names))

for (i in 1:nrow(spp_names)) {
  sci_names[[i]] <- rfishbase::common_to_sci(spp_names[i,])
  spp_names[i,2] <- ifelse(nrow(sci_names[[i]]) == 1, sci_names[[i]][[1]], NA) 
}

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

MaxN <- fish_N %>% 
  group_by(site_ipa, ComName) %>%
  filter(species_count == max(species_count)) %>% 
  ungroup() %>% 
  rename(MaxN = "species_count") %>% 
  select(!date) %>% 
  distinct(site_ipa, ComName, MaxN)
 
fish_MaxN <- MaxN %>% 
  complete(site_ipa, ComName) %>% 
  replace(is.na(.), 0) %>% 
  pivot_wider(names_from = ComName, values_from = MaxN) %>% 
  select(-contains("UnID")) %>% 
  column_to_rownames(var="site_ipa") %>% 
  clean_names() %>% 
  as.matrix()

###############################################################################
# trait data

fork_length <- net_tidy %>% 
  filter(ComName == "Tidepool Sculpin")
  group_by(ComName) %>% 
  summarize(mean_fork_length = mean(mean_length_mm)) %>% 
  # mutate(fork_length = case_when(mean_fork_length < 70 ~ "small", ## does this need to be categorical??
  #                                mean_fork_length > 150  ~ "large",
  #                                TRUE ~ "medium")) %>% 
  filter(ComName %in% spp_names$ComName)

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
  select(ComName, feeding_guild)

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
  select(ComName, feeding_guild) %>% 
  distinct()

feeding_guild <- rbind(feeding_guild1, feeding_guild2) %>% 
  add_row(ComName = "Silverspotted sculpin", feeding_guild = "Zoobenthivorous") %>% #fishbase "Diet"
  add_row(ComName = "Tidepool Snailfish", feeding_guild = "Zoobenthivorous") %>% #don't have a great source for this one
  arrange(ComName)

body_transverse_shape <- morphology(spp_names$Species) %>% 
  select(Species, BodyShapeI) %>% 
  distinct() %>% 
  mutate(BodyShapeI = ifelse(Species == "Psychrolutes paradoxus", "elongated", BodyShapeI)) %>% 
  inner_join(spp_names) %>% 
  select(!Species) %>% 
  arrange(ComName)


milieu <- spp_names %>% 
  mutate(water_position = c("demersal"))

fish_traits <- inner_join(fork_length, body_transverse_shape)
fish_traits <- left_join(fish_traits, feeding_guild)





