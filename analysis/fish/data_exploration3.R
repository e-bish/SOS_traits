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
  mutate(year = factor(year), month = factor(month, levels = c("Apr", "May", "Jun", "Jul", "Aug", "Sept"))) %>% 
  filter(site %in% SOS_core_sites) #filter to core sites

#look at seasonality in catch by abundance
ggplot(net_tidy, aes(x = month, y = species_count, fill = year)) +
  geom_bar(stat = "identity")
#doesn't seem to be a clear break between "early" and "late" 

#look at seasonality in catch by spp richness
net_tidy %>% 
  group_by(year, month) #?? get # species per sampling event 


#look at differences by basin
net_tidy %>% 
  mutate(region = case_when(site %in% c("EDG", "DOK", "SHR") ~ "south",
                   TRUE ~ "north")) # ?? look at number of species by region?