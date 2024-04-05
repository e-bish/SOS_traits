library(tidyverse)
library(here)
# library(readxl)
# library(lubridate)
library(googlesheets4)
library(googledrive)

#load field data
#first get data from google

## AUTHENTICATE ACCOUNT -- select the box that enables access to google (may have to do this 20x to make it actually work, 
# also if selecting an already authenticated account doesn't work then select 0 for a "new" account)
#drive_auth()

###############################################################################
# #use this to download raw data from the shared project drive
# predator_import <- drive_get("https://docs.google.com/spreadsheets/d/1Z_WaqjGCvTRWk1azulz5lBsPZXWBEb9IWa6Z1RNM3RM/edit#gid=0") %>%
#   read_sheet( ) 
# 
# #write to csv
# write_csv(predator_import, here("data", "raw", "predator_import.csv")) 

predator_data <- read_csv("data/raw/predator_import.csv", col_names = TRUE) %>% 
  rename(comm_name = "species") %>% 
  mutate_if(is.numeric, ~replace_na(., 0)) 
#   
# summary <- predator_data %>% 
#   group_by(year, month, day, site, ipa) %>% 
#   summarize(n())

##### look at representative species for unknowns ####
# Larus <- AVONET %>% 
#   filter(Species2 == "Larus glaucescens" | Species2 == "Larus hyperboreus" | Species2 == "Larus argentatus" | Species2 == 
#            "Larus occidentalis")
# #select herring gulls (Laurus argentatus)
# 
# #compare purple martin, violet green, and barn swallow
# swallows <- AVONET %>% 
#   filter(Species2 == "Progne subis" | Species2 == "Tachycineta thalassina" | Species2 =="Hirundo rustica")
# #select barn swallow (Hirundo rustica) as representative species based on reasons outlined in the methods
# 
# #compare terns
# terns <- AVONET %>% 
#   filter(Species2 == "Hydroprogne caspia" | Species2 == "Sterna hirundo")
# #select caspian tern because I think they're most common
# 
# #compare hummingbirds
# humms <- AVONET %>% 
#   filter(Species2 == "Selasphorus rufus" | Species2 == "Calypte anna")
# #super similar so going with the Rufous hummingbird
# 
# #only kind of sparrow is house sparrow Passer domesticus
# #unknown grebes look like Western Grebes Aechmophorus occidentalis
# #selecting pelagic cormorant Urile pelagicus for all unidentified cormorants because only one was observed each time and this spp is more independent than others
# cor <- predator_data %>% 
#   filter(comm_name == "cormorant")

#load species ids for birds observed in the field
spp_id <- read_csv("data/bird_spp_info.csv", col_names = TRUE) %>% 
  unite(Species2, c("genus", "species"), sep = " ") %>% 
  filter(!(comm_name == "unk_bird" |
             comm_name == "unk_terrestrial_bird" |
             comm_name == "unk_duck" |
             comm_name == "Harbor Seal")) %>%
  mutate(spp_code = case_when(comm_name == "Common/Caspian" ~ "CATE", 
                              comm_name == "Sparrow" ~ 'HOSP',
                              comm_name == "grebe" ~ 'WEGR', 
                              comm_name == "Hummingbird" ~ 'RUHU',
                              comm_name == "Cormorant" ~ 'PECO',
                              comm_name == "swallow" ~ "BARS", 
                              comm_name == "Glaucus winged gull" ~ 'GWGU',
                              comm_name == "gull" ~ 'HERG',
                              TRUE ~ spp_code)) %>% 
  mutate(Species2 = case_when(comm_name == "Common/Caspian" ~ 'Hydroprogne caspia', 
                              comm_name == "Double-crested cormorant" ~ 'Nannopterum auritum', #updated species name in 2014
                              comm_name == "Sparrow" ~ 'Passer domesticus',
                              comm_name == "grebe" ~ 'Aechmophorus occidentalis', 
                              comm_name == "Hummingbird" ~ 'Selasphorus rufus',
                              comm_name == "Cormorant" ~ 'Urile pelagicus',
                              comm_name == "swallow" ~ 'Hirundo rustica', 
                              comm_name == "Glaucus winged gull" ~ 'Larus glaucescens',
                              comm_name == "gull" ~ 'Larus argentatus',
                              TRUE ~ Species2)) #%>% select(-c(comm_name))

#format field data and add species ids
birds <- predator_data %>% 
  filter(species_group == "bird") %>% 
  filter(!(comm_name == "unk_bird" | comm_name == "unk_terrestrial_bird" | comm_name == "unk_duck")) %>% 
  mutate(spp_code = recode(comm_name, 
                           "crow" = "COBR", 
                           "pigeon_guillemot" = 'PIGU',
                           'common_murre' = "COMU",
                           'american_robin' = "AMRO",
                           'bald_eagle' = "BAEA",
                           'pigeon' = "RODO",
                           'great_blue_heron' = "GBHE",
                           'osprey' = "OSPR",
                           'kingfisher' = 'BEKI',
                           'canada_geese' = 'CCGO',
                           'purple_martin' = 'PUMA',
                           'raven' = 'CORA',
                           'mallard' = 'MALL',
                           'turkey_vulture' = 'TUVU',
                           'bufflehead' = "BUFF",
                           'stellars_jay' = 'STJA',
                           'starling' = "EUST",
                           'horned_grebe' = 'HOGR',
                           'surf_scoter'='SUSC',
                           'double_crested_cormorant' = 'DCCO',
                           'tern' = "CATE", #representative tern species: Hydroprogne caspia
                           'glaucous_gull' = "GWGU", #don't necessarily trust this
                           'gull' = 'HERG', #representative gull species: Laurus argentatus
                           'swallow' = 'BARS', #representative swallow species: Hirundo rustica
                           'cormorant' = 'PECO', #representative cormorant species: Pelagic Cormorants 
                           'hummingbird' = 'RUHU', #representative hummingbird species: Selasphorus rufus
                           'grebe' = 'WEGR', #representative grebe species: Western Grebes
                           'sparrow' = 'HOSP')) %>% #representative sparrow species
  filter(year == 2022) %>% #same number of sampling events for each site_ipa this year
  mutate(date = make_date(year, month, day)) %>% 
  filter(!site == "TUR") %>% #remove Turn Island samples from this analysis until we figure out what to do with them
  filter(nchar(site) == 3) %>% #only core sites for now
  mutate(site_ipa = paste(site, ipa, sep = "_")) %>% 
  select(-c(year, month, day, ipa, site, weather, beaufort_sea_state, observer, start_time, 
            time_into_survey, comm_name, species_group, behaviour, notes))

#only keep species present in 2022
spp_id <- spp_id %>% 
  filter(spp_code %in% birds$spp_code) %>% 
  arrange(spp_code) 


spp_id <- spp_id %>% mutate(spp_no = seq(1:nrow(spp_id)))

#write.csv(traits, "project/traits_export.csv", row.names = TRUE)
