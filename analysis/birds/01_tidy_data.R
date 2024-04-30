library(tidyverse)
library(here)
library(googlesheets4)
library(googledrive)
library(readxl)
library(auk) #for bird species names

#load field data
#first get data from google

## AUTHENTICATE ACCOUNT -- select the box that enables access to google (may have to do this 2x to make it actually work, 
# also if selecting an already authenticated account doesn't work then select 0 for a "new" account)
#drive_auth()

###############################################################################
# #use this to download raw data from the shared project drive
# predator_import <- drive_get("https://docs.google.com/spreadsheets/d/1Z_WaqjGCvTRWk1azulz5lBsPZXWBEb9IWa6Z1RNM3RM/edit#gid=0") %>%
#   read_sheet() 
# 
# #write to csv
# write_csv(predator_import, here("data", "raw", "predator_import.csv")) 

predator_data <- read_csv("data/raw/predator_import.csv", col_names = TRUE) %>% 
  rename(comm_name = "species") %>% 
  mutate_if(is.numeric, ~replace_na(., 0)) 

#load trait data from Tobias et al. 2021 (Ecology Letters) AVONET: morphological, ecological and geographical data for all birds.
AVONET <- here("data", "AVONET Supplementary dataset 1.xlsx") %>% 
  read_excel(sheet = "AVONET2_eBird")

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
                           'glaucous_gull' = "GWGU", #could be part of the western complex
                           'gull' = 'HERG', #representative gull species: Laurus argentatus
                           'swallow' = 'BARS', #representative swallow species: Hirundo rustica
                           'cormorant' = 'PECO', #representative cormorant species: Pelagic Cormorants 
                           'hummingbird' = 'RUHU', #representative hummingbird species: Selasphorus rufus
                           'grebe' = 'WEGR', #representative grebe species: Western Grebes
                           'sparrow' = 'HOSP')) %>% #representative sparrow species
  select(-c(weather, beaufort_sea_state, observer, start_time, 
            time_into_survey, species_group, behaviour, notes)) 

#create abundance matrix
bird_L_df <- birds %>% 
  filter(site %in% c("FAM", "TUR", "COR", "SHR", "DOK", "EDG")) %>%  #core sites if we're not using jubilee
  mutate(month = recode(month, `04` = "Apr", `05` = "May", `06` = "Jun", `07` = "Jul", `08` = "Aug", `09` = "Sept")) %>% 
  mutate(year = factor(year), month = factor(month, levels = c("Apr", "May", "Jun", "Jul", "Aug", "Sept"))) %>% 
  select(year, month, site, ipa, spp_code, abundance) %>% 
  complete(nesting(year, month, site), ipa, spp_code, fill = list(abundance = 0)) %>% #fill in zeros for shorelines we surveyed
  group_by(year, month, site, spp_code) %>%
  summarize(spp_mean = round(mean(abundance), 2)) %>% 
  ungroup() %>% 
  arrange(spp_code) %>% 
  pivot_wider(names_from = spp_code, values_from = spp_mean, values_fill = 0) 

bird_L_mat <- bird_L_df %>% 
  mutate(sample = paste(year, month, site, sep = "_")) %>% 
  column_to_rownames(var="sample") %>% 
  select(!c(year, month, site))

bird_L_mat.t <- round(sqrt(bird_L_mat), 2)

CV <- function(x) { 100 * sd(x) / mean(x) } #function from Jon Bakker

CV(x = rowSums(bird_L_mat)) #moderate CV McCune & Grace (2002, p.70)

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
                              TRUE ~ Species2)) %>% 
  mutate(spp_no = seq(1:nrow(.)))

### create the trait matrix
bird_traits <- AVONET %>% 
  select(Species2, 
         Trophic.Level, 
         Trophic.Niche, 
         Migration, 
         Primary.Lifestyle, 
         Mass, 
         Secondary1, 
         Wing.Length) %>% 
  mutate(Trophic.Niche = str_replace(Trophic.Niche, " ", "_")) %>% 
  mutate_if(is.character, as.factor) %>% 
  right_join(spp_id) %>% 
  select(-c(Species2, comm_name, spp_no)) %>% 
  select(spp_code, everything()) %>% 
  arrange(spp_code) %>% 
  column_to_rownames(var="spp_code")

# write.csv(traits, "data/bird_traits.csv", row.names = TRUE)

library(moments)
# library(janitor)
# library(StatMatch)

#apply data transformation to continuous traits

skewness(bird_traits$Mass)
range(bird_traits$Mass)

skewness(bird_traits$Secondary1)
range(bird_traits$Secondary1)

skewness(bird_traits$Wing.Length)
range(bird_traits$Wing.Length)

bird_traits.t <- bird_traits %>% 
  mutate_if(is.numeric, log) %>% 
  as.data.frame()

# skewness(bird_traits.t$Mass)
# range(bird_traits.t$Mass)
# 
# skewness(bird_traits.t$Secondary1)
# range(bird_traits.t$Secondary1)
# 
# skewness(bird_traits.t$Wing.Length)
# range(bird_traits.t$Wing.Length)

bird.list <- list("trait" = bird_traits.t, "abund" = bird_L_mat) #untransformed abundance matrix
save(bird.list, file = here("data", "bird.list.Rdata"))
