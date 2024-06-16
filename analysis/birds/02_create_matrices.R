library(tidyverse)
library(auk) #for bird species names
library(moments)

#load tidy fish data frame created in 01_tidy_data
load(here("data", "birds_tidy.Rdata")) 

#load trait data from Tobias et al. 2021 (Ecology Letters) AVONET: morphological, ecological and geographical data for all birds.
AVONET <- here("data", "AVONET Supplementary dataset 1.xlsx") %>% 
  read_excel(sheet = "AVONET2_eBird")

expand_species <- birds_tidy %>% 
  expand(nesting(year, month, site, ipa), spp_code) %>% 
  filter(!is.na(spp_code))

#create abundance matrix
bird_L <- birds_tidy %>% 
  filter(!is.na(comm_name)) %>% 
  full_join(expand_species) %>% 
  mutate(spp_sum = replace_na(spp_sum, 0)) %>% 
  group_by(site, ipa, spp_code) %>% 
  summarize(spp_avg = mean(spp_sum)) %>% #average across months/years
  ungroup() %>% 
  arrange(spp_code) %>% 
  pivot_wider(names_from = spp_code, values_from = spp_avg, values_fill = 0) %>% 
  mutate(sample = paste(site, ipa, sep = "_"), .after = ipa) %>% 
  select(!1:2) %>% 
  column_to_rownames(var = "sample")

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
         Mass,
         Primary.Lifestyle, 
         Trophic.Level, 
         Migration,
         Trophic.Niche) %>% 
  mutate(Trophic.Niche = str_replace(Trophic.Niche, " ", "_")) %>% 
  mutate_if(is.character, as.factor) %>% 
  right_join(spp_id) %>% 
  select(-c(Species2, comm_name, spp_no)) %>% 
  select(spp_code, everything()) %>% 
  arrange(spp_code) %>% 
  column_to_rownames(var="spp_code")

# write.csv(bird_traits, "data/bird_traits.csv", row.names = TRUE)

#inspect continuous traits 
skewness(bird_traits$Mass)
range(bird_traits$Mass)

#apply data transformation to continuous traits

bird_traits.t <- bird_traits %>% 
  mutate_if(is.numeric, log) %>% 
  as.data.frame()

skewness(bird_traits.t$Mass)
range(bird_traits.t$Mass)

bird.list <- list("trait" = bird_traits.t, "abund" = bird_L) 
save(bird.list, file = here("data", "bird.list.Rdata"))
