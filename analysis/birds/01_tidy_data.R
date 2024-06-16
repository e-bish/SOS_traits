library(tidyverse)
library(here)
library(googlesheets4)
library(googledrive)
library(readxl)
library(auk) #for bird species names
library(moments)
library(janitor)
library(vegan)

OS_core_sites <- factor(c("FAM", "TUR", "COR", "SHR", "DOK", "EDG"), 
                         levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG"))

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

##### look at representative species for unknowns ####

#load trait data from Tobias et al. 2021 (Ecology Letters) AVONET: morphological, ecological and geographical data for all birds.
# AVONET <- here("data", "AVONET Supplementary dataset 1.xlsx") %>% 
#   read_excel(sheet = "AVONET2_eBird")
#
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

#expand the data to get sampling events with zeros
sampling_events <- predator_data %>% 
  select(year, month, day, site, ipa) %>% 
  distinct()
  
#format field data and add species ids
birds <- predator_data %>% 
  filter(species_group == "bird") %>% 
  filter(!(comm_name == "unk_bird" | comm_name == "unk_terrestrial_bird" | comm_name == "unk_duck")) %>% 
  full_join(sampling_events) %>% 
  arrange(year, day, month, site, ipa) %>% 
  mutate(abundance =replace_na(abundance, 0)) %>% 
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
                           'sparrow' = 'HOSP')) %>%  #representative sparrow species
  group_by(year, month, day, site, ipa, comm_name, spp_code) %>%
  mutate(year = as.factor(year)) %>% 
  summarize(spp_sum = sum(abundance)) %>% # there were multiple entries per species before because they were coded with behaviors
  ungroup()

##what were the most common species we caught?
birds %>%
  group_by(comm_name) %>%
  summarize(total = sum(spp_sum)) %>% 
  arrange(-total)

##how often did we encounter each species?
times_encountered <- birds %>% 
  filter(!is.na(comm_name)) %>% 
  group_by(comm_name, site) %>% 
  summarize(freq_obs = n())

#by tax group
ggplot(times_encountered, aes(x = comm_name, y = freq_obs, fill = factor(site, levels = SOS_sites))) + 
  geom_bar(position = "stack", stat="identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "", y = "Frequency observed", fill = "Site") 

ggplot(times_encountered, aes(x = factor(site, levels = SOS_sites), y = freq_obs, fill = comm_name)) + 
  geom_bar(position = "stack", stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "", y = "Frequency observed", fill = "tax group") 

##in how many sites does each species occur?
sites_encountered <- times_encountered %>% 
  mutate(pa = ifelse(freq_obs > 0, 1, 0)) %>% 
  group_by(comm_name) %>% 
  summarize(encounters_site = sum(pa)) %>% 
  arrange(desc(encounters_site))

ggplot(sites_encountered, aes(x = reorder(comm_name, -encounters_site), y = encounters_site)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs(x = "", y = "No. of sites encountered")

spp_by_site <- birds %>% 
  filter(!is.na(comm_name)) %>% 
  count(comm_name, site)

ggplot(spp_by_site, aes(x = reorder(comm_name, -n), y = n, fill = site)) +
  geom_col() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#core sites only
ggplot(filter(spp_by_site, site %in% SOS_core_sites), aes(x = reorder(comm_name, -n), y = n, fill = site)) +
  geom_col() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

spp_site_count <- spp_by_site %>% 
  # filter(site %in% SOS_core_sites) %>% 
  group_by(comm_name) %>% 
  summarize(n_sites = n()) %>% 
  arrange(desc(n_sites))

##what is the abundance of each species when it occurs?
spp_by_site %>% 
  # filter(site %in% SOS_core_sites) %>% 
  ggplot(aes(x = reorder(comm_name, -n), y = n)) +
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

##is the mean abundance correlated with the number of sites where it occurs?
mean_count <- spp_by_site %>% 
  group_by(comm_name) %>% 
  summarize(mean_n = mean(n))

cor(mean_count$mean_n, spp_site_count$n_sites) #no

##is the total abundance of fish correlated with the number of species in each site?
n_spp_site <- spp_by_site %>% 
  group_by(site) %>% 
  summarize(n_spp = n())

abund_by_site <- birds %>%
  filter(!is.na(comm_name)) %>% 
  group_by(site) %>%
  summarize(total = sum(spp_sum)) 

cor(n_spp_site$n_spp, abund_by_site$total) #not really

##did we sample enough to adequately describe the community?
make_filtered_spp_mat <- function(site_ID) {
  
  spp_mat <- birds %>% 
    filter(!is.na(comm_name)) %>% 
    mutate(year_month = ym(paste(year, month, sep = "_"))) %>% 
    filter(site == site_ID) %>% 
    group_by(comm_name, year_month) %>% 
    summarize(spp_sum = sum(spp_sum)) %>% 
    pivot_wider(names_from = comm_name, values_from = spp_sum) %>% 
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

curve_list <- lapply(SOS_core_sites, make_spp_curve)
curve_df <- data.frame() 

for (i in 1:length(SOS_core_sites)) {
  sites <- curve_list[[i]]$sites
  richness <- curve_list[[i]]$richness
  tmp <- data.frame(site = SOS_core_sites[i], samples = sites, richness = richness)
  
  curve_df <- rbind(curve_df, tmp)
}

ggplot() +
  geom_point(data = curve_df, aes(x = samples, y = richness, color = factor(site, level = SOS_core_sites))) +
  geom_line(data = curve_df, aes(x = samples, y = richness, color = factor(site, level = SOS_core_sites))) +
  theme_classic() +
  labs(x = "Times sampled", y = "Number of species", color = "Site")

ggsave("docs/figures/bird_specaccum_plot.png")

##is there obvious seasonality in our catch?
#in total catch abundance?
birds %>% 
  filter(!is.na(comm_name)) %>% 
  filter(site %in% SOS_core_sites) %>% #remove jubilee sites to properly compare June
  ggplot(aes(x = month, y = spp_sum, fill = year)) +
  geom_bar(stat = "identity") +
  labs(x = "Month", y = "Abundance", fill = "Year") + 
  theme_classic()
#no clear seasonal trend

#in species richness?
n_spp_by_month <- birds %>% 
  filter(!is.na(comm_name)) %>% 
  filter(site %in% SOS_core_sites) %>% #remove jubilee sites to properly compare June
  group_by(year, month) %>% 
  summarize(n_spp = n_distinct(comm_name)) 

n_spp_by_month %>% 
  ggplot(aes(x = month, y = n_spp, group = year, color = year)) +
  geom_line() +
  geom_point() + 
  theme_classic() + 
  labs(x = "Month", y = "Species Richness", color = "Year")
#again, doesn't seem to be a clear break between "early" and "late" 

##is there a difference between northern and southern sites?
abund_region <- birds %>% 
  filter(!is.na(comm_name)) %>% 
  mutate(year_month = ym(paste(year, month, sep = "-"))) %>% 
  mutate(region = case_when(site %in% c("FAM", "TUR", "COR", "MA", "HO", "WA") ~ "North",
                            TRUE ~ "South")) %>% 
  group_by(region, site, year_month) %>% 
  summarize(n_spp = n_distinct(comm_name))

abund_region %>% 
  ggplot(aes(x = site, y = n_spp, fill = region)) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "Site", y = "Species Richness", fill = "Region")

#### final tidy dataframe for analysis ####
birds_tidy <- birds %>% 
  filter(site %in% SOS_core_sites) #use only core sites because we didn't sample jubilee sites enough to capture the community

save(birds_tidy, file = here("data", "birds_tidy.Rdata"))
