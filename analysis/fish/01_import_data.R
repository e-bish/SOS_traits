### Code to import data from SOS shoreline restoration project phases 1 (2018/2019) and 2 (2021/2022)

# Load libraries
library(tidyverse)
library(here)
library(stringr)
library(googlesheets4)
library(googledrive)

#first get data from google

## AUTHENTICATE ACCOUNT -- select the box that enables access to google (may have to do this 20x to make it actually work, 
# also if selecting an already authenticated account doesn't work then select 0 for a "new" account)
#drive_auth()

###############################################################################
## 2018 & 2019 data
## adapted from net_base.R file (code by Genoa Sullaway)

#use this if you are using googlesheets
#gs_auth(new_user = TRUE) #authorize account
#net_import <- gs_title("Lampara Net Data")
#net_import<- gs_read(net_import)

#use this if you are using googlesheets4
net_import <- drive_get("https://docs.google.com/spreadsheets/d/1OhsndJNLAlHxT0TTHp7dmgbKgPZ7dgd-xNnE-aAxLdc/edit#gid=0") %>%
  read_sheet( ) 

# write_csv(net_import, here("data", "raw", "raw.18.19.csv"))

net_import <-net_import %>%
  separate(date, into = c("year","month", "day"), sep = "-") %>%
  mutate(month = str_pad(month, 2, side = c("left"), pad = "0")) %>%
  mutate(day = str_pad(day, 2, side = c("left"), pad = "0")) %>%
  mutate(tax_group = replace(tax_group, species == "Tube Snout", "Aulorhynchus")) %>% 
  mutate(species = replace(species, species == "Gunnel Sp", "UnID Gunnel")) %>% 
  mutate(species = replace(species, species == "Jellyfish Sp", "UnID Jellyfish")) %>% 
  mutate(ipa = as.factor(ipa)) %>%
  select(-c(transect_notes))

net <- net_import %>%
  mutate(species = replace_na(species, "none"),length_mm = replace_na(length_mm, 0)) %>%
  group_by(year, month, day, site, ipa, station, org_type, tax_group,species) %>% 
  mutate(count = 1) %>%
  summarize(species_count=sum(count), mean_length_mm = mean(length_mm)) %>%
  ungroup() 

#write to csv
write_csv(net, here("data","net_18.19.csv")) 

###############################################################################
## 2021 data

net_import2 <- drive_get("lampara_net_data_21") %>%
  read_sheet( ) 
net_import2_june <- drive_get("june_blitz_2021") %>%
  read_sheet( ) 

net_import2 <- rbind(net_import2, net_import2_june)

# write_csv(net_import2, here("data", "raw", "raw.21.csv"))

net_2021 <- net_import2 %>%
  mutate(month = str_pad(month, width = 2, pad = "0")) %>% 
  mutate(day = str_pad(day, width = 2, pad = "0")) %>% 
  mutate(length_mm = ifelse(is.na(length_mm), 0, length_mm)) %>% 
  mutate(species = ifelse(is.na(species), "none", species)) %>% 
  group_by(year, month, day, site, ipa, station, org_type, tax_group, species) %>%
  summarize(species_count = n(), mean_length_mm = mean(length_mm)) %>% 
  ungroup() 

#write to csv
write_csv(net_2021, here("data","net_2021.csv"))

###############################################################################
## 2022 data
net_import3 <- drive_get("lampara_net_data_22") %>%
  read_sheet( ) 

# write_csv(net_import3, here("data", "raw", "raw.22.csv"))

net_2022 <- net_import3 %>%
  mutate(month = str_pad(month, width = 2, pad = "0")) %>% 
  mutate(day = str_pad(day, width = 2, pad = "0")) %>% 
  mutate(length_mm = ifelse(is.na(length_mm), 0, length_mm)) %>% 
  mutate(species = ifelse(is.na(species), "none", species)) %>% 
  group_by(year, month, day, site, ipa, station, org_type, tax_group, species) %>%
  summarize(species_count = n(), mean_length_mm = mean(length_mm)) %>% 
  ungroup() 

#write to csv
write_csv(net_2022, here("data","net_2022.csv"))
