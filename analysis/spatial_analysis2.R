library(tidyverse)
library(here)
library(sf)
library(terra)

#### load data ####
core_site_names <- c("FAM", "TUR", "COR", "SHR", "DOK", "EDG")

#GPS locations for our survey stations
SOS_core_sites <- here("data","spatial","reupdatingthearmoringgeodatabase", "Shoreline_armoring_shore_origin_sites_UTM.shp") %>% 
  read_sf() %>% #UTM zone 10 32610
  st_transform(crs = 2927)  %>% #transform site coordinates into the same crs as the HUC layer
  mutate(site = case_when(Site_name == "Dockton Park" ~ "DOK",
                          Site_name == "Cornet Bay" ~ "COR",
                          Site_name == "Edgewater Beach" ~ "EDG",
                          Site_name == "Family Tides" ~ "FAM",
                          Site_name == "Seahurst Park" ~ "SHR", 
                          Site_name == "Turn Island" ~ "TUR",
                          Site_name == "Lost Lake" ~ "LL",
                          Site_name == "Titlow Park" ~ "TL",
                          Site_name == "Penrose Point" ~ "PR",
                          Site_name == "Waterman Shoreline Preserve" ~ "WA" ,
                          Site_name == "Howarth Park" ~ "HO" ,
                          Site_name == "Maylor Point" ~ "MA"), .after = "Site_name") %>% 
  filter(site %in% core_site_names)
  mutate(site = factor(site, levels = core_site_names)) %>% 

#load percent armor data 
perc_armor <- here("data", "perc_armor.csv") %>%  #this is pre-2020! Consider also using the updated one with dockton 2020 restoration
  read_csv() %>% 
  filter(site %in% core_site_names) %>% 
  mutate(site = factor(site, levels = core_site_names)) %>% 
  arrange(site) %>% 
  select(site, "500m", "10km")
  