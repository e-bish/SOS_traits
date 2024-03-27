library(tidyverse)
library(here)
library(sf)
library(leaflet)
library(terra)

#### load data ####
national_HU12 <- here("data","spatial","NWBD_HUC_12-Digit_Basins_of_WA","WBDHU12.shp") %>% 
  read_sf()

site_levels <- c("FAM", "TUR", "COR", "MA", "WA", "HO", "SHR", "DOK", "LL", "TL", "PR", "EDG")

#GPS locations for our survey stations with each ipa
SOS_sites <- here("data","spatial","reupdatingthearmoringgeodatabase", "Shoreline_armoring_shore_origin_sites_UTM.shp") %>% 
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
  mutate(site = factor(site, levels = site_levels)) 

#create a buffer around sites that will overlap with the HUCs
site_buffers <- st_buffer(SOS_sites, 1000)

#remove HUCs in the Salish Sea that only represent water
land_HU12 <- national_HU12 %>% 
  filter(!Name %in% c("Puget Sound", "Skagit Bay", "Haro Strait-Strait of Georgia", "South Puget Sound"))

SOS_HUCS_all <- st_filter(land_HU12, site_buffers)
#note this number matches the number of sites because while there are two sites on Vashon, COR overlaps two HUCs

SOS_HUCS <- SOS_HUCS_all %>% 
  select(c(HUC12, Name, geometry))  

#create centroid points to make ordering N to S easier
SOS_HUC_cents <- st_centroid(SOS_HUCS) %>% 
  cbind(st_coordinates(.)) %>% 
  st_drop_geometry() %>% 
  mutate(Y = ifelse(grepl("Tacoma", Name), Y + 50000, Y)) %>% 
  arrange(desc(Y))

#reorder from north to south and add SOS site names
SOS_HUCS <- SOS_HUCS %>% 
  mutate(Name = factor(Name, levels = unique(SOS_HUC_cents$Name))) %>% 
  arrange(Name)

#### map it ####

SOS_sites_transformed <- SOS_sites %>% 
  st_zm() %>% 
  st_transform(crs = 4326) %>% 
  mutate(lon = st_coordinates(.)[,1],
         lat = st_coordinates(.)[,2])

leaflet(SOS_sites_transformed) %>%
  addCircleMarkers(lng = ~lon, lat = ~lat, color = "purple") %>% 
  addProviderTiles(providers$Esri.WorldGrayCanvas, group =  "Esri") %>%
  setView(lng =-122.420429, lat = 47.886010, zoom = 8) %>%
  addPolylines(data = st_transform(national_HU12, crs = 4326), popup = national_HU12$Name) %>%
  addPolylines(data = st_transform(SOS_HUCS, crs = 4326), color = "red", label = ~Name)

#### load C-CAP data ####
#transform HUCs to the projection of the ccap data
SOS_HUCS_proj <- st_transform(SOS_HUCS, crs = 5070)

#load 30m resolution land cover data
ccap_2016lc <- rast(here("data", "spatial", "2016_CCAP_Job987470", "2016_CCAP_J987470.tif"))

#classify land cover data based on NOAA codes
cover <- c("background", "unclassified", "high_intensity_developed",
           "medium_intensity_developed", "low_intensity_developed",
           "developed_open_space", "cultivated_land", "pasture.hay", 
           "grassland", "deciduous_forest", "evergreen_forest",
           "mixed_forest", "shrub.scrub", "palustrine_forested_wetland",
           "palustrine_scrub.shrub_wetland", "palustrine_emergent_wetland",
           "estuarine_forested_wetland", "estuarine_scrub.shrub_wetland",
           "estuarine_emergent_wetland", "unconsolidated_shore", "bare_land",
           "open_water", "palustrine_aquatic_bed", "estuarine_aquatic_bed","tundra", "snow.ice")

levels(ccap_2016lc) <- data.frame(id=0:25, cover=cover)

#extract land cover by HUCs
start.time <- Sys.time()
ccap_lc_HUCs <- terra::extract(ccap_2016lc, SOS_HUCS_proj) 
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken # <5 seconds 

#format land cover data for each HUC
lc_table <- table(ccap_lc_HUCs) %>% 
  as_tibble() %>% 
  filter(!n == 0) %>% 
  pivot_wider(names_from = cover, values_from = n, values_fill = 0) %>% 
  column_to_rownames(var="ID") %>% 
  mutate(rowsum16 = rowSums(.))

#assign land cover values to proper HUCs
env_raw_table <- SOS_HUCS_proj %>% 
  st_drop_geometry() %>% 
  bind_cols(lc_table)

#load 1m resolution impervious surface data
# ccap_2023imperv <- rast(here("data", "spatial", "wa_2021_ccap_v2_hires_impervious_20231119", "wa_2021_ccap_v2_hires_impervious_20231119.tif"))
# #extract imperv by HUCs
# start.time <- Sys.time()
# ccap_imperv_HUCs <- terra::extract(ccap_2023imperv, SOS_HUCS_proj) #this takes forever
# end.time <- Sys.time()
# time.taken <- round(end.time - start.time,2)
# time.taken #36.39 minutes
# imperv_table <- table(ccap_imperv_HUCs)
# save(imperv_table, file = here("data", "imperv.table.Rdata"))
load(here("data", "imperv.table.Rdata")) #load the table instead of running the extract again

#reformat impervious surface data by HUC
imperv_tibble <- as_tibble(imperv_table) %>% 
  pivot_wider(names_from = wa_2021_ccap_v2_hires_impervious_20231119, values_from = n) %>% 
  rename(natural23 = "0", imperv23 = "1") %>% 
  select(!ID) %>% 
  mutate(rowsum23 = rowSums(.)) 

#combine data from both sets of ccap data
env_raw_table <-  bind_cols(env_raw_table, imperv_tibble)

#add up the values for the two HUCs that Cornet spans
Cornet_Hucs <- env_raw_table %>% 
  filter(HUC12 ==  171100190101 | HUC12 == 171100190102) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  mutate(Name = "Cornet_HUC", .before = high_intensity_developed)

#remove HUC12s to combine data by site
SOS_env_table <- env_raw_table %>% 
  select(!HUC12) %>% 
  bind_rows(Cornet_Hucs) %>% 
  mutate(Name = factor(Name, levels = unique(SOS_HUC_cents$Name))) %>% 
  arrange(Name)

#combine data in the correct order
# SOS_env_table2 <- bind_rows(SOS_env_table[1:2,],
#                            Cornet_Hucs,
#                            SOS_env_table[5:8,], 
#                            SOS_env_table[9,], 
#                            SOS_env_table[9,], 
#                            SOS_env_table[10:12,]) %>% 
#   mutate(site = site_levels, .after = Name)
#   

env_table <- SOS_env_table2 %>% 
  rowwise() %>% 
  mutate(perc.developed = sum(high_intensity_developed, medium_intensity_developed, low_intensity_developed, developed_open_space)/rowsum16,
         perc.agriculture = sum(cultivated_land, pasture.hay)/rowsum16,
         perc.forest = sum(deciduous_forest, evergreen_forest, mixed_forest)/rowsum16,
         perc.shrub = sum(shrub.scrub)/rowsum16,
         perc.wetland = sum(palustrine_forested_wetland, palustrine_scrub.shrub_wetland, 
                            palustrine_emergent_wetland, estuarine_emergent_wetland)/rowsum16,
         perc.bare_land = sum(bare_land)/rowsum16) %>%
  mutate(perc.imperv23 = imperv23/(natural23 + imperv23)) %>% 
  mutate(perc.natural23 = 1 - perc.imperv23) %>% 
  select(!3:26)
  


