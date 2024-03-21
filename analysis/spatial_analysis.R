library(tidyverse)
library(here)
library(sf)
library(leaflet)
library(terra)

#### load data ####
national_HU12 <- here("data","spatial","NWBD_HUC_12-Digit_Basins_of_WA","WBDHU12.shp") %>% 
  read_sf()

shoreline <- here("data","spatial","shorezone_shoreline_only","shorezone_shoreline_only.shp") %>% 
  read_sf() %>% 
  st_transform(crs = 2927) 

site_levels <- c("FAM", "TUR", "COR", "MA", "WA", "HO", "SHR", "DOK", "LL", "TL", "PR", "EDG")

#GPS locations for our survey stations with each ipa
SOS_sites <- here("data","spatial","reupdatingthearmoringgeodatabase", "Shoreline_armoring_shore_origin_sites_UTM.shp") %>% 
  read_sf() %>% #UTM zone 10 32610
  st_transform(crs = 2927)  %>% #transform site coordinates into the same crs as the shoreline layer
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
#note this number matches the number of sites because while COR overlaps two HUCs, there are two sites on Vashon

SOS_HUCS <- SOS_HUCS_all %>% 
  select(c(HUC12, Name, geometry)) 

#### map it ####

SOS_site_cents_sn$lon <- st_coordinates(SOS_site_cents_sn)[,1]
SOS_site_cents_sn$lat <- st_coordinates(SOS_site_cents_sn)[,2]

leaflet(st_transform(SOS_site_cents_sn, crs = 4326)) %>%
  addProviderTiles(providers$Esri.WorldGrayCanvas, group =  "Esri") %>%
  setView(lng =-122.420429, lat = 47.886010, zoom = 8) %>%
  addPolylines(data = st_transform(national_HU12, crs = 4326), popup = national_HU12$Name) %>%
  addPolylines(data = st_transform(SOS_HUCS, crs = 4326), color = "red", label = ~Name) %>%
  addCircleMarkers(lng = ~lon, lat = ~lat, color = "purple")

####################
#### format to add SOS sites ####
# #create centroid points to make ordering N to S easier
# SOS_HUC_cents <- st_centroid(SOS_HUCS) %>% 
#   cbind(st_coordinates(.)) 
# 
# #this is only working in base r for some reason
# SOS_HUC_cents <- SOS_HUC_cents[order(SOS_HUC_cents$Y, decreasing = T),]
# 
# #reorder from north to south and add SOS site names
# SOS_HUCS_full <- SOS_HUCS %>% 
#   slice(rep(1:n(), times = c(1,1,1,4,1,1,1,2)))  %>% #multiple sites in some HUC12s 
#   mutate(Name = factor(Name, levels = unique(SOS_HUC_cents$Name))) %>% 
#   arrange(Name) %>% 
#   mutate(site = site_levels, .after = "Name")

#### load C-CAP data ####

#load 1m resolution impervious surface data
# ccap_2023imperv <- rast(here("data", "spatial", "wa_2021_ccap_v2_hires_impervious_20231119", "wa_2021_ccap_v2_hires_impervious_20231119.tif"))
# 
SOS_HUCS_proj<- st_transform(SOS_HUCS, crs = 5070)
# 
# #extract imperv by HUCs
# start.time <- Sys.time()
# ccap_imperv_HUCs <- terra::extract(ccap_2023imperv, SOS_HUCS_proj) #this takes forever
# end.time <- Sys.time()
# time.taken <- round(end.time - start.time,2)
# time.taken #36.39 minutes
# imperv_table <- table(ccap_imperv_HUCs)
# save(imperv_table, file = here("data", "imperv.table.Rdata"))
load(here("data", "imperv.table.Rdata")) #load the table instead of running the extract again

env_table <- as_tibble(imperv_table) %>% 
  pivot_wider(names_from = wa_2021_ccap_v2_hires_impervious_20231119, values_from = n) %>% 
  rename(natural = "0", imperv = "1") %>% 
  cbind(SOS_HUCS_proj) %>% 
  select(!c(ID,geometry)) %>% 
  select(HUC12, Name, natural, imperv) %>%
  mutate(perc.imperv23 = round(100*(imperv/(natural + imperv)),2)) %>% 
  mutate(perc.natural23 = 100 - perc.imperv23)

#load 30m resolution land cover data
#load 1m resolution impervious surface data
ccap_2016lc <- rast(here("data", "spatial", "2016_CCAP_Job987470", "2016_CCAP_J987470.tif"))

cover <- c("background", "unclassified", "high intensity developed",
           "medium intensity developed", "low intensity developed",
           "developed open space", "cultivated land", "pasture/hay", 
           "grassland", "deciduous forest", "evergreen forest",
           "mixed forest", "shrub/scrub", "palustrine forested wetland",
           "palustrine scrub/shrub wetland", "palustrine emergent wetland",
           "estuarine forested wetland", "estuarine scrub/shrub wetland",
           "estuarine emergent wetland", "unconsolidated shore", "bare land",
           "open water", "palustrine aquatic bed", "estuarine aquatic bed","tundra", "snow/ice")

levels(ccap_2016lc) <- data.frame(id=0:25, cover=cover)

start.time <- Sys.time()
ccap_lc_HUCs <- terra::extract(ccap_2016lc, SOS_HUCS_proj) 
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken # <5 seconds 

lc_table <- table(ccap_lc_HUCs) %>% 
  as_tibble() %>% 
  filter(!n == 0) %>% 
  pivot_wider(names_from = cover, values_from = n, values_fill = 0) %>% 
  column_to_rownames(var="ID") %>% 
  mutate(rowsum = rowSums(.))


