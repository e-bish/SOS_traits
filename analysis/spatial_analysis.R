library(tidyverse)
library(here)
library(sf)
library(leaflet)

national_HU12 <- here("data","spatial","NWBD_HUC_12-Digit_Basins_of_WA","WBDHU12.shp") %>% 
  read_sf()

shoreline <- here("data","spatial","shorezone_shoreline_only","shorezone_shoreline_only.shp") %>% 
  read_sf() %>% 
  st_transform(crs = 2927) 

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
                          Site_name == "Maylor Point" ~ "MA"), .after = "Site_name")

#set a center point for each site 
SOS_site_cents <- SOS_sites %>%
  group_by(site) %>%
  summarize(geometry = st_union(geometry)) %>%
  st_centroid()

#snap it to the closest shoreline
closest_points <- SOS_site_cents %>%
  rowwise() %>%
  mutate(
    nearest_segment = shoreline[st_nearest_feature(geometry,
                                                   shoreline),],
    line_to_point = st_nearest_points(geometry, nearest_segment),
    closest_point = st_cast(line_to_point, 'POINT')[2],
    snapped_point_cond = st_sfc(st_geometry(closest_point),
                                crs = st_crs(shoreline)))

SOS_site_cents_sn <- closest_points %>% 
  select(site, snapped_point_cond) %>% 
  st_drop_geometry() %>% 
  st_as_sf()

#this extracts the right HUC12s but would have to manually assign them to sites
SOS_HUCS <- st_filter(national_HU12, SOS_site_cents_sn) %>%
  bind_rows(national_HU12 %>% filter(Name %in% "Vashon Island")) 


#this assigns the correct information but not the right geometry
# Vashon_HUCS <- national_HU12 %>% 
#   filter(Name %in% "Vashon Island") %>% 
#   bind_rows(., .) %>% 
#   mutate(site = c("DOK", "LL"), .after = "Name")
# 
# SOS_HUCS <- st_intersection(SOS_site_cents_sn[,2], national_HU12, sparse = F) %>% 
#   mutate(site = SOS_site_cents$site, .after = "Name") %>% 
#   filter(!site %in% c("DOK", "LL")) %>% 
#   bind_rows(Vashon_HUCS)
  

SOS_site_cents_sn$lon <- st_coordinates(SOS_site_cents_sn)[,1]
SOS_site_cents_sn$lat <- st_coordinates(SOS_site_cents_sn)[,2]

leaflet(st_transform(SOS_site_cents_sn, crs = 4326)) %>%
  addProviderTiles(providers$Esri.WorldGrayCanvas, group =  "Esri") %>%
  setView(lng =-122.420429, lat = 47.886010, zoom = 8) %>%
  addPolylines(data = st_transform(national_HU12, crs = 4326), label = ~Name) %>% 
  addPolylines(data = st_transform(SOS_HUCS, crs = 4326), color = "red", label = ~Name) %>% 
  addCircleMarkers(lng = ~lon, lat = ~lat, color = "purple")

# WA_HU12_l <- st_intersects(national_HU12, shoreline[,1], sparse = F) 
# WA_HU12 <- national_HU12[WA_HU12_l,]




