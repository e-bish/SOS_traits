# Load libraries
library(tidyverse)
library(here)
library(stringr)
library(googlesheets4)
library(googledrive)


#import data from drive
# wq_import <- drive_get("https://docs.google.com/spreadsheets/d/1yBiUKuI4AqL9VqItuFHUWVE_-_ph7fpp6jn2PibduAU/edit#gid=0") %>%
#   read_sheet() %>%
#   rename(notes = "...10")
# 
# #write to csv
# write_csv(wq_import, here("data","wq_import.csv")) 

#read csv
wq_import <- here::here("data","wq_import.csv") %>% read_csv()

wq_df <- wq_import %>% 
  mutate(ipa = replace(ipa, site == "TUR" & ipa == "Restored", "Natural")) %>%  #no restoration at Turn Island
  mutate(site = ifelse(site == "MAY", "MA", site)) %>% 
  mutate(secchi_depth_m = replace(secchi_depth_m, secchi_depth_m =="NULL", NA)) %>% 
  mutate(secchi_depth_m = as.numeric(unlist(secchi_depth_m))) %>% 
  suppressWarnings() %>% 
  select(!notes) %>% 
  pivot_longer(c(secchi_depth_m, do_mg_l, salinity_ppm, temperature), names_to = "metric") %>% 
  na.omit()

#rough n of wq measurements by site
wq_df %>% filter(metric == "salinity_ppm" & ipa == "Armored") %>% count(site)

#inspect variation by site and ipa
wq_df %>% 
  ggplot() +
  geom_boxplot(aes(x = site, y = value)) +
  facet_wrap(~metric)
#salinity is the most interesting metric here

wq_df %>% 
  ggplot() +
  geom_boxplot(aes(x = ipa, y = value)) +
  facet_wrap(~metric)
#nothing to see here

wq_metrics <- unique(wq_df$metric)

plot_metric <- function(metric_id) {
  wq_df %>% 
    filter(metric == metric_id) %>% 
    ggplot() +
    geom_boxplot(aes(x = site, y = value, fill = ipa)) +
    labs(title = metric_id) +
    theme_bw()
}

map(wq_metrics, plot_metric)
#some differences between some sites, particularly jubilee sites (maybe only should look at those with other june data ~ spring outflow?)

plot_time_series <- function(metric_id, jubilee) {
  
  make_plots <- function(df) {
    df %>% 
      filter(metric == metric_id) %>% 
      ggplot() +
      geom_point(aes(x = month, y = value, color = site, shape = ipa)) +
      labs(title = metric_id) +
      theme_bw() +
      facet_wrap(~year)
  }

  if (jubilee == T) {
    wq_df %>% 
      filter(month == 6) %>% 
      mutate(month = as.factor(month)) %>% 
      make_plots()
  } else {
    wq_df %>% 
      make_plots()
  }
}

map(wq_metrics, plot_time_series, jubilee = T)


################ format the matrix

## this represents the mean but may want to change it depending on what abundance measure is used
## for example, if MaxN would we just take the values from those specific days?

env.data <- wq_df %>% 
  mutate(month = str_pad(month, 2, pad = "0")) %>% 
  mutate(month = if_else(site == "MA", "06", month)) %>% # we did a July 1st survey at Maylor that we want to count as a June survey
  mutate(site_month = paste(site, month, sep = "_")) %>% 
  group_by(site_month, metric) %>% 
  summarize(mean = mean(value, na.rm = TRUE)) %>% 
  pivot_wider(names_from = metric, values_from = mean) %>% 
  column_to_rownames(var = "site_month")

save(env.data, file = here("data", "env.data.Rdata"))

# filter to june, this version simply takes the mean of all june data for each site
env.data.june <- wq_df %>% 
  mutate(month = str_pad(month, 2, pad = "0")) %>% 
  mutate(month = if_else(site == "MA", "06", month)) %>% # we did a July 1st survey at Maylor that we want to count as a June survey
  filter(month == "06") %>% 
  group_by(site, metric) %>% 
  summarize(mean = mean(value, na.rm = TRUE)) %>% 
  pivot_wider(names_from = metric, values_from = mean) %>% 
  # mutate(veg = recode(site, 
  #                     'COR' = "Present", 
  #                     'TUR'="Present", 
  #                     'FAM'="Absent", 
  #                     'DOK'="Absent", 
  #                     'EDG'="Absent", 
  #                     'SHR'="Present", 
  #                     'HO' = "Present", 
  #                     'LL' = "Absent", # overall absent though armored site may have had some eelgrass, based on field notes
  #                     'MA' = "Absent", 
  #                     'PR' = "Present",
  #                     'TL' = "Present", 
  #                     'WA' = "Present")) %>% 
  column_to_rownames(var = "site")

save(env.data.june, file = here("data", "env.data.june.Rdata"))

