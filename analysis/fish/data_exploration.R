library(tidyverse)
library(here)
library(lme4)
library(sjPlot)

load_data <- function() {
  #Load data
  #load field data (raw data downloaded/formatted in "import_data.R")
  net_2018.19 <- here::here("data","net_18.19.csv") %>% 
    read_csv()
  
  net_2021 <- here::here("data","net_2021.csv") %>% 
    read_csv()
  
  net_2022 <- here::here("data","net_2022.csv") %>% 
    read_csv()
  
  net_tidy <- bind_rows(net_2018.19, net_2021, net_2022) %>% 
    filter(!ipa == "Armored_2") %>% #remove second armored site from Titlow in 2021
    mutate(ipa = replace(ipa, site == "TUR" & ipa == "Restored", "Natural")) %>% #no restoration at Turn Island
    # mutate(site_ipa = paste(site, ipa, sep = "_")) %>% 
    mutate(date = make_date(year, month, day)) %>% 
    mutate(species_count = as.numeric(species_count)) %>%
    mutate(species_count = ifelse(is.na(org_type), 0, species_count)) %>% 
    filter(org_type == "Fish") %>% #### this is where you lose zero counts, if those are needed later
    rename(ComName = "species") %>% 
    mutate(ComName = case_when(tax_group == "Salmon" ~ paste(ComName, "salmon"), 
                               ComName == "Sand Sole" ~ "Pacific sand sole",
                               ComName == "Herring" ~ "Pacific herring", 
                               ComName == "Snake Prickleback" ~ "Pacific snake prickleback",
                               ComName == "Speckled Sand Dab" ~ "Speckled sanddab",
                               ComName == "Poacher" ~ "UnID Poacher",
                               ComName == "Striped Perch" ~ "Striped Seaperch",
                               ComName == "Staghorn Sculpin" ~ "Pacific staghorn sculpin",
                               ComName == "Stickleback" ~ "Three-spined stickleback",
                               ComName == "Tom Cod" ~ "Pacific tomcod",
                               ComName == "Sand Lance" ~ "Pacific sand lance",
                               ComName == "Pipefish" ~ "Bay pipefish",
                               ComName == "Irish Lord" ~ "UnID Irish Lord",
                               ComName == "Midshipman" ~ "Plainfin midshipman",
                               ComName == "White Spotted Greenling" ~ "Whitespotted greenling",
                               ComName == "Silver Spotted Sculpin" ~ "Silverspotted sculpin",
                               ComName == "Tube Snout" ~ "Tube-snout",
                               TRUE ~ ComName)) %>% 
    mutate(ComName = case_when(ComName == "Steelhead salmon" ~ "Steelhead trout",
                               ComName == "Cutthroat salmon" ~ "Cutthroat trout",
                               TRUE ~ ComName)) 
  return(net_tidy)
}

net_tidy <- load_data()

spp_richness <- net_tidy %>% 
  mutate(species_count = ifelse(is.na(species_count), 0, species_count)) %>% 
  group_by(year, month, site, ipa, station) %>% 
  summarize(richness = n_distinct(ComName, na.rm = FALSE)) %>% 
  ungroup() %>% 
  complete(nesting(year, month, site, ipa), station) %>% #use nesting to only include combinations already present in the data
  replace(is.na(.),0) %>% 
  mutate(ipa = factor(ipa, levels = c("Natural", "Armored", "Restored")), 
         station = factor(station, labels = c("shallow", "mid", "deep"))) %>% 
  mutate(season = ifelse(month %in% c("04", "05", "06"), "spring", "summer"))


ggplot(spp_richness) + 
  geom_violin(aes(x = station, y = richness))

ggplot(spp_richness) + 
  geom_violin(aes(x = ipa, y = richness))

ggplot(spp_richness) + 
  geom_violin(aes(x = site, y = richness))

ggplot(spp_richness) + 
  geom_violin(aes(x = season, y = richness))

summary(aov(richness ~ station, data = spp_richness))
summary(aov(richness ~ ipa, data = spp_richness))
summary(aov(richness ~ site, data = spp_richness))

station_mod <- lmer(richness ~ station + (1|ipa:site) + (1|site), data = spp_richness)
summary(station_mod)

plot_model(station_mod)
tab_model(station_mod)

season_mod <- lmer(richness ~ season + (1|site), data = spp_richness)
summary(season_mod)

plot_model(season_mod)
tab_model(season_mod)
