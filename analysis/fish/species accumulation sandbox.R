net_2018.19 <- here::here("data","net_18.19.csv") %>% 
  read_csv()

net_2021 <- here::here("data","net_2021.csv") %>% 
  read_csv()

net_2022 <- here::here("data","net_2022.csv") %>% 
  read_csv()

net_tidy <- bind_rows(net_2018.19, net_2021, net_2022) %>% 
  filter(!(site == "TUR" & ipa == "Natural")) %>% 
  mutate(ipa = replace(ipa, site == "TUR" & ipa == "Restored", "Natural")) %>% #no restoration at Turn Island, we sampled at 2 natural sites. Keeping the "second" one ("restored" in the field data) because it was an easier site to sample with less variability
  mutate(month = recode(month, `04` = "Apr", `05` = "May", `06` = "Jun", `07` = "Jul", `08` = "Aug", `09` = "Sept")) %>% 
  mutate(year = factor(year), 
         month = factor(month, levels = c("Apr", "May", "Jun", "Jul", "Aug", "Sept")),
         site = factor(site, levels = c("FAM", "TUR", "COR", "MA", "WA", "HO", "SHR", "DOK", "LL", "TL", "PR", "EDG"))) %>% 
  mutate(species_count = as.numeric(species_count)) %>%
  mutate(species_count = ifelse(is.na(org_type), 0, species_count)) 

sampling_events <- net_tidy %>% 
  select(year, month, day, site, ipa, station) %>% 
  distinct() %>% 
  filter(!(site == "COR" & year == "2019" & month == "Apr" & day == "30")) %>% #these COR seem like data entry errors because they were only sampled at one station and we already had complete sampling at COR in april and may. No fish recorded in either entry
  filter(!(site == "COR" & year == "2019" & month == "May" & day == "01")) %>% 
  filter(!(site == "TUR" & year == "2018" & month == "Jul" & day == "11")) %>%  #this is an incomplete sampling event. Complete sampling at TUR occured on 7/12/18
  filter(!(site == "TUR" & year == "2018" & month == "Sept" & day == "11")) #another month with repeat sampling; keeping only the second september TUR sampling event

#### determine which site/days were unevenly sampled
unbalanced_events <- sampling_events %>% 
  group_by(year, month, day, site) %>% 
  summarize(n_samples = n()) %>% 
  filter(n_samples < 9) %>% #sampling days when we missed either a depth or a shoreline (ipa)
  filter(!site == "TUR") #all of the TUR samples are balanced, but there are only 6 because there is no Restored ipa

balanced_events <- sampling_events %>%
  anti_join(unbalanced_events, join_by("year", "month", "day", "site"))

net_tidy2 <- net_tidy %>% 
  semi_join(sampling_events, join_by("year", "month", "day", "site")) %>%  #keep only the legit sampling events (days)
  anti_join(unbalanced_events, join_by("year", "month", "day", "site")) %>% #then remove unbalanced survey days
  filter(org_type == "Fish") %>% #### this is where you lose zero counts, if those are needed later (can't just add |is.na(org_type) because some tows had just jellies)
  full_join(balanced_events) %>% #add back in sampling events with zeros so they're represented in the species accumulation chart
  arrange(year, month, day, site, ipa, station) %>% 
  mutate(species_count=replace_na(species_count, 0)) %>% 
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
                             TRUE ~ ComName)) %>% 
  filter(!grepl("UnID", ComName))   #remove unidentified species

sampling_event18 <- sampling_events %>% 
  anti_join(unbalanced_events, join_by("year", "month", "day", "site")) %>%
  arrange(year, month, day) %>% 
  distinct(year, month, site, ipa) %>% 
  group_by(site, ipa) %>% 
  filter(row_number() < 19) #keep the first 18 samples from each ipa

view_cropped_effort <- sampling_event18 %>% 
  group_by(site, ipa) %>% 
  count

net_tidy18 <- net_tidy_b %>%
  semi_join(sampling_event18, join_by("year", "month", "site", "ipa"))

check_balance <- net_tidy18 %>%
  distinct(year, month, site, ipa) %>%
  group_by(site, ipa) %>%
  count

### species accumulation curve

make_filtered_spp_mat <- function(site_ID) {
  
  spp_mat <- net_tidy18 %>% 
    mutate(year_month = ym(paste(year, month, sep = "-"))) %>% 
    mutate(site_ipa = paste(site, ipa, sep = "_")) %>% 
    filter(site_ipa == site_ID) %>% 
    group_by(ComName, year_month) %>% 
    summarize(spp_sum = sum(species_count)) %>% 
    pivot_wider(names_from = ComName, values_from = spp_sum) %>% 
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

sample_ids <- paste(SOS_core_sites, rep(c("Armored", "Natural", "Restored"), each = 6), sep = "_") %>% 
  as_tibble() %>% 
  filter(!value == "TUR_Restored") %>% 
  pull()
  
curve_list <- lapply(sample_ids, make_spp_curve)
curve_df <- data.frame() 

for (i in 1:length(sample_ids)) {
  sites <- curve_list[[i]]$sites
  richness <- curve_list[[i]]$richness
  tmp <- data.frame(site = sample_ids[i], samples = sites, richness = richness)
  
  curve_df <- rbind(curve_df, tmp)
}

ggplot() +
  geom_point(data = curve_df, aes(x = samples, y = richness, color = site)) +
  geom_line(data = curve_df, aes(x = samples, y = richness, color = site)) +
  theme_classic() +
  labs(x = "Times sampled", y = "Number of species", color = "Site")

