library(tidyverse)
library(here)
library(janitor)
library(vegan)

###### This creates a list of matrices for trait and abundance data only by sampling event ######

# If common_to_sci gives issues, install the duckdb package
# options(timeout=100)
# install.packages("duckdb", repos = c("https://duckdb.r-universe.dev", "https://cloud.r-project.org"))

SOS_sites <- factor(c("FAM", "TUR", "COR", "MA", "WA", "HO", "SHR", "DOK", "LL", "TL", "PR", "EDG"), 
                    levels = c("FAM", "TUR", "COR", "MA", "WA", "HO", "SHR", "DOK", "LL", "TL", "PR", "EDG"))
SOS_core_sites <- factor(c("FAM", "TUR", "COR", "SHR", "DOK", "EDG"), 
                         levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG"))
SOS_jubilee_sites <- factor(c("MA", "WA", "HO", "LL", "TL", "PR"),
                            levels = c("MA", "WA", "HO", "LL", "TL", "PR"))

load_data <- function() {
  #Load data
  #load field data (raw data downloaded/formatted in "01_import_data.R")
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
  
  net_tidy_balanced <- net_tidy %>% 
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
  
  net_tidy18 <- net_tidy_balanced %>%
    semi_join(sampling_event18, join_by("year", "month", "site", "ipa"))
  
  return(list(net_tidy_balanced, net_tidy18))
}

net_tidy_list <- load_data()

net_tidy_balanced <- net_tidy_list[[1]]
net_tidy18 <- net_tidy_list[[2]]

#### explore the data ####

##what were the most common species we caught?
net_tidy18 %>%
  group_by(ComName) %>%
  summarize(total = sum(species_count)) %>% 
  arrange(-total)

##how often did we encounter each species?
times_encountered <- net_tidy18 %>% 
  group_by(ComName, tax_group, site) %>% 
  summarize(freq_obs = n()) %>% 
  filter(!is.na(tax_group))

#by tax group
ggplot(times_encountered, aes(x = tax_group, y = freq_obs, fill = factor(site, levels = SOS_sites))) + 
  geom_bar(position = "stack", stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "", y = "Frequency observed", fill = "Site") 

ggplot(times_encountered, aes(x = factor(site, levels = SOS_sites), y = freq_obs, fill = tax_group)) + 
  geom_bar(position = "stack", stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "", y = "Frequency observed", fill = "tax group") 

##in how many sites does each species occur?
sites_encountered <- times_encountered %>% 
  mutate(pa = ifelse(freq_obs > 0, 1, 0)) %>% 
  group_by(ComName) %>% 
  summarize(encounters_site = sum(pa)) %>% 
  arrange(desc(encounters_site))

ggplot(sites_encountered, aes(x = reorder(ComName, -encounters_site), y = encounters_site)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs(x = "", y = "No. of sites encountered")

spp_by_site <- net_tidy18 %>% 
  count(ComName, site) %>% 
  filter(!is.na(ComName))

ggplot(spp_by_site, aes(x = reorder(ComName, -n), y = n, fill = site)) +
  geom_col() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#core sites only
ggplot(filter(spp_by_site, site %in% SOS_core_sites), aes(x = reorder(ComName, -n), y = n, fill = site)) +
  geom_col() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

spp_site_count <- spp_by_site %>% 
  # filter(site %in% SOS_core_sites) %>% 
  group_by(ComName) %>% 
  summarize(n_sites = n()) %>% 
  arrange(desc(n_sites))

##what is the abundance of each species when it occurs?
spp_by_site %>% 
  # filter(site %in% SOS_core_sites) %>% 
  ggplot(aes(x = reorder(ComName, -n), y = n)) +
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

##is the mean abundance correlated with the number of sites where it occurs?
mean_count <- spp_by_site %>% 
  group_by(ComName) %>% 
  summarize(mean_n = mean(n))

cor(mean_count$mean_n, spp_site_count$n_sites) #no

##is the total abundance of fish correlated with the number of species in each site?
n_spp_site <- spp_by_site %>% 
  group_by(site) %>% 
  summarize(n_spp = n())

abund_by_site <- net_tidy18 %>%
  group_by(site) %>%
  summarize(total = sum(species_count)) 

cor(n_spp_site$n_spp, abund_by_site$total) #yes

##did we sample enough to adequately describe the community?
make_filtered_spp_mat <- function(site_ID, df) {
  
  spp_mat <- df %>% 
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

make_spp_curve <- function(site_ID, df) {
  
  site_mat <- make_filtered_spp_mat(site_ID, df)
  site_curve <- specaccum(site_mat, method = "collector", permutations = 100) #collector method preserves the order
  return(site_curve)
  
}

sample_ids <- paste(SOS_core_sites, rep(c("Armored", "Natural", "Restored"), each = 6), sep = "_") %>% 
  as_tibble() %>% 
  filter(!value == "TUR_Restored") %>% 
  pull()

## for the full balanced dataset ##
curve_list <- lapply(sample_ids, make_spp_curve, df = net_tidy_balanced)
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

# ggsave("docs/figures/specaccum_plot_bal.png")

## for the first 18 samples ##
curve_list <- lapply(sample_ids, make_spp_curve, df = net_tidy18)
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

# ggsave("docs/figures/specaccum_plot18.png")

##is there obvious seasonality in our catch?
#in total catch abundance?
net_tidy18 %>% 
  filter(site %in% SOS_core_sites) %>% #remove jubilee sites to properly compare June, though there still isn't total balance between months
  ggplot(aes(x = month, y = species_count, fill = year)) +
  geom_bar(stat = "identity") +
  labs(x = "Month", y = "Abundance", fill = "Year") + 
  theme_classic()
#doesn't seem to be a clear break between "early" and "late" 

#in species richness?
n_spp_by_month <- net_tidy18 %>% 
  filter(site %in% SOS_core_sites) %>% #remove jubilee sites to properly compare June, though there still isn't total balance between months
  group_by(year, month) %>% 
  summarize(n_spp = n_distinct(ComName)) 

net_tidy_balanced %>% 
  filter(site %in% SOS_core_sites) %>% #remove jubilee sites to properly compare June, though there still isn't total balance between months
  group_by(year, month) %>% 
  summarize(n_spp = n_distinct(ComName))  %>% 
  ggplot(aes(x = month, y = n_spp, group = year, color = year)) +
  geom_line() +
  geom_point() + 
  theme_classic() + 
  labs(x = "Month", y = "Species Richness", color = "Year")
#again, doesn't seem to be a clear break between "early" and "late" except that things die down in August

##is there a difference between northern and southern sites?
abund_region <- net_tidy18 %>% 
  mutate(year_month = ym(paste(year, month, sep = "-"))) %>% 
  mutate(region = case_when(site %in% c("FAM", "TUR", "COR", "MA", "HO", "WA") ~ "North",
                            TRUE ~ "South")) %>% 
  group_by(region, site, year_month) %>% 
  summarize(n_spp = n_distinct(ComName))

abund_region %>% 
  ggplot(aes(x = site, y = n_spp, fill = region)) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "Site", y = "Species Richness", fill = "Region")

#### final tidy dataframe for analysis ####
net_tidy.balanced <- net_tidy18 %>% 
  filter(site %in% SOS_core_sites) #use only core sites because we didn't sample jubilee sites enough to capture the community

save(net_tidy.balanced, file = here("data", "net_tidy.balanced.Rdata"))



