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
    mutate(month = if_else(site == "MA", "06", month)) %>% # we did a July 1st survey at Maylor that we want to count as a June survey
    mutate(ipa = replace(ipa, site == "TUR" & ipa == "Restored", "Natural")) %>% #no restoration at Turn Island
    mutate(date = make_date(year, month, day)) %>% 
    mutate(month = recode(month, `04` = "Apr", `05` = "May", `06` = "Jun", `07` = "Jul", `08` = "Aug", `09` = "Sept")) %>% 
    mutate(year = factor(year), 
           month = factor(month, levels = c("Apr", "May", "Jun", "Jul", "Aug", "Sept")),
           site = factor(site, levels = c("FAM", "TUR", "COR", "MA", "WA", "HO", "SHR", "DOK", "LL", "TL", "PR", "EDG"))) %>% 
    mutate(species_count = as.numeric(species_count)) %>%
    mutate(species_count = ifelse(is.na(org_type), 0, species_count)) 
  
  #check that we're keeping only the data we want
  sampling_events <- net_tidy %>% 
    select(year, month, day, site, ipa, station) %>% 
    distinct() %>% 
    filter(!(site == "COR" & year == "2019" & month == "Apr" & day == "30")) %>% #these COR seem like data entry errors because they were only sampled at one station and we already had complete sampling at COR in april and may. No fish recorded in either entry
    filter(!(site == "COR" & year == "2019" & month == "May" & day == "01")) %>% 
    filter(!(site == "TUR" & year == "2018" & month == "Jul" & day == "11")) %>%  #this is an incomplete sampling event. Complete sampling at TUR occured on 7/12/18
    filter(!(site == "TUR" & year == "2018" & month == "Sept" & day == "11")) #another month with repeat sampling; keeping only the second september TUR sampling event
  
  net_tidy2 <- net_tidy %>% 
    semi_join(sampling_events, join_by("year", "month", "day", "site")) %>%  #keep only the legit sampling events (days)
    filter(org_type == "Fish") %>% #### this is where you lose zero counts, if those are needed later
    full_join(sampling_events) %>% #add back in sampling events with zeros so they're represented in the species accumulation chart
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
    filter(!grepl("UnID", ComName)) #remove unidentified species
  
  return(net_tidy2)
}

net_tidy <- load_data()

#### explore the data ####

##what were the most common species we caught?
net_tidy %>%
  group_by(ComName) %>%
  summarize(total = sum(species_count)) %>% 
  arrange(-total)

##how often did we encounter each species?
times_encountered <- net_tidy %>% 
  group_by(ComName, tax_group, site) %>% 
  summarize(freq_obs = n())

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

spp_by_site <- net_tidy %>% 
  count(ComName, site)

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

abund_by_site <- net_tidy %>%
  group_by(site) %>%
  summarize(total = sum(species_count)) 

cor(n_spp_site$n_spp, abund_by_site$total) #yes

##did we sample enough to adequately describe the community?
make_filtered_spp_mat <- function(site_ID) {
  
  spp_mat <- net_tidy %>% 
    mutate(year_month = ym(paste(year, month, sep = "-"))) %>% 
    filter(site == site_ID) %>% 
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

curve_list <- lapply(SOS_sites, make_spp_curve)
curve_df <- data.frame() 

for (i in 1:length(SOS_sites)) {
  sites <- curve_list[[i]]$sites
  richness <- curve_list[[i]]$richness
  tmp <- data.frame(site = SOS_sites[i], samples = sites, richness = richness)
  
  curve_df <- rbind(curve_df, tmp)
}

ggplot() +
  geom_point(data = curve_df, aes(x = samples, y = richness, color = factor(site, level = SOS_sites))) +
  geom_line(data = curve_df, aes(x = samples, y = richness, color = factor(site, level = SOS_sites))) +
  theme_classic() +
  labs(x = "Times sampled", y = "Number of species", color = "Site")

ggsave("docs/figures/specaccum_plot.png")

##is there obvious seasonality in our catch?
#in total catch abundance?
net_tidy %>% 
  filter(site %in% SOS_core_sites) %>% #remove jubilee sites to properly compare June
  ggplot(aes(x = month, y = species_count, fill = year)) +
  geom_bar(stat = "identity") +
  labs(x = "Month", y = "Abundance", fill = "Year") + 
  theme_classic()
#doesn't seem to be a clear break between "early" and "late" 

#in species richness?
n_spp_by_month <- net_tidy %>% 
  filter(site %in% SOS_core_sites) %>% #remove jubilee sites to properly compare June
  group_by(year, month) %>% 
  summarize(n_spp = n_distinct(ComName)) 

n_spp_by_month %>% 
  ggplot(aes(x = month, y = n_spp, group = year, color = year)) +
  geom_line() +
  geom_point() + 
  theme_classic() + 
  labs(x = "Month", y = "Species Richness", color = "Year")
#again, doesn't seem to be a clear break between "early" and "late" 

##is there a difference between northern and southern sites?
abund_region <- net_tidy %>% 
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
net_tidy <- net_tidy %>% 
  filter(site %in% SOS_core_sites) #use only core sites because we didn't sample jubilee sites enough to capture the community

save(net_tidy, file = here("data", "net_tidy.Rdata"))



