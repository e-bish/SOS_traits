library(tidyverse)
library(here)
library(janitor)
library(rfishbase)

############################## this was a test but ultimately doesn't work because some species have all zeros
############################## which throws an error when you run FD

# If common_to_sci gives issues, install the duckdb package
# options(timeout=100)
# install.packages("duckdb", repos = c("https://duckdb.r-universe.dev", "https://cloud.r-project.org"))

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

# sum_counts <- net_tidy %>% 
#   group_by(ComName) %>% 
#   summarize(total = sum(species_count)) %>% 
#   arrange(-total)

################################################################################
#prepare field abundance data for analysis
create_fish_matrices <- function(net_tidy) {
  
fish_N <- net_tidy %>% 
  mutate(month = factor(month, labels = c("Apr", "May", "Jun", "Jul", "Aug", "Sept"))) %>% 
  select(month, site, ComName, species_count) %>% #fish counts summed by site/day
  arrange(ComName) %>% 
  mutate(month = replace(month, site == "MA", "Jun")) # we did a July 1st survey at Maylor that we want to count as a June survey

spp_names <- fish_N %>% 
  distinct(ComName) %>% 
  mutate(Species = NA) 

sci_names <- vector(mode = 'list', length = length(spp_names))

for (i in 1:nrow(spp_names)) {
  sci_names[[i]] <- rfishbase::common_to_sci(spp_names[i,])
  spp_names[i,2] <- ifelse(nrow(sci_names[[i]]) == 1, sci_names[[i]][[1]], NA)
}

spp_names <- spp_names %>% 
  mutate(Species = case_when(ComName == "Pacific herring" ~ "Clupea pallasii", 
                           ComName == "Pacific Cod" ~ "Gadus macrocephalus",
                           ComName == "Northern Anchovy" ~ "Engraulis mordax",
                           ComName == "Striped Seaperch" ~ "Embiotoca lateralis",
                           ComName == "Rock Sole" ~ "Lepidopsetta bilineata",
                           ComName == "Steelhead trout" ~ "Oncorhynchus mykiss",
                           ComName == "Cutthroat trout" ~ "Oncorhynchus clarkii",
                           ComName == "Tube-snout" ~ "Aulorhynchus flavidus",
                           TRUE ~ Species)) %>% 
  arrange(ComName) %>% 
  filter(!str_detect(ComName, 'UnID'))

spp_names$Clean <- make_clean_names(spp_names$ComName)

# spp_names <- milieu %>% 
#   select(Species, SpecCode) %>% 
#   inner_join(spp_names) %>% 
#   mutate(Species2 = str_to_sentence(ComName))

Max_N <- fish_N %>% 
  group_by(month, site, ComName) %>%
  filter(species_count == max(species_count)) %>% 
  ungroup() %>% 
  rename(MaxN = "species_count") %>% 
  distinct(month, site, ComName, MaxN) %>% 
  filter(!str_detect(ComName, 'UnID'))

MaxN_list <- split(Max_N, Max_N$month)

prep_count_dfs <- function(df){
  
fish_counts <- df %>% 
  select(!month) %>% 
  add_row(site = NA, ComName = setdiff(spp_names$ComName, df$ComName), MaxN = 0) %>% 
  complete(site, ComName) %>% 
  filter(!is.na(site)) %>% 
  replace(is.na(.), 0) %>% 
  pivot_wider(names_from = ComName, values_from = MaxN) %>% 
  column_to_rownames(var="site") %>% 
  clean_names() %>% 
  as.matrix()

return(fish_counts)
}

prepped_MaxN_list <- lapply(MaxN_list, FUN = prep_count_dfs)

# trait data

fork_length <- net_tidy %>% 
  group_by(ComName) %>% 
  summarize(mean_length_mm = mean(mean_length_mm)) %>% 
  # mutate(fork_length = case_when(mean_fork_length < 70 ~ "small", ## does this need to be categorical??
  #                                mean_fork_length > 150  ~ "large",
  #                                TRUE ~ "medium")) %>% 
  # filter(ComName %in% spp_names$ComName) %>% 
  inner_join(spp_names) %>% 
  mutate(mean_length_mm = ifelse(ComName == "Tidepool Sculpin", 89.0, mean_length_mm)) #no length in our df so taking the max length from fishbase

schooling <- net_tidy %>% 
  mutate(schooling = ifelse(species_count > 10, "school", "nonschool")) %>% #arbitrary schooling cutoff
  filter(ComName %in% spp_names$ComName) %>% 
  group_by(ComName, schooling) %>% 
  summarize(count = n()) %>%  #see which spp had schooling and nonschooling results
  ungroup() %>% 
  pivot_wider(names_from = schooling, values_from = count) %>% 
  replace(is.na(.), 0) %>% 
  mutate(schooling = ifelse(school > 3, "schooling", "nonschooling")) %>% #arbitrary cuttoff
  mutate(schooling = ifelse(ComName == "Tube-snout", "schooling", schooling)) %>% 
  select(-c(nonschool, school))

milieu <- species(spp_names$Species) %>% #could also do length
  select(Species, BodyShapeI, DemersPelag, AnaCat) %>% 
  mutate(DemersPelag = ifelse(DemersPelag == "pelagic-neritic", "pelagic", DemersPelag)) %>% #simplify this category because there is only one pelagic and two pelagic-neritic species
  mutate(migrations = ifelse(is.na(AnaCat), "non-migratory", AnaCat)) %>% #presumed non migratory if no information is available
  select(!AnaCat)

feeding_guild1 <- fooditems(spp_names$Species) %>% 
  select(Species, FoodI, FoodII, PredatorStage) %>% 
  group_by(Species, FoodI) %>% 
  summarize(count = n()) %>% 
  group_by(Species) %>% 
  mutate(per =  100 *count/sum(count)) %>% 
  filter(per >= 60) %>% #classify main food source if more than 60% of diet is in one category
  inner_join(spp_names) %>% 
  arrange(ComName) %>% 
  mutate(feeding_guild = case_when(FoodI == "zoobenthos" ~ "Zoobenthivorous",
                                   FoodI == "zooplankton" ~ "Planktivorous", 
                                   FoodI == "nekton" ~ "Piscivorous",
                                   TRUE ~ FoodI)) %>% 
  ungroup() %>% 
  select(Species, feeding_guild)

feeding_guild2 <- fooditems(spp_names$Species) %>% 
  select(Species, FoodI, FoodII, PredatorStage) %>% 
  group_by(Species, FoodI) %>% 
  summarize(count = n()) %>% 
  group_by(Species) %>% 
  mutate(per = 100 *count/sum(count)) %>% 
  mutate(category = ifelse(per >= 60, FoodI, NA)) %>% 
  mutate(onlyNA = all(is.na(category))) %>% 
  filter(onlyNA == TRUE) %>% 
  inner_join(spp_names) %>% 
  select(!c(category, onlyNA)) %>% 
  mutate(feeding_guild = "Omnivorous") %>% 
  ungroup() %>% 
  select(Species, feeding_guild) %>% 
  distinct()

feeding_guild <- rbind(feeding_guild1, feeding_guild2) %>% 
  add_row(Species = "Blepsias cirrhosus", feeding_guild = "Zoobenthivorous") %>% #fishbase "Diet"
  add_row(Species = "Liparis florae", feeding_guild = "Zoobenthivorous")  %>% #don't have a great source for this one
  mutate(feeding_guild = str_to_lower(feeding_guild))

# morphs <- morphology(spp_names$Species) sparse information on mouth position
# swim <- swimming(spp_names$Species) less than half of species represented

fish_traits <- full_join(fork_length, milieu) %>% 
  select(3,1,2,4,5,6) %>% 
  left_join(schooling) %>% 
  left_join(feeding_guild) %>% 
  mutate(ComName = replace(ComName, ComName == "Pacific Sandfish", "Pacific sandfish")) %>%
  arrange(ComName)

# write_csv(fish_traits, "data/fish_traits.csv")

fish_trait_mat <- fish_traits %>% 
  select(-c(Species, ComName)) %>% 
  mutate_if(is.character, as.factor) %>% 
  clean_names() %>% 
  column_to_rownames(var = "clean") %>% 
  as.matrix()
  
return(list("traits" = fish_trait_mat, "abund" = prepped_MaxN_list))

}

fish.list <- create_fish_matrices(net_tidy)

save(fish.list, file = here("data", "fish.list.Rdata"))


