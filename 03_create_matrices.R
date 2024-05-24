################################################################################
#prepare field abundance data for analysis
create_fish_matrices <- function(net_tidy) {
  
  spp_names <- net_tidy %>% 
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
  
  fish_L_df <- net_tidy %>% #L is referring to the RLQ analysis
    group_by(year, month, site, ipa, ComName) %>%
    summarize(spp_sum = sum(species_count)) %>% 
    ungroup() %>% 
    complete(nesting(year, month, site), ipa, ComName, fill = list(spp_sum = 0)) %>% #fill in zeros for shorelines we sampled
    group_by(year, month, site, ComName) %>% 
    summarize(spp_mean = round(mean(spp_sum), 2)) %>% 
    pivot_wider(names_from = ComName, values_from = spp_mean, values_fill = 0) %>% 
    clean_names() %>% 
    ungroup()
  
  fish_L_mat <- fish_L_df %>% 
    mutate(sample = paste(year, month, site, sep = "_"), .after = site) %>% 
    select(!1:3) %>% 
    column_to_rownames(var = "sample")
  
  # trait data
  
  fork_length <- net_tidy %>% 
    group_by(ComName) %>% 
    summarize(mean_length_mm = mean(mean_length_mm)) %>% 
    inner_join(spp_names) %>% 
    mutate(mean_length_mm = ifelse(ComName == "Tidepool Sculpin", 89.0, mean_length_mm)) #no length in our df so taking the max length from fishbase
  
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
  
  fish_traits <- full_join(fork_length, milieu) %>% 
    select(3,1,2,4,5,6) %>% 
    left_join(feeding_guild) %>% 
    mutate(ComName = replace(ComName, ComName == "Pacific Sandfish", "Pacific sandfish")) %>%
    arrange(ComName)
  
  # write_csv(fish_traits, "data/fish_traits.csv")
  
  fish_trait_mat <- fish_traits %>% 
    select(-c(Species, ComName)) %>% 
    mutate_if(is.character, as.factor) %>% 
    clean_names() %>% 
    as.data.frame()
  
  rownames(fish_trait_mat) <- colnames(fish_L_mat)
  
  fish_trait_mat.t <- fish_trait_mat %>% mutate_if(is.numeric, log)
  
  return(list("trait" = fish_trait_mat.t, "abund" = fish_L_mat)) #no transformation on the abundance matrix
  
}

fish.list <- create_fish_matrices(net_tidy)

save(fish.list, file = here("data", "fish.list.Rdata"))