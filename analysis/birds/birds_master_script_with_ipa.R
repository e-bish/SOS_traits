
### Load libraries ###
## tidy data
library(tidyverse)
library(here)
library(readxl)
library(auk) #for bird species names
library(moments)

##FD analysis
library(ggrepel)
library(FD)
library(PNWColors)
library(patchwork)
library(vegan)

##ISA
library(indicspecies)

#### Create matrices ####
#load tidy fish data frame created in 01_tidy_data
load(here("data", "birds_tidy.Rdata")) 

#load trait data from Tobias et al. 2021 (Ecology Letters) AVONET: morphological, ecological and geographical data for all birds.
AVONET <- here("data", "AVONET Supplementary dataset 1.xlsx") %>% 
  read_excel(sheet = "AVONET2_eBird")

expand_species <- birds_tidy %>% 
  expand(nesting(year, month, site, ipa), spp_code) %>% 
  filter(!is.na(spp_code))

#create abundance matrix
bird_L <- birds_tidy %>% 
  filter(!is.na(comm_name)) %>% 
  full_join(expand_species) %>% 
  mutate(spp_sum = replace_na(spp_sum, 0)) %>% 
  group_by(site, ipa, year, spp_code) %>% 
  summarize(sp_avg = mean(spp_sum)) %>% #average across months
  ungroup() %>% 
  arrange(spp_code) %>% 
  pivot_wider(names_from = spp_code, values_from = sp_avg, values_fill = 0) %>% 
  mutate(sample = paste(site, ipa, year, sep = "_"), .after = year) %>% 
  select(!1:3) %>% 
  column_to_rownames(var = "sample") %>% 
  filter(rowSums(across(where(is.numeric)))!=0) #cornet natural 2021 didn't see any species

#load species ids for birds observed in the field
spp_id <- read_csv("data/bird_spp_info.csv", col_names = TRUE) %>% 
  unite(Species2, c("genus", "species"), sep = " ") %>% 
  filter(!(comm_name == "unk_bird" |
             comm_name == "unk_terrestrial_bird" |
             comm_name == "unk_duck" | 
             comm_name == "Harbor Seal")) %>%
  mutate(spp_code = case_when(comm_name == "Common/Caspian" ~ "CATE", 
                              comm_name == "Sparrow" ~ 'HOSP',
                              comm_name == "grebe" ~ 'WEGR', 
                              comm_name == "Hummingbird" ~ 'RUHU',
                              comm_name == "Cormorant" ~ 'PECO',
                              comm_name == "swallow" ~ "BARS", 
                              comm_name == "Glaucus winged gull" ~ 'GWGU',
                              comm_name == "gull" ~ 'HERG',
                              TRUE ~ spp_code)) %>% 
  mutate(Species2 = case_when(comm_name == "Common/Caspian" ~ 'Hydroprogne caspia', 
                              comm_name == "Double-crested cormorant" ~ 'Nannopterum auritum', #updated species name in 2014
                              comm_name == "Sparrow" ~ 'Passer domesticus',
                              comm_name == "grebe" ~ 'Aechmophorus occidentalis', 
                              comm_name == "Hummingbird" ~ 'Selasphorus rufus',
                              comm_name == "Cormorant" ~ 'Urile pelagicus',
                              comm_name == "swallow" ~ 'Hirundo rustica', 
                              comm_name == "Glaucus winged gull" ~ 'Larus glaucescens',
                              comm_name == "gull" ~ 'Larus argentatus',
                              TRUE ~ Species2)) %>% 
  mutate(spp_no = seq(1:nrow(.)))

### create the trait matrix
bird_traits <- AVONET %>% 
  select(Species2, 
         Mass,
         Primary.Lifestyle, 
         Trophic.Level, 
         Migration,
         Trophic.Niche) %>% 
  mutate(Trophic.Niche = str_replace(Trophic.Niche, " ", "_")) %>% 
  right_join(spp_id) %>% 
  mutate_if(is.character, as.factor) %>% 
  select(-c(Species2, comm_name, spp_no)) %>% 
  select(spp_code, everything()) %>% 
  arrange(spp_code) %>% 
  column_to_rownames(var="spp_code")

# write.csv(bird_traits, "data/bird_traits.csv", row.names = TRUE)

#inspect continuous traits 
# skewness(bird_traits$Mass)
# range(bird_traits$Mass)

#apply data transformation to continuous traits

bird_traits.t <- bird_traits %>% 
  mutate_if(is.numeric, log) %>% 
  as.data.frame()

# skewness(bird_traits.t$Mass)
# range(bird_traits.t$Mass)

bird.list <- list("trait" = bird_traits.t, "abund" = bird_L) 
# save(bird.list, file = here("data", "bird.list.Rdata"))

#### FD Analysis ####
# with the FD package 
birdFD <- dbFD(x = bird.list$trait, #must be a df where character columns are factors
               a = bird.list$abund,
               ord = "podani",
               corr = "none", 
               m = 5,
               calc.FDiv = TRUE, 
               print.pco = FALSE)

FD_values <- cbind(birdFD$nbsp, birdFD$FRic, birdFD$FEve, birdFD$FDiv,
                   birdFD$FDis) #extract indices

colnames(FD_values)[1:5] <- c("Species_Richness", "FRic", "FEve", "FDiv", "FDis")

FD_results <- FD_values %>% 
  as_tibble(rownames = "sample") %>% 
  separate_wider_delim(sample, delim = "_", names = c("site", "ipa", "year"), cols_remove = FALSE) %>% 
  relocate(sample) %>% 
  mutate(site = factor(site, levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG"))) %>% 
  mutate(region = ifelse(site %in% c("FAM", "TUR", "COR"), "North", "South"), .after = site) %>% 
  replace(is.na(.), 0) %>% 
  mutate(ipa = ifelse(ipa == "Natural2", "Natural", ipa))

plot_site_index <- function (index){
  ggplot(data = FD_results, aes(x = site, 
                                y = .data[[index]], 
                                color = site)) +
    geom_boxplot() +
    geom_point(size = 3) + 
    theme_classic() +
    ylab(index) +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
}

index_plots <- lapply(names(FD_results[6:10]), plot_site_index)

index_plots[[1]] + index_plots[[5]] + guide_area()+ index_plots[[2]] + index_plots[[3]] + index_plots[[4]] +  
  plot_layout(ncol = 3) + plot_layout(guides = "collect")

adonis2(FD_results[,c("Species_Richness","FDis", "FEve", "FRic", "FDiv")] ~ site, 
        strata = FD_results$year, data = FD_results, method = "euc")
#differences between sites

adonis2(FD_results$Species_Richness ~ site, strata = FD_results$year, data = FD_results, method = "euc")
adonis2(FD_results$FEve ~ site, strata = FD_results$year, data = FD_results, method = "euc")
adonis2(FD_results$FRic ~ site, strata = FD_results$year, data = FD_results, method = "euc")
adonis2(FD_results$FDiv ~ site, strata = FD_results$year, data = FD_results, method = "euc")
adonis2(FD_results$FDis ~ site, strata = FD_results$year, data = FD_results, method = "euc")
#pattern driven by species richness

plot_index <- function (index, by){
  ggplot(data = FD_results, aes(x = .data[[by]], 
                                y = .data[[index]], 
                                fill = .data[[by]])) +
    geom_boxplot() + 
    theme_classic() +
    ylab(index) +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
}

index_plots <- lapply(names(FD_results[6:10]), plot_index, by = "ipa")

index_plots[[1]] + index_plots[[5]] + guide_area() + index_plots[[2]] + index_plots[[3]] + index_plots[[4]] + 
  plot_layout(ncol = 3) + plot_layout(guides = "collect")

adonis2(FD_results[,c("Species_Richness","FDis", "FEve", "FRic", "FDiv")] ~ ipa, 
        strata = FD_results$site, data = FD_results, method = "euc")
#no differences between ipas

#### alpha diversity ####
bird_L_long <- bird.list$abund %>% 
  rownames_to_column("sample") %>% 
  pivot_longer(!sample, names_to = "species", values_to = "avg_n")

alpha_div <- bird_L_long %>% 
  group_by(sample) %>% 
  summarize(richness = specnumber(avg_n),
            shannon = diversity(avg_n, index = "shannon"),
            simpson = diversity(avg_n, index = "simpson"),
            invsimpson = diversity(avg_n, index = "invsimpson"),
            n = sum(avg_n)) %>% 
  ungroup() %>% 
  separate_wider_delim(sample, delim = "_", names = c("site", "ipa", "year"), cols_remove = FALSE) %>% #######
mutate(site = factor(site, levels = SOS_core_sites),
       veg = ifelse(site %in% c("TUR", "COR", "SHR"), "present", "absent"))

adonis2(alpha_div[,c("richness","shannon", "simpson", "invsimpson")] ~ site, 
        strata = alpha_div$year, data = alpha_div, method = "euc")
#significant differences bewteen sites

adonis2(alpha_div$richness ~ site, strata = alpha_div$year, data = alpha_div, method = "euc")
adonis2(alpha_div$shannon ~ site, strata = alpha_div$year, data = alpha_div, method = "euc")
adonis2(alpha_div$simpson ~ site, strata = alpha_div$year, data = alpha_div, method = "euc")
adonis2(alpha_div$invsimpson ~ site, strata = alpha_div$year, data = alpha_div, method = "euc")
#differences driven by species richness

#### ISA ####
site.ISA <- multipatt(x = bird_L, cluster = FD_results$site, duleg = TRUE)
summary(site.ISA)
#crows, gulls, and terns associated with SHR
#purple martins and great blue herons at DOK
#pigeon guillemots at EDG

#### compute null model ####
#load tidy fish data frame created in 02_tidy_data
# load(here("data", "net_tidy.Rdata"))
# load(here("data", "fish.list.Rdata")) #object created in 03_create_matrices

null_L <- list(bird.list$abund)
n_iter <- 1000 #observed data + 999 permutations

# #frequency is C4 null model in Gotzenberger et al. that works well for detection of environmental filtering
# #independentswap seems to be the method that other FD papers have used (e.g., Zhang) and is recommended by Swenson

for (i in 2:n_iter) {
  set.seed(i)
  null_L[[i]] <- randomizeMatrix(bird.list$abund, null.model ="independentswap", iterations = 1000)
}

# # calculate FD indices with the FD package
FD_null_output <- list()
for (i in 1:n_iter){
  
  null_birdFD <- dbFD(x = bird.list$trait, #must be a df where character columns are factors or a distance matrix
                      a = null_L[[i]],
                      corr = "cailliez",
                      m = 5, #might only keep 2 anyway
                      calc.FDiv = TRUE,
                      print.pco = FALSE)
  
  null_FD_values <- cbind(null_birdFD$nbsp, null_birdFD$FRic, null_birdFD$FEve, null_birdFD$FDiv,
                          null_birdFD$FDis) #extract indices
  
  FD_null_df <- null_FD_values %>%
    as_tibble(rownames = "sample")
  
  FD_null_output <- rbind(FD_null_output, FD_null_df)
}

colnames(FD_null_output)[2:6] <- c("Species_Richness", "FRic", "FEve", "FDiv", "FDis") 

FD_null_results <- FD_null_output %>%
  separate_wider_delim(sample, delim = "_", names = c("site", "ipa", "year"), cols_remove =  TRUE) %>%
  mutate(site = factor(site, levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG")),
         year = factor(year, levels = c("2018", "2019", "2021", "2022"))) %>%
  mutate(region = ifelse(site %in% c("FAM", "TUR", "COR"), "North", "South"), .after = site) %>% 
  replace(is.na(.), 0) 

# save(FD_null_results, file = "data/FD_null_results.Rda")

## load null dataset 
# load("data/FD_null_results.Rda")  

FD_null_summary <- FD_null_results %>% 
  group_by(site) %>% 
  summarize(across(where(is.numeric), list(mean = mean, sd = sd)))


#### calculate SES samples ####
# load("data/FD_results.Rda")  

FD_means <- FD_results %>% 
  group_by(site) %>% 
  summarize(across(where(is.numeric), mean)) %>% 
  ungroup()

SES_tab <- as.data.frame(FD_means$site)
SES_tab[,2] <- FD_means$Species_Richness 
SES_tab[,3] <- (FD_means$FRic - FD_null_summary$FRic_mean) / FD_null_summary$FRic_sd
SES_tab[,4] <- (FD_means$FEve - FD_null_summary$FEve_mean) / FD_null_summary$FEve_sd
SES_tab[,5] <- (FD_means$FDiv - FD_null_summary$FDiv_mean) / FD_null_summary$FDiv_sd
SES_tab[,6] <- (FD_means$FDis - FD_null_summary$FDis_mean) / FD_null_summary$FDis_sd

names(SES_tab) <- names(FD_means)

# calculate p values for SES
#cant get this to work in a function!!! had to create a different function for each metric
pull_FRic_ntiles <- function(site_ID) {
  #commented out code left here for checking that the proper order has been extracted
  lower <- FD_null_results %>% 
    filter(site == site_ID) %>% 
    # rownames_to_column(var = "index") %>% #
    arrange(FRic) %>% 
    # rownames_to_column(var = "arranged") %>% #
    slice(250) %>% #5th percentile of 5000 observations
    select(FRic)
  # select(index, arranged, FRic) #
  
  upper <- FD_null_results %>%
    filter(site == site_ID) %>%
    # rownames_to_column(var = "index") %>% #
    arrange(FRic) %>% 
    # rownames_to_column(var = "arranged") %>% #
    slice(4750) %>% #95th percentile of 5000 observations
    select(FRic)
  # select(index, arranged, FRic) #
  
  names <- c("site", 
             # "index_lower",
             # "arranged_lower",
             "FRic_lower",
             # "index_upper",
             # "arranged_upper",
             "FRic_upper")
  
  df <- data.frame(site_ID)
  df <- cbind(df, lower, upper)
  names(df) <- names
  
  return(df)
}

pull_FEve_ntiles <- function(site_ID) {
  
  lower <- FD_null_results %>% 
    filter(site == site_ID) %>% 
    arrange(FEve) %>% 
    slice(250) %>%
    select(FEve)
  
  upper <- FD_null_results %>%
    filter(site == site_ID) %>%
    arrange(FEve) %>% 
    ungroup() %>% 
    slice(4750) %>%
    select(FEve)
  
  names <- c("site", "FEve_lower","FEve_upper")
  
  df <- data.frame(site_ID)
  df <- cbind(df, lower, upper)
  names(df) <- names
  
  return(df)
}

pull_FDiv_ntiles <- function(site_ID) {
  
  lower <- FD_null_results %>% 
    filter(site == site_ID) %>% 
    arrange(FDiv) %>% 
    slice(250) %>%
    select( FDiv)
  
  upper <- FD_null_results %>%
    filter(site == site_ID) %>%
    arrange(FDiv) %>% 
    ungroup() %>% 
    slice(4750) %>%
    select(FDiv)
  
  names <- c("site", "FDiv_lower","FDiv_upper")
  
  df <- data.frame(site_ID)
  df <- cbind(df, lower, upper)
  names(df) <- names
  
  return(df)
}

pull_FDis_ntiles <- function(site_ID) {
  
  lower <- FD_null_results %>% 
    filter(site == site_ID) %>% 
    arrange(FDis) %>% 
    slice(250) %>%
    select(FDis)
  
  upper <- FD_null_results %>%
    filter(site == site_ID) %>%
    arrange(FDis) %>% 
    ungroup() %>% 
    slice(4750) %>%
    select(FDis)
  
  names <- c("site",  "FDis_lower","FDis_upper")
  
  df <- data.frame(site_ID)
  df <- cbind(df, lower, upper)
  names(df) <- names
  
  return(df)
}

SOS_core_sites <- factor(c("FAM", "TUR", "COR", "SHR", "DOK", "EDG"), 
                         levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG"))

#extract the 95% confidence intervals for each metric
FRic_CI <- lapply(SOS_core_sites, pull_FRic_ntiles)
FEve_CI <- lapply(SOS_core_sites, pull_FEve_ntiles)
FDiv_CI <- lapply(SOS_core_sites, pull_FDiv_ntiles)
FDis_CI <- lapply(SOS_core_sites, pull_FDis_ntiles)

#combine the lists using mapply
combine_CI <- t(mapply(c, FRic_CI, FEve_CI, FDiv_CI, FDis_CI)) %>% as.data.frame()

#turn the confidence intervals into long format
CI_uppers_lowers <- combine_CI %>% 
  select(!starts_with("site")) %>% 
  mutate(site = SOS_core_sites) %>% 
  pivot_longer(!site, names_to = "metric") %>% 
  separate_wider_delim(metric, delim = "_", names = c("metric", "range"), cols_remove = TRUE) %>% 
  mutate(value = unlist(value))

#turn the mean values into long format
FD_means_long <- FD_means %>% 
  pivot_longer(!c(site, Species_Richness), names_to = "metric", values_to = "value") %>% 
  select(!Species_Richness) %>% 
  mutate(range = "mean", .before = value)

#combine the intervals with the means into one df and check significance
df_for_p_vals <- bind_rows(CI_uppers_lowers,FD_means_long) %>% 
  pivot_wider(names_from = range, values_from = value) %>% 
  mutate(significant = ifelse(mean > lower & mean < upper, "no", "yes")) 

#format significance checks into a table that's the same format as the SES table
p_vals_tbl <- df_for_p_vals %>% 
  select(site, metric, significant) %>% 
  pivot_wider(names_from = metric, values_from = significant)
