library(here)
library(tidyverse)
library(picante)
library(FD)

#### create null dataset ####
# code to create null dataframe with 999 permutations. This takes a long time so skip below to load the data


# #load tidy fish data frame created in 02_tidy_data
# load(here("data", "net_tidy.Rdata")) 
# load(here("data", "fish.list.Rdata")) #object created in 03_create_matrices
# 
# null_L <- list(fish.list$abund.ipa)
# n_iter <- 1000 #observed data + 999 permutations
# 
# #frequency is C4 null model in Gotzenberger et al. that works well for detection of environmental filtering
# #independentswap seems to be the method that other FD papers have used (e.g., Zhang) and is recommended by Swenson
# for (i in 2:n_iter) {
#   set.seed(i)
#   null_L[[i]] <- randomizeMatrix(fish.list$abund.ipa, null.model ="independentswap", iterations = 1000) 
# }
# 
# # calculate FD indices with the FD package
# FD_null_results <- list()
# for (i in 1:n_iter){
#   
#   null_fishFD <- dbFD(x = fish.list$trait.t, #must be a df where character columns are factors or a distance matrix
#                  a = null_L[[i]],
#                  corr = "cailliez", 
#                  m = 5,
#                  calc.FDiv = TRUE, 
#                  print.pco = FALSE)
#   
#   null_FD_values <- cbind(null_fishFD$nbsp, null_fishFD$FRic, null_fishFD$FEve, null_fishFD$FDiv,
#                           null_fishFD$FDis) #extract indices
#   
#   FD_null_df <- null_FD_values %>%
#     as_tibble(rownames = "sample")
#   
#   FD_null_results <- rbind(FD_null_results, FD_null_df)
# }
#   
# colnames(FD_null_results)[2:6] <- c("Species_Richness", "FRic", "FEve", "FDiv", "FDis")
#   
# FD_null_results <- FD_null_results %>% 
#     separate_wider_delim(sample, delim = "_", names = c("site", "ipa"), cols_remove = FALSE) %>% 
#     relocate(sample) %>% 
#     mutate(site = factor(site, levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG"))) %>% 
#     mutate(region = ifelse(site %in% c("FAM", "TUR", "COR"), "North", "South"), .after = site)

# save(FD_null_results, file = "data/FD_null_results.Rda")

#### load null dataset ####
load("data/FD_null_results.Rda")  

FD_null_summary <- FD_null_results %>% 
  group_by(site) %>% 
  summarize(across(where(is.numeric), list(mean = mean, sd = sd)))

#### calculate SES for bootstrapped samples ####
load("data/FD_boot_results.Rda")  

FD_boot_means <- FD_boot_results %>% 
  replace(is.na(.), 0) %>% #replace na with zero for EDG armored where there was only two species caught at one ipa so cant calculate FRic, FEve, FDiv
  group_by(site) %>% 
  summarize(across(where(is.numeric), mean))

SES_boot_tab <- as.data.frame(FD_boot_means$site)
SES_boot_tab[,2] <- FD_boot_means$Species_Richness 
SES_boot_tab[,3] <- (FD_boot_means$FRic - FD_null_summary$FRic_mean) / FD_null_summary$FRic_sd
SES_boot_tab[,4] <- (FD_boot_means$FEve - FD_null_summary$FEve_mean) / FD_null_summary$FEve_sd
SES_boot_tab[,5] <- (FD_boot_means$FDiv - FD_null_summary$FDiv_mean) / FD_null_summary$FDiv_sd
SES_boot_tab[,6] <- (FD_boot_means$FDis - FD_null_summary$FDis_mean) / FD_null_summary$FDis_sd

names(SES_boot_tab) <- names(FD_boot_means)

#### calculate significance with p-values ####
SOS_core_sites <- c("FAM", "TUR", "COR", "SHR", "DOK", "EDG")

#cant get this to work in a function!!! had to create a different function for each metric
pull_FRic_ntiles <- function(site_ID) {
  #commented out code left here for checking that the proper order has been extracted
  lower <- FD_null_results %>% 
    filter(site == site_ID) %>% 
    # rownames_to_column(var = "index") %>% 
    arrange(FRic) %>% 
    # rownames_to_column(var = "arranged") %>% 
    slice(25) %>%
    select(FRic)
    # select(index, arranged, FRic)
  
  upper <- FD_null_results %>%
    filter(site == site_ID) %>%
    # rownames_to_column(var = "index") %>% 
    arrange(FRic) %>% 
    # rownames_to_column(var = "arranged") %>% 
    ungroup() %>% 
    slice(975) %>%
    select(FRic)
    # select(index, arranged, FRic)
  
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
    slice(25) %>%
    select(FEve)
  
  upper <- FD_null_results %>%
    filter(site == site_ID) %>%
    arrange(FEve) %>% 
    ungroup() %>% 
    slice(975) %>%
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
    slice(25) %>%
    select( FDiv)
  
  upper <- FD_null_results %>%
    filter(site == site_ID) %>%
    arrange(FDiv) %>% 
    ungroup() %>% 
    slice(975) %>%
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
    slice(25) %>%
    select(FDis)
  
  upper <- FD_null_results %>%
    filter(site == site_ID) %>%
    arrange(FDis) %>% 
    ungroup() %>% 
    slice(975) %>%
    select(FDis)
  
  names <- c("site",  "FDis_lower","FDis_upper")
  
  df <- data.frame(site_ID)
  df <- cbind(df, lower, upper)
  names(df) <- names
  
  return(df)
}

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

#turn the mean boostrapped values into long format
FD_boot_means_long <- FD_boot_means %>% 
  pivot_longer(!site, names_to = "metric", values_to = "value") %>% 
  mutate(range = "boot_mean", .before = value)

#get the sd from the bootstrapped values
FD_boot_sd <- FD_boot_results %>% 
  replace(is.na(.), 0) %>% #replace na with zero for EDG armored where there was only two species caught at one ipa so cant calculate FRic, FEve, FDiv
  group_by(site) %>% 
  summarize(across(where(is.numeric), sd)) %>% 
  pivot_longer(!site, names_to = "metric", values_to = "sd") %>% 
  mutate(lower_sd = FD_boot_means_long$value - sd, upper_sd = FD_boot_means_long$value + sd)

FD_boot_sd2 <- FD_boot_sd %>% 
  filter(!metric == "Species_Richness")

#combine the intervals with the bootstrapped means into one df and check significance
df_for_p_vals <- bind_rows(CI_uppers_lowers,  filter(FD_boot_means_long, !metric == "Species_Richness")) %>% 
  pivot_wider(names_from = range, values_from = value) %>% 
  mutate(significant = ifelse(boot_mean > lower & 
                                boot_mean < upper, "no", "yes")) %>% 
  mutate(significant_lower = ifelse(significant == "yes" & 
                                 FD_boot_sd2$lower_sd < lower |
                                 FD_boot_sd2$lower_sd > upper, "yes", "no")) %>% 
  mutate(significant_upper = ifelse(significant == "yes" & 
                                      FD_boot_sd2$upper_sd < lower |
                                      FD_boot_sd2$upper_sd > upper, "yes", "no")) %>% 
  mutate(significance_total = ifelse(significant == "yes" & significant_lower == "yes" & significant_upper == "yes", "yes", "no"))

#format significance checks into a table that's the same format as the FD boot means table
p_vals_tbl <- df_for_p_vals %>% 
  select(site, metric, significance_total) %>% 
  pivot_wider(names_from = metric, values_from = significance_total)


#test one way anova
# library(wPerm)
# perm.oneway.anova(FD_results$Species_Richness, FD_null_results$Species_Richness)
