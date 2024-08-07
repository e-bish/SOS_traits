library(here)
library(tidyverse)
library(picante)
library(FD)
library(mFD)

#load tidy fish data frame created in 01_tidy_data
load(here("data", "birds_tidy.Rdata")) 
load(here("data", "bird.list.Rdata")) #object created in 02_create_matrices

null_L <- list()
n_iter <- 999

#randomize the abundance matrix
for (i in 1:n_iter) {
  set.seed(i)
  null_L[[i]] <- randomizeMatrix(bird.list$abund, null.model ="independentswap", iterations = 1000) 
}

# calculate FD indices with the FD package
FD_null_results <- list()

for (i in 1:n_iter){
  
  null_birdFD <- dbFD(x = bird.list$trait, #must be a df where character columns are factors or a distance matrix
                      a = null_L[[i]],
                      ord = "podani",
                      corr = "none", 
                      m = 5,
                      calc.FDiv = TRUE, 
                      print.pco = FALSE)
  
  null_FD_values <- cbind(null_birdFD$nbsp, null_birdFD$FRic, null_birdFD$FEve, null_birdFD$FDiv,
                          null_birdFD$FDis) #extract indices
  
  FD_null_df <- null_FD_values %>%
    as_tibble(rownames = "sample")
  
  FD_null_results <- rbind(FD_null_results, FD_null_df)
}

colnames(FD_null_results)[2:6] <- c("Species_Richness", "FRic", "FEve", "FDiv", "FDis")

FD_null_results <- FD_null_results %>% 
  separate_wider_delim(sample, delim = "_", names = c("site", "ipa"), cols_remove = FALSE) %>% 
  relocate(sample) %>% 
  mutate(site = factor(site, levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG"))) %>% 
  mutate(region = ifelse(site %in% c("FAM", "TUR", "COR"), "North", "South"), .after = site)

FD_null_bird_results <- FD_null_results
save(FD_null_results, file = "data/FD_null_bird_results.Rda")
load("data/FD_null_bird_results.Rda")  

FD_null_summary <- FD_null_results %>% 
  group_by(site) %>% 
  summarize(across(where(is.numeric), list(mean = mean, sd = sd))) 

## use the equation in Zhang et al. to calculate SES in original dataset
load("data/FD_bird_results.Rda")  

FD_means <- FD_bird_results %>% 
  group_by(site) %>% 
  summarize(across(where(is.numeric), mean))

SES_tab <- as.data.frame(FD_means$site)
SES_tab[,2] <- (FD_means$Species_Richness - FD_null_summary$Species_Richness_mean) / FD_null_summary$Species_Richness_sd
SES_tab[,3] <- (FD_means$FRic - FD_null_summary$FRic_mean) / FD_null_summary$FRic_sd
SES_tab[,4] <- (FD_means$FEve - FD_null_summary$FEve_mean) / FD_null_summary$FEve_sd
SES_tab[,5] <- (FD_means$FDiv - FD_null_summary$FDiv_mean) / FD_null_summary$FDiv_sd
SES_tab[,6] <- (FD_means$FDis - FD_null_summary$FDis_mean) / FD_null_summary$FDis_sd

names(SES_tab) <- names(FD_means[1:6])

## calculate SES for bootstrapped samples 
load("data/FD_bird_boot_results.Rda")  

FD_boot_means <- FD_bird_boot_results[1:6] %>% 
  group_by(site) %>% 
  summarize(across(where(is.numeric), mean))

SES_boot_tab <- as.data.frame(FD_boot_means$site)
names(SES_boot_tab) <- "site"
SES_boot_tab[,2] <- (FD_boot_means$Species_Richness - FD_null_summary$Species_Richness_mean) / FD_null_summary$Species_Richness_sd
SES_boot_tab[,3] <- (FD_boot_means$FRic - FD_null_summary$FRic_mean) / FD_null_summary$FRic_sd
SES_boot_tab[,4] <- (FD_boot_means$FEve - FD_null_summary$FEve_mean) / FD_null_summary$FEve_sd
SES_boot_tab[,5] <- (FD_boot_means$FDiv - FD_null_summary$FDiv_mean) / FD_null_summary$FDiv_sd
SES_boot_tab[,6] <- (FD_boot_means$FDis - FD_null_summary$FDis_mean) / FD_null_summary$FDis_sd

names(SES_boot_tab) <- names(FD_boot_means)

SES_all <- full_join(SES_tab, SES_boot_tab, by = "site", suffix = c(".o", ".boot")) #.o for original
#Negative SES values indicated trait convergence, while positive values indicated trait divergence (Götzenberger et al., 2012).

## calculate significance with p-values

SOS_core_sites <- c("FAM", "TUR", "COR", "SHR", "DOK", "EDG")

ggplot() + 
  geom_histogram(data = FD_bird_boot_results, mapping = aes(Species_Richness)) +
  facet_wrap(~site)

calc_p_values <- function(HA_df) {
  
  p_tab <- as.data.frame(matrix(nrow = 6, ncol = 6))
  names(p_tab) <- names(FD_boot_means)
  p_tab[,1] <- FD_boot_means$site
  
  for(i in 1:length(SOS_core_sites)) {
    
    tmp_null <- FD_null_results %>% filter(site == SOS_core_sites[i])
    tmp_HA <- HA_df %>% filter(site == SOS_core_sites[i])
    
    p_tab[i,2] <- sum(tmp_null$Species_Richness >= tmp_HA$Species_Richness) / nrow(tmp_null)
    p_tab[i,3] <- sum(tmp_null$FRic >= tmp_HA$FRic) / nrow(tmp_null)
    p_tab[i,4] <- sum(tmp_null$FEve >= tmp_HA$FEve) / nrow(tmp_null)
    p_tab[i,5] <- sum(tmp_null$FDiv >= tmp_HA$FDiv) / nrow(tmp_null)
    p_tab[i,6] <- sum(tmp_null$FDis >= tmp_HA$FDis) / nrow(tmp_null)
    
  }
  
  return(p_tab)
}

P_vals <- calc_p_values(FD_means)
P_boot_vals <- calc_p_values(FD_boot_means)
