library(here)
library(tidyverse)
library(picante)
library(FD)
library(mFD)

#load tidy fish data frame created in 02_tidy_data
load(here("data", "net_tidy.Rdata")) 
load(here("data", "fish.list.Rdata")) #object created in 03_create_matrices

#picante
null_L <- list()
n_iter <- 999

for (i in 1:n_iter) {
  set.seed(i)
  null_L[[i]] <- randomizeMatrix(fish.list$abund, null.model ="independentswap", iterations = 1000)   
}

# with the FD package
FD_null_results <- list()

for (i in 1:n_iter){
  
  null_fishFD <- dbFD(x = fish.list$trait.t, #must be a df where character columns are factors or a distance matrix
                 a = null_L[[i]],
                 corr = "none", 
                 m = 5,
                 calc.FDiv = TRUE, 
                 print.pco = FALSE)
  
  null_FD_values <- cbind(null_fishFD$nbsp, null_fishFD$FRic, null_fishFD$FEve, null_fishFD$FDiv,
                          null_fishFD$FDis, null_fishFD$RaoQ) #extract indices
  
  FD_null_df <- null_FD_values %>%
    as_tibble(rownames = "sample")
  
  FD_null_results <- rbind(FD_null_results, FD_null_df)
}
  
colnames(FD_null_results)[2:7] <- c("Species_Richness", "FRic", "FEve", "FDiv", "FDis", "Rao")
  
FD_null_results <- FD_null_results %>% 
    separate_wider_delim(sample, delim = "_", names = c("site", "ipa"), cols_remove = FALSE) %>% 
    relocate(sample) %>% 
    mutate(site = factor(site, levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG"))) %>% 
    mutate(region = ifelse(site %in% c("FAM", "TUR", "COR"), "North", "South"), .after = site)

save(FD_null_results, file = "data/null_df.Rda")
load("data/null_df.Rda")  

FD_null_summary <- FD_null_results %>% 
  group_by(site) %>% 
  summarize(across(where(is.numeric), list(mean = mean, sd = sd)))

## use the equation in Zhang et al. to calculate SES
## OR figure out how to calculate significance with p-values
 