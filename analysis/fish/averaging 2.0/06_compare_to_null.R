library(here)
library(tidyverse)
library(picante)
library(FD)

set.seed(1993)

#load tidy fish data frame created in 02_tidy_data
load(here("data", "net_tidy.Rdata")) 
load(here("data", "fish.list.Rdata")) #object created in 03_create_matrices

#picante - creates a single matrix
null_mat <- randomizeMatrix(fish.list$abund, null.model ="independentswap", iterations = 1000)   

#vegan
null_object <- nullmodel(fish.list$abund, method = "abuswap_r") #limited to algorithms that accept non-integer values
test <- simulate(null_object, nsim = 2)

#another vegan function
permatswap(fish.list$abund, "quasiswap")

#using the FD package
dist_mat <- gowdis(fish.list$trait.t) #"classic" method matches mFD, which treats categorical variables as continuous

n_axes_to_retain

# with the FD package 
null_fishFD <- dbFD(x = trait_space, #must be a df where character columns are factors or a distance matrix
               a = null_mat,
               # stand.x = FALSE, #standardization by range is automatic for Gower's distances
               corr = "none", 
               m = n_axes_to_retain,
               calc.FDiv = TRUE, 
               scale.RaoQ = TRUE, #scale Rao's Q 0-1 to make comparable
               print.pco = FALSE)

null_FD_values <- cbind(null_fishFD$nbsp, null_fishFD$FRic, null_fishFD$FEve, null_fishFD$FDiv,
                        null_fishFD$FDis, null_fishFD$RaoQ) #extract indices

colnames(null_FD_values) <- c("Species_Richness", "FRic", "FEve", "FDiv", "FDis", "Rao")

nulll_FD_results <- null_FD_values %>% 
  as_tibble(rownames = "sample") %>% 
  separate_wider_delim(sample, delim = "_", names = c("site", "ipa"), cols_remove = FALSE) %>% 
  relocate(sample) %>% 
  mutate(site = factor(site, levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG"))) %>% 
  mutate(region = ifelse(site %in% c("FAM", "TUR", "COR"), "North", "South"), .after = site)
