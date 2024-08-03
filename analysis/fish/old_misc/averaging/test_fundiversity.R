library(here)
library(tidyverse)
library(fundiversity)
library(FD)

#load tidy fish data frame created in 02_tidy_data
load(here("data", "net_tidy.Rdata")) 
load(here("data", "fish.list.Rdata")) #object created in 03_create_matrices

distance_mat <- gowdis(fish.list$trait) #traits are standardized (given equal weight) in this step
all_vectors <- ape::pcoa(distance_mat) #no correction for negative eigenvalues, by default they are removed
trait_space <- all_vectors$vectors

n_axes_to_retain <- 5

FRic <- fd_fric(trait_space[,1:n_axes_to_retain], fish.list$abund, stand = TRUE)
