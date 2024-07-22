library(here)
library(tidyverse)
library(picante)

set.seed(1993)

#load tidy fish data frame created in 02_tidy_data
load(here("data", "net_tidy.Rdata")) 
load(here("data", "fish.list.Rdata")) #object created in 03_create_matrices

null_mat <- randomizeMatrix(fish.list$abund, null.model ="independentswap", iterations = 1000)   
