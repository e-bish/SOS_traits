library(picante)
library(rcompanion)

load(here("data", "fish.list.Rdata")) #object created in 02_tidy_data

apply(fish.list$abund, 1, sum) #rows
apply(fish.list$abund, 2, sum) #columns

confirm_proper_permutations <- replicate(5, randomizeMatrix(fish.list$abund, null.model = "independentswap"))

apply(confirm_proper_permutations, 3, rowSums) #site abundance varies from observed
apply(confirm_proper_permutations, 3, colSums) #species abundances are fixed

null_matrices <- replicate(999, randomizeMatrix(fish.list$abund, null.model = "independentswap"))

#### calculate the SES and p-values for each site

# #Calculate standardized effect sizes
# #SES = mean(obs) - mean(sim) / sd(sim)
# ses.all <- (fish.list$abund - apply(null_matrices, 1, mean))/apply(null_matrices, 1, sd) 
# #Negative SES values indicated trait convergence, while positive values indicated trait divergence (Götzenberger et al., 2012).
# 
# #Test significance from null expectation (using p-value)
# # p-value =quantile.obs/(total.interation+1)
# p.val.all <- apply(cbind(fish.list$abund, null_matrices), 1, rank)[1, ]/1000

