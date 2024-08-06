library(tidyverse)
library(vegan)

load("data/fish.list.Rdata") #object created in 03_create_matrices
SOS_core_sites <- c("FAM", "TUR", "COR", "SHR", "DOK", "EDG")


fish_L <- fish.list$abund %>% 
  as_tibble(rownames = "sample") %>% 
  separate_wider_delim(sample, delim = "_", names = c("site", "ipa"), cols_remove = FALSE) 
  
  
#SIMPER code borrowed from Jon Bakker
simper.results <- simper(fish_L[4:45], fish_L$site)
simper.results
summary(simper.results)

comparisons <- names(simper.results)
simper.df <- c()

for(i in 1:length(comparisons)) {
  temp <- summary(simper.results)[as.character(comparisons[i])] %>%
    as.data.frame()
  colnames(temp) <- gsub(
    paste(comparisons[i],".", sep = ""), "", colnames(temp))
  temp <- temp %>%
    mutate(Comparison = comparisons[i],
           Position = row_number()) %>%
    rownames_to_column(var = "Species")
  simper.df <- rbind(simper.df, temp)
}


simper.df %>%
  filter(p <= 0.05) %>%
  select(Species, average, Comparison, Position) %>% 
  arrange(Species)

#### ISA ####
library(indicspecies)

test.ISA <- multipatt(x = fish_L[4:45], cluster = fish_L$site, duleg = TRUE)
summary(test.ISA)

#### TITAN ####
library(TITAN2)
test.titan <- titan(env = FD_results$FRic, txa = fish_L[4:45])
#wont run if a species occurs <3 times
