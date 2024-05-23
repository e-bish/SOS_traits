library(tidyverse)
library(here)
library(janitor)
library(vegan)

load(here("data", "fish.list.Rdata")) #object created in 02_tidy_data

### Conduct data quality and integrity check using steps outlined in Palacio et al. 2022 
#### Step 1: Plot the community data matrix to assess the prevalence of zeros ####

fish.list$abund %>% 
  rownames_to_column(var = "sample") %>% 
  pivot_longer(cols = 2:43, names_to = "species") %>% 
  ggplot() +
  geom_tile(aes(x = species, y = sample, fill = value))

#lots of zeros

#### Step 2: Check species sampling coverage (e.g. rarefaction) ####
rarecurve(floor(fish.list$abund))


rare_mat <- fish.list$abund %>% 
  rownames_to_column(var = "sample") %>% 
  pivot_longer(cols = 2:43, names_to = "species") %>% 
  separate_wider_delim(sample, delim = "_", names = c("year", "month", "site"), cols_remove = FALSE) %>% 
  group_by(site) %>% 
  arrange(year, factor(month, levels = c("Apr", "May", "Jun", "Jul", "Aug", "Sept"))) 

group_map(rare_mat, rarecurve)

#### Step 3: Plot the distribution of continuous traits to check for outliers and continuous traits to check the balance of levels ####

#### Step 4: Evaluate multicollinearity among continuous traits and associations with categorical traits ####

#### Step 5: Identify missing trait data ####








