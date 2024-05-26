library(tidyverse)
library(here)
library(FD)

load(here("data", "fish.list.Rdata")) #object created in 02_tidy_data

fish.traits <- fish.list$trait %>% 
  as.data.frame() %>% 
  mutate(mean_length_mm = as.numeric(mean_length_mm)) %>% 
  mutate_if(is.character, as.factor)

#### inspect continuous traits ####

par(mfrow=c(1,2))
hist(fish.traits$mean_length_mm, main="Histogram", xlab="Mean Length (mm)")
boxplot(fish.traits$mean_length_mm, main="Boxplot", ylab="Mean Length (mm)")
#trait is well distributed with no outliers

#check the coefficient of variation
CV <- function(x) { 100 * sd(x) / mean(x) }

CV(fish.traits$mean_length_mm)
# <50 is small, so there does not appear to be evidence that we need to relativize the length data

#### inspect categorical traits ####
par(mfrow=c(1,1))

barplot(table(fish.traits$body_shape_i, useNA="ifany"))
#body shape is well distributed with no missing data

barplot(table(fish.traits$demers_pelag, useNA="ifany"))
table(fish.traits$demers_pelag, useNA="ifany")
#only one pelagic species and 2 pelagic-neritic species
#ended up deciding to group the pelagic categories

barplot(table(fish.traits$migrations, useNA="ifany"))
table(fish.traits$migrations, useNA="ifany")
#some categories should be grouped

barplot(table(fish.traits$schooling, useNA="ifany"))
#well distributed

barplot(table(fish.traits$feeding_guild, useNA="ifany"))
table(fish.traits$feeding_guild, useNA="ifany")
#distributed well enough