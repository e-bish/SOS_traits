library(ade4)
library(tidyverse)

here("analysis", "fish", "02_tidy_data.R") %>% source()

fish.list <- create_fish_matrices(net_tidy)

fish.traits <- fish.list$trait