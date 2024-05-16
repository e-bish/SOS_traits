library(tidyverse)
library(here)
library(ggrepel)
library(FD)
library(ggordiplots)
library(PNWColors)
library(patchwork)
library(vegan)
library(mFD)

load(here("data", "fish.list.Rdata")) #object created in 02_tidy_data
here("analysis", "general_functions", "geb12299-sup-0002-si.r") %>% source()
here("analysis", "general_functions", "scree.r") %>% source()

#create species x species matrix
fish.gowdist <- gowdis(fish.list$trait)

#view pcoa
pcoa <- cmdscale(sqrt(fish.gowdist)) 
colnames(pcoa) <- c("pcoa1", "pcoa2")

pcoa %>%
  as_tibble(rownames = "species") %>%
  ggplot(aes(x = pcoa1, y = pcoa2)) +
  geom_point(size = 3) +
  labs(x="PCoA Axis 1 (27.53%)",
       y="PCoA Axis 2 (17.37%)") +
  theme_bw(base_size = 16)

pcoa %>%
  as_tibble(rownames = "species") %>%
  ggplot(aes(x = pcoa1, y = pcoa2)) +
  geom_text_repel(
    label=rownames(pcoa),
    max.time = 1,
    max.overlaps = Inf) +
  geom_point() +
  labs(x="PCoA Axis 1 (27.53%)",
       y="PCoA Axis 2 (17.37%)") +
  theme_bw(base_size = 16)

#### FD indices calculation ####################################################
fishFD <- dbFD(x = fish.list$trait, #must be a df where character columns are factors
               a = fish.list$abund,
               corr = "sqrt", #this is the default option to deal with negative eigenvalues
               m = 4, #seems to keep 4 axis no matter how many you specify?
               calc.FGR = TRUE, 
               clust.type = "ward.D2",
               calc.FDiv = TRUE, #won't return a value if there are only categorical traits
               print.pco = TRUE)
#it will prompt to see if you want to separate by groups (g) and ask how many (7)

important.indices <- cbind(fishFD$nbsp, fishFD$FRic, fishFD$FEve, fishFD$FDiv,
                           fishFD$FDis, fishFD$RaoQ) #extract indices

colnames(important.indices) <- c("NumbSpecies", "FRic", "FDiv", "FEve", "FDis", "Rao")

FD_results <- important.indices %>% 
  as_tibble(rownames = "sample") %>% 
  separate_wider_delim(sample, delim = "_", names = c("year", "month", "site"), cols_remove = FALSE) %>% 
  relocate(sample) %>% 
  mutate(month = factor(month, labels = c("Apr", "May", "Jun", "Jul", "Aug", "Sept"))) %>% 
  mutate(site = factor(site, levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG"))) %>% 
  replace(is.na(.), 0) #can't calculate indices with <3 species

FD_results %>% 
  filter(NumbSpecies >2) %>%
  group_by(site) %>% 
  summarize(across(where(is.numeric), mean))

DOK_df <- fish.list[[2]] %>% 
  as_tibble(rownames = "sample") %>% 
  separate_wider_delim(sample, delim = "_", names = c("year", "month", "site"), cols_remove = FALSE) %>% 
  filter(site == "DOK") %>% 
  select(-c(1:4)) %>% 
  select_if(colSums(.) != 0)

DOK_spp <- colnames(DOK_df)

pcoa %>% 
  as_tibble(rownames = "species") %>%
  filter(species %in% DOK_spp) %>% 
  ggplot(aes(x = pcoa1, y = pcoa2)) +
    geom_text_repel(
      label=DOK_spp,
      max.time = 1,
      max.overlaps = Inf) +
    geom_point() +
    labs(x="PCoA Axis 1 (27.53%)",
       y="PCoA Axis 2 (17.37%)") +
    theme_bw(base_size = 16)

COR_df <- fish.list[[2]] %>% 
  as_tibble(rownames = "sample") %>% 
  separate_wider_delim(sample, delim = "_", names = c("year", "month", "site"), cols_remove = FALSE) %>% 
  filter(site == "COR") %>% 
  select(-c(1:4)) %>% 
  select_if(colSums(.) != 0)

COR_spp <- colnames(COR_df)

pcoa %>% 
  as_tibble(rownames = "species") %>%
  filter(species %in% COR_spp) %>% 
  ggplot(aes(x = pcoa1, y = pcoa2)) +
  geom_text_repel(
    label=COR_spp,
    max.time = 1,
    max.overlaps = Inf) +
  geom_point() +
  labs(x="PCoA Axis 1 (27.53%)",
       y="PCoA Axis 2 (17.37%)") +
  theme_bw(base_size = 16)


