library(tidyverse)
library(here)
library(FD)

load(here("data", "fish.list.Rdata")) #object created in 02_tidy_data
here("analysis", "general_functions", "geb12299-sup-0002-si.r") %>% source()
here("analysis", "general_functions", "scree.r") %>% source()

fish.traits <- fish.list$trait %>% 
  as.data.frame() %>% 
  mutate(mean_length_mm = as.numeric(mean_length_mm)) %>% 
  mutate_if(is.character, as.factor)

#create species x species matrix
fish.gowdist <- gowdis(fish.traits)

#view pcoa
# pcoa <- cmdscale(fish.gowdist)
# colnames(pcoa) <- c("pcoa1", "pcoa2")
# 
# pcoa %>%
#   as_tibble(rownames = "species") %>%
#   ggplot(aes(x = pcoa1, y = pcoa2)) +
#   geom_text(
#     label=rownames(pcoa),
#     check_overlap = T
#   )

#### test the quality of the functional space ####
# space_qual <- qual_funct_space(mat_funct = fish.gowdist, nbdim = 6, metric = "Gower")
# space_qual$meanSD # very little difference between 4 and >4 dimensions (we're shooting for 0 here)
# #"The mSD is 0 when the functional space perfectly represents the initial distance and increases as pairs of species become less represented in the functional space (Maire et al., 2015)."
# 
# #### test to see how many groups to use ####
# test1 <- hclust(fish.gowdist, method = "average") #meh
# test2 <- hclust(fish.gowdist, method = "complete") #not terrible
# test3 <- hclust(fish.gowdist, method = "median") #terrible
# fish.hclust <- hclust(fish.gowdist, method = "ward.D2") #better!
# plot(fish.hclust, hang = -1,
#      main = "Method = Ward’s minimum variance",
#      xlab = "Species",
#      sub = "BC dissimilarity; ward.D2 linkage")
# scree(hclust.obj = fish.hclust) # not a solid answer here but 5 is generally at the bottom of the elbow, so we'll try 5 groups

#### FD indices calculation ####################################################


fishFD <- dbFD(fish.traits, #must be a df where character columns are factors
                fish.list[[2]], 
                 m = 4, #seems to keep 4 axis no matter how many you specify?
                 calc.FGR = TRUE, 
                 clust.type = "ward.D2",
                 calc.FDiv = TRUE, #won't return a value if there are only categorical traits
                 print.pco = TRUE)
#it will prompt to see if you want to separate by groups (g) and ask how many (5)
  
important.indices <- cbind(fishFD$nbsp, fishFD$FRic, fishFD$FEve, fishFD$FDiv,
                             fishFD$FDis, fishFD$RaoQ) #extract indices
  
colnames(important.indices) <- c("NumbSpecies", "FRic", "FDiv", "FEve", "FDis", "Rao")

FD_results <- important.indices %>% 
  as_tibble(rownames = "site_month") %>% 
  separate_wider_delim(site_month, delim = "_", names = c("site", "month"), cols_remove = FALSE) %>% 
  relocate(site_month) %>% 
  mutate(month = factor(month, labels = c("Apr", "May", "Jun", "Jul", "Aug", "Sept"))) %>% 
  mutate_if(is.character, as.factor)

plot_index <- function (index, index_name){
  FD_results %>%
    filter(!site %in% c("MA", "WA", "HO", "TL", "LL", "PR")) %>% 
    ggplot(aes(x = site, y = index, fill = site)) +
    geom_boxplot(outlier.shape = NA) + 
    theme_classic() +
    ylab(index_name) +
    theme(axis.title.x = element_blank(), 
          text = element_text(size = 18),
          #legend.position = "none",
          axis.text = element_text(size = 16))
}

index_plots <- list()

index_plots[1:5] <- FD_results %>% 
  select(NumbSpecies, FRic, FEve, FDis, Rao) %>% 
  map2(.y = names(.), ~ plot_index(.x, .y))

index_plots$FEve <- FD_results %>%
  filter(!is.na(FEve)) %>% 
  ggplot(aes(x = factor(site, levels = c("FAM", "TUR", "COR", "MA", "WA", "HO", "SHR", "DOK", "LL", "TL", "PR", "EDG")), y = FEve, fill = site)) +
  geom_boxplot(outlier.shape = NA) + 
  theme_classic() +
  ylab("Number of Species") +
  theme(axis.title.x = element_blank(), 
        text = element_text(size = 16),
        axis.text = element_text(size = 16))

# index_plots$FRic <- FD_results %>%
#   # filter(!is.na(FEve)) %>% 
#   ggplot(aes(x = shoreline, y = FRic, fill = shoreline)) +
#   geom_boxplot(outlier.shape = NA) + 
#   theme_classic() +
#   ylab("Number of Species") +
#   theme(axis.title.x = element_blank(), 
#         text = element_text(size = 16),
#         axis.text = element_text(size = 16))

index_plots[[1]]


#### RLQ analysis ###############################################################


            