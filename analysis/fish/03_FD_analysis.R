library(tidyverse)
library(here)
library(ggrepel)
library(FD)
library(PNWColors)
library(patchwork)
library(vegan)

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
pcoa <- cmdscale(fish.gowdist)
colnames(pcoa) <- c("pcoa1", "pcoa2")

fish_pcoa <- pcoa %>%
  as_tibble(rownames = "species") %>%
  ggplot(aes(x = pcoa1, y = pcoa2)) +
  geom_text_repel(
    label=rownames(pcoa),
    max.time = 1,
    max.overlaps = Inf) +
  geom_point(color = "darkblue") +
  theme_bw()

ggsave(plot = fish_pcoa, "docs/figures/fish_pcoa.png")

#### test the quality of the functional space ####
space_qual <- qual_funct_space(mat_funct = fish.gowdist, nbdim = 6, metric = "Gower")
space_qual$meanSD # very little difference between 4 and >4 dimensions (we're shooting for 0 here)
#"The mSD is 0 when the functional space perfectly represents the initial distance and increases as pairs of species become less represented in the functional space (Maire et al., 2015)."

#### test to see how many groups to use ####
test1 <- hclust(fish.gowdist, method = "average") #meh
test2 <- hclust(fish.gowdist, method = "complete") #not terrible
test3 <- hclust(fish.gowdist, method = "median") #terrible
fish.hclust <- hclust(fish.gowdist, method = "ward.D2") #better!
plot(fish.hclust, hang = -1,
     main = "Method = Ward’s minimum variance",
     xlab = "Species",
     sub = "BC dissimilarity; ward.D2 linkage")
scree(hclust.obj = fish.hclust) # not a solid answer here but 7 is generally at the bottom of the elbow, so we'll try 5 groups

#### FD indices calculation ####################################################
fishFD <- dbFD(fish.traits, #must be a df where character columns are factors
                fish.list[[2]], 
                 m = 4, #seems to keep 4 axis no matter how many you specify?
                 calc.FGR = TRUE, 
                 clust.type = "ward.D2",
                 calc.FDiv = TRUE, #won't return a value if there are only categorical traits
                 print.pco = TRUE)
#it will prompt to see if you want to separate by groups (g) and ask how many (7)

################################################################################
## More advanced PCoA
library(ggordiplots)

plot(pcoa)

spp_PCoA <- fishFD$x.axes #extract coordinates from PCoA
spp_groups <- fishFD$spfgr %>% as.factor() #extract groups

spp_groups_df <- as_tibble(spp_groups) %>% 
  mutate(species = labels(spp_groups))

#variance per axis
var_axis <- apply(spp_PCoA, 2, var)
## Proportional variance
var_axis/sum(var_axis)

# create basic ordiplot object to pull into advanced plot later
PCoA_ordiplot <- gg_ordiplot(ord = spp_PCoA,
                             groups = spp_groups,
                             ellipse = TRUE,
                             hull = FALSE,
                             spiders = FALSE)

#extract groups
PCoA_ellipses <- PCoA_ordiplot$df_ellipse 
# PCoA_spiders <- PCoA_ordiplot$df_spiders %>% 
#   rename(A1 = x, A2 = y) %>%
#   as_tibble() 

# basic plot building block
ord_gg <- ggplot() + theme_bw(16) + 
  geom_vline(xintercept = 0, lty=3, color="darkgrey") +
  geom_hline(yintercept = 0, lty=3, color="darkgrey") +
  theme(panel.grid=element_blank())

# pal <- pnw_palette("Bay", type = "discrete")

# advanced plot with spiders, sans vectors
fish_pcoa2 <- ord_gg +
  # geom_text_repel(data=spp_PCoA,
  #           aes(x=A1, y=A2,
  #               colour=spp_groups,
  #               label= rownames(pcoa)),
  #           #label.padding=unit(0.3,"lines"),
  #           fontface="bold")  +
  labs(x="PCoA Axis 1 (29.88%)",
       y="PCoA Axis 2 (15.46%)") +
  geom_point(data=spp_PCoA,
             aes(x=A1, y=A2,
                 colour=spp_groups))  +
  geom_polygon(data = PCoA_ellipses, 
               aes(x=x, y=y,
                   group = Group,
                   color = Group),
               fill = NA,
               linetype = "dashed",
               show.legend=FALSE);fish_pcoa2

ggsave(plot = fish_pcoa2, "docs/figures/fish_pcoa2.png")

###############################################################################
  
important.indices <- cbind(fishFD$nbsp, fishFD$FRic, fishFD$FEve, fishFD$FDiv,
                             fishFD$FDis, fishFD$RaoQ) #extract indices
  
colnames(important.indices) <- c("NumbSpecies", "FRic", "FDiv", "FEve", "FDis", "Rao")

FD_results <- important.indices %>% 
  as_tibble(rownames = "site_month") %>% 
  separate_wider_delim(site_month, delim = "_", names = c("site", "month"), cols_remove = FALSE) %>% 
  relocate(site_month) %>% 
  mutate(month = factor(month, labels = c("Apr", "May", "Jun", "Jul", "Aug", "Sept"))) %>% 
  mutate(site = factor(site, levels = c("FAM", "TUR", "COR", "MA", "WA", "HO", "SHR", "DOK", "LL", "TL", "PR", "EDG"))) %>% 
  replace(is.na(.), 0) #we only ever caught one penpoint gunnel at SHR in Apr (plus jellyfish)

FD_results.core <- FD_results %>% filter(!site %in% c("MA", "WA", "HO", "TL", "LL", "PR"))

#### plot indices for core sites ####

plot_index <- function (index){
  ggplot(data = FD_results.core, aes(x = site, y = .data[[index]], fill = site)) +
    geom_boxplot(outlier.shape = NA) + 
    theme_classic() +
    ylab(index) +
    theme(axis.title.x = element_blank())
}

index_plots <- lapply(names(FD_results[4:9]), plot_index)

index_plots[[1]] + index_plots[[2]] + index_plots[[3]] + index_plots[[4]] + index_plots[[5]] + index_plots[[6]] +
  plot_layout(ncol = 3, guides = 'collect')

ggsave("docs/figures/fish_FDcorepatch.png")

#### plot indices for all sites ####

plot_index_all <- function (index){
  FD_results %>% 
    filter(month == "Jun") %>% 
    ggplot(aes(x = site, y = .data[[index]], color = site)) +
      geom_point() + 
      geom_segment(aes(x = site, xend = site, y = 0, yend = .data[[index]])) +
      theme_classic() +
      ylab(index) +
      theme(axis.title.x = element_blank())
}

index_plots_all <- lapply(names(FD_results[4:9]), plot_index_all)

index_plots_all[[1]] + index_plots_all[[2]] + index_plots_all[[3]] + index_plots_all[[4]] + index_plots_all[[5]] + index_plots_all[[6]] +
  plot_layout(ncol = 3, guides = 'collect')

ggsave("docs/figures/fish_FDallpatch.png")

#### plot seasonality ####

plot_seasonal_index <- function (index){
ggplot(data = FD_results.core, aes(x = month, y = .data[[index]], group = site, color = site)) +
    geom_line() + 
    geom_point() +
    theme_classic() +
    ylab(index) +
    theme(axis.title.x = element_blank())
}

seasonal_index_plots <- lapply(names(FD_results[4:9]), plot_seasonal_index)

seasonal_index_plots[[1]] + seasonal_index_plots[[2]] + seasonal_index_plots[[3]] + seasonal_index_plots[[4]] + seasonal_index_plots[[5]] + seasonal_index_plots[[6]] +
  plot_layout(ncol = 2, guides = 'collect')

ggsave("docs/figures/fish_seasonalFDcorepatch.png")
            