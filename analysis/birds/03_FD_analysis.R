library(tidyverse)
library(here)
library(ggrepel)
library(mFD)
library(FD)
library(PNWColors)
library(patchwork)
library(vegan)

load(here("data", "bird.list.Rdata")) #object created in 02_create_matrices

traits.cat <- data.frame(trait_name = colnames(bird.list$trait),
                         trait_type = c("Q", "N", "N", "N", "N"))

#Species trait summary
traits_summary <- sp.tr.summary(tr_cat = traits.cat, 
                                sp_tr = bird.list$trait, 
                                stop_if_NA = T)
traits_summary

#create the trait space
dist_mat <- funct.dist(sp_tr = bird.list$trait, 
                       tr_cat = traits.cat,
                       metric = "gower",
                       weight_type = "equal",
                       stop_if_NA = TRUE)

#examine the quality of the potential functional spaces 
space_quality <- quality.fspaces(sp_dist = dist_mat,
                                 maxdim_pcoa = 10,
                                 deviation_weighting = "absolute",
                                 fdist_scaling = FALSE,
                                 fdendro = "ward.D2")

round(space_quality$"quality_fspaces",3) #lowest value is the best (<.1 is good), meaning species pairs are accurately represented
#aka the distances in euclidean space are accurately reflecting the gowers distances

#plot the trait space with the first 4 axes
plot_object <- space_quality$"details_fspaces"$"sp_pc_coord"

funct.space.plot(sp_faxes_coord = plot_object[ , c("PC1", "PC2", "PC3", "PC4")],
                 faxes = c("PC1", "PC2", "PC3", "PC4"))

#test for correlation between functional axes and traits
trait_axes <- traits.faxes.cor(
  sp_tr = bird.list$trait, 
  sp_faxes_coord = plot_object[ , c("PC1", "PC2", "PC3", "PC4")], 
  plot = TRUE)

trait_axes$"tr_faxes_plot"

# with the FD package 
birdFD <- dbFD(x = bird.list$trait, #must be a df where character columns are factors
               a = bird.list$abund,
               stand.x = TRUE, # we already log transformed the mean length variable so it doesn't need additional standardization
               corr = "none", #mFD package gives an explanation of why sqrt is misleading, and just removing the negative eigenvalues is preferred 
               m = 4,
               # calc.FGR = TRUE, 
               # clust.type = "ward.D2",
               calc.FDiv = TRUE, 
               print.pco = TRUE)

# with the mFD package
alpha_indices <- alpha.fd.multidim(sp_faxes_coord = plot_object[ , c("PC1", "PC2", "PC3", "PC4")],
                                   asb_sp_w = data.matrix(bird.list$abund),
                                   ind_vect = c("fdis", "feve", "fric", "fdiv"),
                                   scaling = TRUE,
                                   check_input = TRUE,
                                   details_returned = TRUE)
#the mFD package uses ape::pcoa() which automatically removes negative eigenvalues rather than applying a correction

FD_values <- alpha_indices$"functional_diversity_indices"

colnames(FD_values)[1:5] <- c("Species_Richness", "FDis", "FEve", "FRic", "FDiv")

FD_results <- FD_values %>% 
  as_tibble(rownames = "sample") %>% 
  separate_wider_delim(sample, delim = "_", names = c("site", "ipa"), cols_remove = FALSE) %>% 
  relocate(sample) %>% 
  mutate(site = factor(site, levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG"))) %>% 
  mutate(region = ifelse(site %in% c("FAM", "TUR", "COR"), "North", "South"), .after = site)

plot_site_index <- function (index){
  ggplot(data = FD_results, aes(x = site, 
                                y = .data[[index]], 
                                color = site,
                                shape = ipa)) +
    geom_point(size = 3) + 
    theme_classic() +
    ylab(index) +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
}

index_plots <- lapply(names(FD_results[5:9]), plot_site_index)

index_plots[[1]] + index_plots[[2]] + index_plots[[3]] + index_plots[[4]] + index_plots[[5]] + guide_area() + 
  plot_layout(ncol = 3) + plot_layout(guides = "collect")

plot_index <- function (index, by){
  ggplot(data = FD_results, aes(x = .data[[by]], 
                                y = .data[[index]], 
                                fill = .data[[by]])) +
    geom_boxplot() + 
    theme_classic() +
    ylab(index) +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
}

index_plots <- lapply(names(FD_results[5:9]), plot_index, by = "region")

index_plots[[1]] + index_plots[[2]] + index_plots[[3]] + index_plots[[4]] + index_plots[[5]] + guide_area() + 
  plot_layout(ncol = 3) + plot_layout(guides = "collect")

index_plots <- lapply(names(FD_results[5:9]), plot_index, by = "ipa")

index_plots[[1]] + index_plots[[2]] + index_plots[[3]] + index_plots[[4]] + index_plots[[5]] + guide_area() + 
  plot_layout(ncol = 3) + plot_layout(guides = "collect")

#### rank sum test ####

kruskal.test(Species_Richness ~ site, data = FD_results) #same
kruskal.test(Species_Richness ~ region, data = FD_results) #different
kruskal.test(Species_Richness ~ ipa, data = FD_results) #same


kruskal.test(FRic ~ site, data = FD_results) #same 
kruskal.test(FRic ~ region, data = FD_results) #same
kruskal.test(FRic ~ ipa, data = FD_results) #same


kruskal.test(FEve ~ site, data = FD_results) #different
kruskal.test(FEve ~ region, data = FD_results) #different
kruskal.test(FEve ~ ipa, data = FD_results) #same


kruskal.test(FDiv ~ site, data = FD_results) #same
kruskal.test(FDiv ~ region, data = FD_results) #different
kruskal.test(FDiv ~ ipa, data = FD_results) #same


#view pcoa
pcoa <- cmdscale(bird.gowdist)
colnames(pcoa) <- c("pcoa1", "pcoa2")

bird_pcoa <- pcoa %>%
  as_tibble(rownames = "species") %>%
  ggplot(aes(x = pcoa1, y = pcoa2)) +
  geom_text_repel(
    label=rownames(pcoa),
    max.time = 1,
    max.overlaps = Inf) +
  geom_point(color = "darkblue") +
  theme_bw()

ggsave(plot = bird_pcoa, "docs/figures/bird_pcoa.png")

#### test the quality of the functional space ####
space_qual <- qual_funct_space(bird.traits, nbdim = 6, metric = "Gower")
#needs to be df with no NAs
space_qual$meanSD #  (we're shooting for 0 here)
#"The mSD is 0 when the functional space perfectly represents the initial distance and increases as pairs of species become less represented in the functional space (Maire et al., 2015)."

#### test each of these by plotting to see how many groups to use ####
test1 <- hclust(bird.gowdist, method = "average") #meh
test2 <- hclust(bird.gowdist, method = "complete") #pretty good
test3 <- hclust(bird.gowdist, method = "median") #terrible
bird.hclust <- hclust(bird.gowdist, method = "ward.D2") #better!
plot(bird.hclust, hang = -1,
     main = "Method = Ward’s minimum variance",
     xlab = "Species",
     sub = "BC dissimilarity; ward.D2 linkage")
scree(hclust.obj = bird.hclust) # not a solid answer here but 6-8 is generally at the bottom of the elbow, so we'll try 6 groups

#### FD indices calculation ####################################################
birdFD <- dbFD(bird.traits, #must be a df where character columns are factors
               bird.abund, 
                 m = 3, #seems to keep 4 axis no matter how many you specify?
                 calc.FGR = TRUE, 
                 clust.type = "ward.D2",
                 calc.FRic = FALSE, #turn off if you have less species than traits
                 calc.FDiv = TRUE, #won't return a value if there are only categorical traits
                 print.pco = TRUE)
#it will prompt to see if you want to separate by groups (g) and ask how many (6)

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
            