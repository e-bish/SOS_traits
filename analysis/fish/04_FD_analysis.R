library(tidyverse)
library(here)
library(ggrepel)
library(FD)
# library(mFD)
library(ggordiplots)
library(PNWColors)
library(patchwork)
library(vegan)

load(here("data", "fish.list.Rdata")) #object created in 03_create_matrices

#test the quality of the trait space
dist_mat <- gowdis(fish.list$trait.t, ord = "podani") #"classic" method matches mFD, which treats categorical variables as continuous
cor_dist_mat <- cailliez(dist_mat) #correction for negative eigenvalues

# #examine the quality of the potential functional spaces using the mFD package
# space_quality <- quality.fspaces(sp_dist = dist_mat,
#                                  maxdim_pcoa = 10,
#                                  deviation_weighting = c("absolute", "squared"), #setting this to squared and dist scaling to TRUE aligns with the original Maire et al. 2015 method
#                                  fdist_scaling = TRUE,
#                                  fdendro = "ward.D2")
# 
# round(space_quality$"quality_fspaces",3) #lowest value is the best (<.1 is good), meaning species pairs are accurately represented
# #aka the distances in euclidean space are accurately reflecting the gowers distances

n_axes_to_retain <- 4

#plot the trait space with the first 4 axes
# trait_space <- space_quality$"details_fspaces"$"sp_pc_coord"
# 
# funct.space.plot(sp_faxes_coord = trait_space[ , c("PC1", "PC2", "PC3", "PC4")],
#                   faxes = c("PC1", "PC2", "PC3", "PC4"))

# ggsave("docs/figures/fish_funct.space.plot.png")

# trait_space %>%
#   as_tibble(rownames = "species") %>%
#   ggplot(aes(x = PC1, y = PC2)) +
#   geom_text_repel(
#     label=rownames(trait_space),
#     max.time = 1,
#     max.overlaps = Inf) +
#   geom_point(color = "darkblue") +
#   theme_bw()

# ggsave("docs/figures/fish_pcoa.png")

#test for correlation between functional axes and traits
# trait_axes <- traits.faxes.cor(
#   sp_tr = fish.list$trait, 
#   sp_faxes_coord = trait_space[ , 1:n_axes_to_retain], 
#   plot = TRUE)
# 
# # Print traits with significant effect:
# trait_axes$"tr_faxes_stat"[which(trait_axes$"tr_faxes_stat"$"p.value" < 0.05), ]
# 
# trait_axes$"tr_faxes_plot"
# # ggsave("docs/figures/fish_trait_axes_plot.png")
# 
# #plot the functional space
# functional_space_plot <- mFD::funct.space.plot(
#   sp_faxes_coord  = trait_space[ , 1:4], #this function won't let you plot more than 4
#   faxes           = c("PC1", "PC2", "PC3", "PC4"),
#   name_file       = NULL,
#   faxes_nm        = NULL,
#   range_faxes     = c(NA, NA),
#   color_bg        = "grey95",
#   color_pool      = "darkgreen",
#   fill_pool       = "white",
#   shape_pool      = 21,
#   size_pool       = 1,
#   plot_ch         = TRUE,
#   color_ch        = "black",
#   fill_ch         = "white",
#   alpha_ch        = 0.5,
#   plot_vertices   = TRUE,
#   color_vert      = "blueviolet",
#   fill_vert       = "blueviolet",
#   shape_vert      = 23,
#   size_vert       = 1,
#   plot_sp_nm      = NULL,
#   nm_size         = 3,
#   nm_color        = "black",
#   nm_fontface     = "plain",
#   check_input     = TRUE)

#### calculate diversity indices ####

# with the FD package 
fishFD.ipa <- dbFD(x = fish.list$trait.t, #must be a df where character columns are factors or a distance matrix
               a = fish.list$abund.year,
               ord = "podani",
               corr = "cailliez", 
               m = n_axes_to_retain,
               calc.FDiv = TRUE, 
               print.pco = TRUE)

FD_values.ipa <- cbind(fishFD.ipa$nbsp, fishFD.ipa$FRic, fishFD.ipa$FEve, fishFD.ipa$FDiv, fishFD.ipa$FDis) #extract indices
colnames(FD_values.ipa) <- c("Species_Richness", "FRic", "FEve", "FDiv", "FDis")

FD_results.ipa <- FD_values.ipa %>% 
  as_tibble(rownames = "sample") %>% 
  separate_wider_delim(sample, delim = "_", names = c("site", "ipa"), cols_remove = FALSE) %>% 
  relocate(sample) %>% 
  mutate(site = factor(site, levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG"))) %>% 
  mutate(region = ifelse(site %in% c("FAM", "TUR", "COR"), "North", "South"), .after = site)

# save(FD_results, file = "data/FD_results.ipa.Rda")

plot_index <- function (index, by){
  ggplot(data = FD_results.ipa, aes(x = .data[[by]], 
                                y = .data[[index]], 
                                fill = .data[[by]])) +
    geom_boxplot() +
    geom_point(size = 2, color = "black") + 
    theme_classic() +
    ylab(index) +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
}

index_plots <- lapply(names(FD_results.ipa[5:9]), plot_index, by = "ipa")

index_plots[[1]] + index_plots[[2]] + index_plots[[3]] + index_plots[[4]] + index_plots[[5]] + guide_area() + 
  plot_layout(ncol = 3) + plot_layout(guides = "collect")

# ggsave("docs/figures/fish_FDipapatch.png")

#### test for differences in means ####

#test for multivariate normality
# library(mvn)
# mvn_test <- mvn(FD_results[5:9], mvnTest="mardia", multivariatePlot="qq")
# mvn_test$multivariateNormality #skewness does not agree; reject the null hypothesis
# #does not meet assumptions of multivariate normality

#permanova - nonparametric 
adonis2(FD_results.ipa[,c("Species_Richness","FDis", "FEve", "FRic", "FDiv")] ~ ipa, data = FD_results.ipa, method = "euc")
#cannot reject the null

#### recompute by site ####
fishFD <- dbFD(x = fish.list$trait.t, #must be a df where character columns are factors or a distance matrix
                   a = fish.list$abund,
                   ord = "podani",
                   corr = "cailliez", 
                   m = n_axes_to_retain,
                   calc.FDiv = TRUE, 
                   print.pco = TRUE)

FD_values <- cbind(fishFD$nbsp, fishFD$FRic, fishFD$FEve, fishFD$FDiv, fishFD$FDis) #extract indices
colnames(FD_values) <- c("Species_Richness", "FRic", "FEve", "FDiv", "FDis")

FD_results <- FD_values %>% 
  as_tibble(rownames = "sample") %>% 
  separate_wider_delim(sample, delim = "_", names = c("site", "year"), cols_remove = TRUE) %>% 
  mutate(site = factor(site, levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG")),
         year = factor(year, levels = c("2018", "2019", "2021", "2022"))) %>% 
  mutate(region = ifelse(site %in% c("FAM", "TUR", "COR"), "North", "South"), .after = site)

# save(FD_results, file = "data/FD_results.Rdata")  

FD_results %>% 
  pivot_longer(!c(site, region, year), names_to = "metric", values_to = "value") %>% 
  ggplot(aes(x = site, y = value, fill = region)) +
  geom_boxplot() +
  geom_point(show.legend = FALSE) +
  theme_classic() +
  facet_wrap(~metric, scales = "free_y")

#check to see if sites are different while controlling for year
adonis2(FD_results[,c("Species_Richness","FDis", "FEve", "FRic", "FDiv")] ~ site, 
        strata = FD_results$year, data = FD_results, method = "euc")
#they are different

#check to see if years are the same
adonis2(FD_results[,c("Species_Richness","FDis", "FEve", "FRic", "FDiv")] ~ year, 
        strata = FD_results$site, data = FD_results, method = "euc")
#they are

#check to see if there are differences by region
adonis2(FD_results[,c("Species_Richness","FDis", "FEve", "FRic", "FDiv")] ~ region, 
        strata = FD_results$year, data = FD_results, method = "euc")

adonis2(FD_results[,"FDis"] ~ site, strata = FD_results$year, data = FD_results, method = "euc")
adonis2(FD_results[,"FEve"] ~ site, strata = FD_results$year, data = FD_results, method = "euc")
adonis2(FD_results[,"FDiv"] ~ site, strata = FD_results$year, data = FD_results, method = "euc")
adonis2(FD_results[,"FRic"] ~ site, strata = FD_results$year, data = FD_results, method = "euc")

adonis2(FD_results[,"FEve"] ~ region, strata = FD_results$year, data = FD_results, method = "euc")


#### beta diversity ####
library(betapart)

abund_pa <- decostand(fish.list$abund, method = "pa")

# beta_obj <- functional.betapart.core(as.matrix(abund_pa), 
#                                      as.matrix(fishFD$x.axes[,1:4]), 
#                                      parallel = FALSE, #gives a weird error if you try to do this
#                                      warning.time = TRUE,
#                                      progress = TRUE)
# #this takes a reallllllly long time with 4 axes

beta_pairs <- functional.beta.pair(as.matrix(abund_pa), 
                                   as.matrix(fishFD$x.axes[,1:3]), #need 5 species for all samples to use 4 axes
                                   index.family = "jaccard")
