library(tidyverse)
library(here)
library(ggrepel)
library(FD)
library(mFD)
library(fundiversity)
library(ggordiplots)
library(PNWColors)
library(patchwork)
library(vegan)

load(here("data", "fish.list.Rdata")) #object created in 03_create_matrices

# traits.cat <- data.frame(trait_name = colnames(fish.list$trait),
#                          trait_type = c("Q", "N", "N", "N", "N"))

#Species trait summary
# traits_summary <- sp.tr.summary(tr_cat = traits.cat, 
#                                 sp_tr = fish.list$trait, 
#                                 stop_if_NA = T)
# traits_summary

#create the trait space
#don't need to specify how to deal with ordinal variables because ours are factors
# dist_mat <- funct.dist(sp_tr = fish.list$trait,
#                        tr_cat = traits.cat,
#                        metric = "gower",
#                        weight_type = "equal",
#                        stop_if_NA = TRUE)

#using the FD package
# dist_mat <- gowdis(fish.list$trait, ord = "classic")
# funct_space <- ape::pcoa(dist_mat)$vectors
# 
# euc_distances1 <- vegdist(funct_space[,1], method = "euc")
# euc_distances2 <- vegdist(funct_space[,1:2], method = "euc")
# euc_distances3 <- vegdist(funct_space[,1:3], method = "euc")
# euc_distances4 <- vegdist(funct_space[,1:4], method = "euc")
# euc_distances5 <- vegdist(funct_space[,1:5], method = "euc")
# euc_distances6 <- vegdist(funct_space[,1:6], method = "euc")
# 
# euc_list <- list(euc_distances1, euc_distances2, euc_distances3, euc_distances4, euc_distances5, euc_distances6)
# 
# calculate_space_qual <- function(euc_distances) {
#   std_funct_space <- euc_distances / max(euc_distances2) * max(dist_mat)
#   
#   S <- 42
#   
#   # sum(abs(dist_mat - std_funct_space)) / ((S * S-1)/2) #maire 2015 used squared deviation but mFD uses absolute differences
#   mean(abs(dist_mat - std_funct_space)) #maire 2015 used squared deviation but mFD uses absolute differences
#   
# }
# 
# lapply(euc_list, calculate_space_qual)

#examine the quality of the potential functional spaces using the mFD package
space_quality <- quality.fspaces(sp_dist = dist_mat,
                                 maxdim_pcoa = 10,
                                 deviation_weighting = "absolute",
                                 fdist_scaling = FALSE,
                                 fdendro = "ward.D2")

round(space_quality$"quality_fspaces",3) #lowest value is the best (<.1 is good), meaning species pairs are accurately represented
#aka the distances in euclidean space are accurately reflecting the gowers distances

n_axes_to_retain <- 5

#plot the trait space with the first 4 axes
trait_space <- space_quality$"details_fspaces"$"sp_pc_coord"

# funct.space.plot(sp_faxes_coord = trait_space[ , c("PC1", "PC2", "PC3", "PC4")],
#                  faxes = c("PC1", "PC2", "PC3", "PC4"))

# ggsave("docs/figures/fish_funct.space.plot.png")

trait_space %>%
  as_tibble(rownames = "species") %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_text_repel(  
    label=rownames(trait_space),
    max.time = 1,
    max.overlaps = Inf) +
  geom_point(color = "darkblue") +
  theme_bw()

# ggsave("docs/figures/fish_pcoa.png")

#test for correlation between functional axes and traits
trait_axes <- traits.faxes.cor(
  sp_tr = fish.list$trait, 
  sp_faxes_coord = plot_object[ , 1:n_axes_to_retain], 
  plot = TRUE)

# Print traits with significant effect:
trait_axes$"tr_faxes_stat"[which(trait_axes$"tr_faxes_stat"$"p.value" < 0.05), ]

trait_axes$"tr_faxes_plot"
# ggsave("docs/figures/fish_trait_axes_plot.png")

#plot the functional space
functional_space_plot <- mFD::funct.space.plot(
  sp_faxes_coord  = plot_object[ , 1:4], #this function won't let you plot more than 4
  faxes           = c("PC1", "PC2", "PC3", "PC4"),
  name_file       = NULL,
  faxes_nm        = NULL,
  range_faxes     = c(NA, NA),
  color_bg        = "grey95",
  color_pool      = "darkgreen",
  fill_pool       = "white",
  shape_pool      = 21,
  size_pool       = 1,
  plot_ch         = TRUE,
  color_ch        = "black",
  fill_ch         = "white",
  alpha_ch        = 0.5,
  plot_vertices   = TRUE,
  color_vert      = "blueviolet",
  fill_vert       = "blueviolet",
  shape_vert      = 23,
  size_vert       = 1,
  plot_sp_nm      = NULL,
  nm_size         = 3,
  nm_color        = "black",
  nm_fontface     = "plain",
  check_input     = TRUE)

#### calculate diversity indices ####
#using fundiversity
# FRic <- fd_fric(trait_space[,1:n_axes_to_retain], fish.list$abund, stand = TRUE)
# FEve <- fd_feve(trait_space[,1:n_axes_to_retain], fish.list$abund)
# FDiv <- fd_fdiv(trait_space[,1:n_axes_to_retain], fish.list$abund)

#fdis doesn't work for some reason
# FDis <- fd_fdis(trait_space[,1:n_axes_to_retain], fish.list$abund)
# 
# FD_values <- full_join(FRic, FEve) %>% 
#   left_join(FDiv) %>% 
#   full_join(FDis)

alpha_indices <- alpha.fd.multidim(sp_faxes_coord = plot_object[ , 1:n_axes_to_retain],
                                   asb_sp_w = data.matrix(fish.list$abund),
                                   ind_vect = c("fdis", "feve", "fric", "fdiv"),
                                   scaling = TRUE,
                                   check_input = TRUE,
                                   details_returned = TRUE)
# the mFD package uses ape::pcoa() which automatically removes negative eigenvalues rather than applying a correction

FD_values <- alpha_indices$"functional_diversity_indices"

colnames(FD_values)[1:5] <- c("Species_Richness", "FDis", "FEve", "FRic", "FDiv")

FD_results <- FD_values %>% 
  as_tibble(rownames = "sample") %>% 
  separate_wider_delim(sample, delim = "_", names = c("site", "ipa"), cols_remove = FALSE) %>% 
  relocate(sample) %>% 
  mutate(site = factor(site, levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG"))) %>% 
  mutate(region = ifelse(site %in% c("FAM", "TUR", "COR"), "North", "South"), .after = site)

# ggsave("docs/figures/fish_FDpatch.png")

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

index_plots <- lapply(names(FD_results[5:9]), plot_index, by = "ipa")

index_plots[[1]] + index_plots[[2]] + index_plots[[3]] + index_plots[[4]] + index_plots[[5]] + guide_area() + 
  plot_layout(ncol = 3) + plot_layout(guides = "collect")

ggsave("docs/figures/fish_FDipapatch.png")

#### test for equal variances ####
ipa_test <- manova(cbind(Species_Richness, FDis, FEve, FRic, FDiv) ~ ipa, data = FD_results)
summary(ipa_test) #cannot reject the null
summary.aov(ipa_test)

site_test <- manova(cbind(Species_Richness, FDis, FEve, FRic, FDiv) ~ site, data = FD_results)
summary(site_test) #reject the null
summary.aov(site_test)

region_test <- manova(cbind(Species_Richness, FDis, FEve, FRic, FDiv) ~ region, data = FD_results)
summary(region_test) #cannot reject the null
summary.aov(region_test)


