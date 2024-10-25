library(tidyverse)
library(here)
library(ggrepel)
library(FD)
library(mFD)
library(ggordiplots)
library(PNWColors)
library(patchwork)
library(vegan)

load(here("data", "fish.list.balanced.Rdata")) #object created in 02_tidy_data
#current version is setup by site_ipa

fish.list <- fish.list.balanced

traits.cat <- data.frame(trait_name = colnames(fish.list$trait),
                         trait_type = c("Q", "N", "N", "N", "N"))

#Species trait summary
traits_summary <- sp.tr.summary(tr_cat = traits.cat, 
                                sp_tr = fish.list$trait, 
                                stop_if_NA = T)
traits_summary

#create the trait space
dist_mat <- funct.dist(sp_tr = fish.list$trait, 
                       tr_cat = traits.cat,
                       metric = "gower",
                       ordinal_var   = "podani",
                       weight_type = "equal",
                       stop_if_NA = TRUE)

#examine the quality of the potential functional spaces 
space_quality <- quality.fspaces(sp_dist = dist_mat,
                                 maxdim_pcoa = 10,
                                 deviation_weighting = c("absolute","squared"),
                                 fdist_scaling = FALSE,
                                 fdendro = "ward.D2")

round(space_quality$"quality_fspaces",3) #lowest value is the best (<.1 is good), meaning species pairs are accurately represented
#aka the distances in euclidean space are accurately reflecting the gowers distances

#plot the trait space with the first 4 axes
plot_object <- space_quality$"details_fspaces"$"sp_pc_coord"

funct.space.plot(sp_faxes_coord = plot_object[ , c("PC1", "PC2", "PC3", "PC4")],
                 faxes = c("PC1", "PC2", "PC3", "PC4"))

# ggsave("docs/figures/fish_funct.space.plot.png")

plot_object %>%
  as_tibble(rownames = "species") %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_text_repel(  
    label=rownames(plot_object),
    max.time = 1,
    max.overlaps = Inf) +
  geom_point(color = "darkblue") +
  theme_bw()

# ggsave("docs/figures/fish_pcoa.png")

# recreate the PCoA manually
fish.gowdist <- gowdis(fish.list$trait)
pcoa <- cmdscale(fish.gowdist, add = TRUE) #cailliez correction to recreate their plot, though the mFD package drops negative eigenvalues before calculating FD indices
colnames(pcoa$points) <- c("pcoa1", "pcoa2")

pcoa$points %>%
  as_tibble(rownames = "species") %>%
  ggplot(aes(x = pcoa1, y = pcoa2)) +
  geom_text_repel(
    label=rownames(pcoa$points),
    max.time = 1,
    max.overlaps = Inf) +
  geom_point(color = "darkblue") +
  theme_bw()

#test for correlation between functional axes and traits
trait_axes <- traits.faxes.cor(
  sp_tr = fish.list$trait, 
  sp_faxes_coord = plot_object[ , c("PC1", "PC2", "PC3", "PC4")], 
  plot = TRUE)

# Print traits with significant effect:
trait_axes$"tr_faxes_stat"[which(trait_axes$"tr_faxes_stat"$"p.value" < 0.05), ]

trait_axes$"tr_faxes_plot"
# ggsave("docs/figures/fish_trait_axes_plot.png")

#plot the functional space
functional_space_plot <- mFD::funct.space.plot(
  sp_faxes_coord  = plot_object[ , c("PC1", "PC2", "PC3", "PC4")],
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

# here("analysis", "general_functions", "geb12299-sup-0002-si.r") %>% source()
# space_qual <- qual_funct_space(fish_L_filtered, nbdim = 3, metric = "Gower")

# with the FD package 
fishFD <- dbFD(x = fish.list$trait, #must be a df where character columns are factors
               a = fish.list$abund,
               stand.x = TRUE, # revisit this
               ord = "podani",
               corr = "none", #mFD package gives an explanation of why sqrt is misleading, and just removing the negative eigenvalues is preferred 
               m = 4,
               # calc.FGR = TRUE, 
               # clust.type = "ward.D2",
               calc.FDiv = TRUE, 
               print.pco = TRUE)

# with the mFD package
alpha_indices <- alpha.fd.multidim(sp_faxes_coord = plot_object[ , c("PC1", "PC2", "PC3", "PC4")], #even though 4d was higher quality, 3d allows us to retain more samples and is still of high enough quality
                                   asb_sp_w = data.matrix(fish_L),
                                   ind_vect = c("fdis", "feve", "fric", "fdiv"),
                                   scaling = TRUE,
                                   check_input = TRUE,
                                   details_returned = TRUE)
#the mFD package uses ape::pcoa() which automatically removes negative eigenvalues rather than applying a correction

FD_values <- alpha_indices$"functional_diversity_indices"

colnames(FD_values)[1:5] <- c("Species_Richness", "FDis", "FEve", "FRic", "FDiv")

FD_results <- FD_values %>% 
  as_tibble(rownames = "site_ipa") %>% 
  separate_wider_delim(site_ipa, delim = "_", names = c("site", "ipa"), cols_remove = FALSE) %>% 
  relocate(site_ipa) %>% 
  # mutate(month = factor(month, labels = c("Apr", "May", "Jun", "Jul", "Aug", "Sept"))) %>% 
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

# ggsave("docs/figures/fish_FDpatch18.png")

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

# ggsave("docs/figures/fish_FDregionalpatch18.png")

index_plots <- lapply(names(FD_results[5:9]), plot_index, by = "ipa")

index_plots[[1]] + index_plots[[2]] + index_plots[[3]] + index_plots[[4]] + index_plots[[5]] + guide_area() + 
  plot_layout(ncol = 3) + plot_layout(guides = "collect")

# ggsave("docs/figures/fish_FDipapatch18.png")

#### PERMANOVA ####
set.seed(1993)
#block: region, the highest level of permutation. samples are never swapped between blocks
#plot: site, level below blocks. plots can be permuted while blocks cannot
#strata: ipa, describes the grouping of variables at the plot level

#are sites different within regions?
#because TUR has an uneven number of IPAs we need to drop that site
FD_results2 <- FD_results %>% filter(!site == "TUR")

CTRL <- permute::how(within = Within(type = "none"), #this type allows us to permute plots (sites) without permuting their ipas
                     plots = Plots(strata = FD_results2$site, type = "free"), #this is the level we're permuting
                     nperm = 9999,
                     observed = TRUE)

adonis2(FD_results2[5:9] ~ FD_results2$region, FD_results2$site, #the blocking factor is typically listed first in type 1 SS to account for the most variation
        data = FD_results2,
        method = "euclidean",
        permutations = CTRL)


# # are regions different (tested by permuting sites while controlling for ipas)?
# CTRL <- permute::how(within = Within(type = "none"), #this permutes the ordering within plots
#                      blocks = FD_results$ipa, #permutations are restricted between ipas
#                      nperm = 9999,
#                      observed = TRUE)

adonis2(FD_results[5:9] ~ FD_results$region,
        data = FD_results,
        method = "euclidean",
        permutations = CTRL)
#region explains 36.2% of the variation between sites with a significance level of p = 0.0076





kruskal.test(Species_Richness ~ site, data = FD_results) #different
kruskal.test(Species_Richness ~ region, data = FD_results) #different
kruskal.test(Species_Richness ~ ipa, data = FD_results) #same

kruskal.test(FRic ~ site, data = FD_results) #same
kruskal.test(FRic ~ region, data = FD_results) #different
kruskal.test(FRic ~ ipa, data = FD_results) #same

kruskal.test(FEve ~ site, data = FD_results) #same 
kruskal.test(FEve ~ region, data = FD_results) #same
kruskal.test(FEve ~ ipa, data = FD_results) #same

kruskal.test(FDiv ~ site, data = FD_results) #same 
kruskal.test(FDiv ~ region, data = FD_results) #different
kruskal.test(FDiv ~ ipa, data = FD_results) #same

kruskal.test(FDis ~ site, data = FD_results) #different with 4d
kruskal.test(FDis ~ region, data = FD_results) #same, but borderline
kruskal.test(FDis ~ ipa, data = FD_results) #same

#### PCoA by region ####
North_indices <- which(FD_results$region == "North")
South_indices <- which(FD_results$region == "South")

North_df <- fish_L_filtered %>% 
  filter(row_number() %in% North_indices) %>% 
  select_if(colSums(.) != 0)

North_spp <- colnames(North_df)

North_coords <- plot_object %>% 
  as_tibble(rownames = "species") %>% 
  filter(species %in% North_spp) 

South_df <- fish_L_filtered %>% 
  filter(row_number() %in% South_indices) %>% 
  select_if(colSums(.) != 0)

South_spp <- colnames(South_df)

South_coords <- plot_object %>% 
  as_tibble(rownames = "species") %>% 
  filter(species %in% South_spp) 

all_coords <- bind_rows(North_coords, South_coords, .id = "Region") %>% 
  mutate(Region = ifelse(Region == 1, "North", "South")) 

hulls <- all_coords %>% 
  group_by(Region) %>% 
  do(.[chull(.[3:4]), ])

ggplot(all_coords, aes(x = PC1, y = PC2, col = factor(Region), fill = factor(Region)))+
  geom_polygon(data = hulls, alpha = 0.3) +
  geom_point()

#### test for differences in beta diversity ####
beta_fd_indices <- mFD::beta.fd.multidim(
  sp_faxes_coord   = plot_object[ , c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_occ       = decostand(fish_L, method = "pa"),
  check_input      = TRUE,
  beta_family      = c("Jaccard"),
  details_returned = TRUE)

jac_diss <- beta_fd_indices[[1]][[1]]
jac_turn <- beta_fd_indices[[1]][[2]]
jac_nest <- beta_fd_indices[[1]][[3]]

plot(jac_diss)

beta_plot_fruits <- mFD::beta.multidim.plot(
  output_beta_fd_multidim = beta_fd_indices,
  plot_asb_nm             = c("EDG_Armored", "COR_Natural"),
  beta_family             = c("Jaccard"),
  plot_sp_nm              = NULL,
  faxes                   = paste0("PC", 1:4),
  # name_file               = NULL,
  # faxes_nm                = NULL,
  # range_faxes             = c(NA, NA),
  # color_bg                = "grey95",
  # shape_sp                = c("pool" = 3.0, asb1 = 22, asb2 = 21),
  # size_sp                 = c("pool" = 0.8, asb1 =  1, asb2 =  1),
  # color_sp                = c("pool" = "grey50", asb1 = "blue", asb2 = "red"),
  # fill_sp                 = c("pool" = NA, asb1 = "white", asb2 = "white"),
  # fill_vert               = c("pool" = NA, asb1 = "blue", asb2 = "red"),
  # color_ch                = c("pool" = NA, asb1 = "blue", asb2 = "red"),
  # fill_ch                 = c("pool" = "white", asb1 = "blue", asb2 = "red"),
  # alpha_ch                = c("pool" = 1, asb1 = 0.3, asb2 = 0.3),
  # nm_size                 = 3,
  # nm_color                = "black",
  # nm_fontface             = "plain",
  check_input             = TRUE)
