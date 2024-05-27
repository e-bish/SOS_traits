library(tidyverse)
library(here)
library(ggrepel)
library(FD)
library(mFD)
library(ggordiplots)
library(PNWColors)
library(patchwork)
library(vegan)

load(here("data", "fish.list.Rdata")) #object created in 02_tidy_data

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
                       weight_type = "equal",
                       stop_if_NA = TRUE)

#examine the quality of the potential functional spaces 
space_quality <- quality.fspaces(sp_dist = dist_mat,
                                 maxdim_pcoa = 10,
                                 deviation_weighting = "absolute",
                                 fdist_scaling = FALSE,
                                 fdendro = "ward.D2")

round(space_quality$"quality_fspaces",3) #lowest value is the best (<.1 is good), meaning species pairs are accurately represented

#plot the trait space with the first 3 axes
plot_object <- space_quality$"details_fspaces"$"sp_pc_coord"

funct.space.plot(sp_faxes_coord = plot_object[ , c("PC1", "PC2", "PC3")],
                 faxes = c("PC1", "PC2", "PC3"))


#test for correlation between functional axes and traits
trait_axes <- traits.faxes.cor(
  sp_tr = fish.list$trait, 
  sp_faxes_coord = plot_object[ , c("PC1", "PC2", "PC3")], 
  plot = TRUE)

# Print traits with significant effect:
trait_axes$"tr_faxes_stat"[which(trait_axes$"tr_faxes_stat"$"p.value" < 0.05), ]

trait_axes$"tr_faxes_plot"

#plot the functional space
functional_space_plot <- mFD::funct.space.plot(
  sp_faxes_coord  = plot_object[ , c("PC1", "PC2", "PC3")],
  faxes           = c("PC1", "PC2", "PC3"),
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
               stand.x = TRUE, # we already log transformed the mean length variable so it doesn't need additional standardization
               corr = "none", #mFD package gives an explanation of why sqrt is misleading, and just removing the negative eigenvalues is preferred 
               m = 3,
               # calc.FGR = TRUE, 
               # clust.type = "ward.D2",
               calc.FDiv = TRUE, 
               print.pco = TRUE)

# with the mFD package
low_n_samples <- fish.list$abund %>% 
  decostand(method = "pa") %>% 
  mutate(n_spp = rowSums(.)) %>% 
  filter(n_spp < 4) %>% 
  rownames()

fish_L_filtered <- fish.list$abund %>% 
  filter(!rownames(.) %in% low_n_samples)

alpha_indices <- alpha.fd.multidim(sp_faxes_coord = plot_object[ , c("PC1", "PC2", "PC3")],
                                   asb_sp_w = data.matrix(fish_L_filtered),
                                   ind_vect = c("fdis", "feve", "fric", "fdiv"),
                                   scaling = TRUE,
                                   check_input = TRUE,
                                   details_returned = TRUE)
#the mFD package uses ape::pcoa() which automatically removes negative eigenvalues rather than applying a correction

FD_values <- alpha_indices$"functional_diversity_indices"

colnames(FD_values)[1:5] <- c("Species_Richness", "FDis", "FEve", "FRic", "FDiv")

FD_results <- FD_values %>% 
  as_tibble(rownames = "sample") %>% 
  separate_wider_delim(sample, delim = "_", names = c("year", "month", "site"), cols_remove = FALSE) %>% 
  relocate(sample) %>% 
  mutate(month = factor(month, labels = c("Apr", "May", "Jun", "Jul", "Aug", "Sept"))) %>% 
  mutate(site = factor(site, levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG"))) %>% 
  mutate(region = ifelse(site %in% c("FAM", "TUR", "COR"), "North", "South"), .after = site)

plot_index <- function (index, by){
    ggplot(data = FD_results, aes(x = .data[[by]], y = .data[[index]], fill = .data[[by]])) +
      geom_boxplot(outlier.shape = NA) + 
      theme_classic() +
      ylab(index) +
      theme(axis.title.x = element_blank(), 
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
      guides(fill="none")
}

index_plots <- lapply(names(FD_results[6:10]), plot_index, by = "site")

index_plots[[1]] + index_plots[[2]] + index_plots[[3]] + index_plots[[4]] + index_plots[[5]] +
  plot_layout(ncol = 3)

ggsave("docs/figures/fish_FDpatch.png")

index_plots <- lapply(names(FD_results[6:10]), plot_index, by = "region")

index_plots[[1]] + index_plots[[2]] + index_plots[[3]] + index_plots[[4]] + index_plots[[5]] +
  plot_layout(ncol = 3)

ggsave("docs/figures/fish_FDregionalpatch.png")


#### PERMANOVA ####
adonis2(FD_results[7:10] ~ FD_results$region, strata = FD_results$year, method = "euc")

CTRL <- permute::how(within = Within(type = "free"),
            blocks = year, #permutations are restricted within years
            plots = Plots(strata = FD_results$region, type = "free"),
            nperm = 9999,
            observed = TRUE)
#not working because there's an unbalanced design?

adonis2(FD_results[7:10] ~ year + FD_results$region,
        data = FD_results,
        method = "euclidean",
        permutations = CTRL)

kruskal.test(Species_Richness ~ site, data = FD_results)
kruskal.test(Species_Richness ~ region, data = FD_results)
#different

kruskal.test(FRic ~ site, data = FD_results)
kruskal.test(FRic ~ region, data = FD_results)
#different

kruskal.test(FEve ~ site, data = FD_results)
kruskal.test(FEve ~ region, data = FD_results)
#same

kruskal.test(FDiv ~ site, data = FD_results)
kruskal.test(FDiv ~ region, data = FD_results)
#same

kruskal.test(FDis ~ site, data = FD_results)
kruskal.test(FDis ~ region, data = FD_results)
#same

#### PCoA ####
North_indices <- which(FD_results$region == "North")
South_indices <- which(FD_results$site == "South")

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
  # geom_polygon(data = hulls, alpha = 0.3) +
  geom_point()
 