library(tidyverse)
library(here)
library(ggrepel)
#library(FD)
library(mFD)
library(ggordiplots)
library(PNWColors)
library(patchwork)
library(vegan)

load(here("data", "fish.list.Rdata")) #object created in 02_tidy_data
here("analysis", "general_functions", "geb12299-sup-0002-si.r") %>% source()
here("analysis", "general_functions", "scree.r") %>% source()

traits.cat <- data.frame(trait_name = colnames(fish.list$trait),
                         trait_type = c("Q", "N", "N", "N", "N"))

#Species trait summary
traits_summary <- sp.tr.summary(tr_cat = traits.cat, sp_tr = fish.list$trait, stop_if_NA = T)

#create the trait space
dist_mat <- funct.dist(sp_tr = fish.list$trait, 
                       tr_cat = traits.cat,
                       metric = "gower",
                       scale_euclid = "scale_center",
                       ordinal_var = "classic",
                       weight_type = "equal",
                       stop_if_NA = TRUE)

#examine the quality of the potential functional spaces 
space_quality <- quality.fspaces(sp_dist = dist_mat,
                                 maxdim_pcoa = 10,
                                 deviation_weighting = "absolute",
                                 fdist_scaling = FALSE,
                                 fdendro = "ward.D2")

round(space_quality$"quality_fspaces",3) #lowest value is the best

#plot the trait space with the first 4 axes
plot_object <- space_quality$"details_fspaces"$"sp_pc_coord"

funct.space.plot(sp_faxes_coord = plot_object[ , c("PC1", "PC2", "PC3", "PC4")],
                 faxes = c("PC1", "PC2", "PC3", "PC4"))

alpha_indices <- alpha.fd.multidim(sp_faxes_coord = plot_object[ , c("PC1", "PC2", "PC3", "PC4")],
                                   asb_sp_w = data.matrix(fish.list$abund),
                                   ind_vect = c("fdis", "feve", "fric", "fdiv"),
                                   scaling = TRUE,
                                   check_input = TRUE,
                                   details_returned = TRUE)
#need to remove 2018 Jun Dok or decrease the number of axes
