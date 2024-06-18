library(here)
library(tidyverse)
library(mFD)
library(ggridges)
library(patchwork)

#load tidy fish data frame created in 01_tidy_data
load(here("data", "birds_tidy.Rdata")) 
load(here("data", "bird.list.Rdata")) #object created in 02_create_matrices

SOS_core_sites <- c("FAM", "TUR", "COR", "SHR", "DOK", "EDG")

#### bootstrap tidy dataframe ####

expand_species <- birds_tidy %>% 
  expand(nesting(year, month, site, ipa), spp_code) %>% 
  filter(!is.na(spp_code))

#create abundance matrix
mat_to_boot <- birds_tidy %>% 
  filter(!is.na(comm_name)) %>% 
  select(!c(day, comm_name)) %>% 
  full_join(expand_species) %>% 
  mutate(site = factor(site, levels = SOS_core_sites)) %>% 
  mutate(spp_sum = replace_na(spp_sum, 0)) %>% 
  arrange(spp_code)

#bootstrap
boot_list <- list()
nboot <- 999

for (i in seq_len(nboot)) {
  boot_list[[i]] <- mat_to_boot[sample(seq_len(nrow(mat_to_boot)), nrow(mat_to_boot), replace = TRUE),]
}

  
format_bird_L <- function(df) {
  
  bird_L <- df %>% 
    group_by(site, spp_code) %>% 
    summarize(spp_total_sum = sum(spp_sum)) %>% #sum across months/years
    ungroup() %>% 
    pivot_wider(names_from = spp_code, values_from = spp_total_sum, values_fill = 0) %>% 
    column_to_rownames(var = "site")
  
  return(bird_L)
}

boot_L <- lapply(boot_list, format_bird_L)

#### prep the functional trait space ####
traits.cat <- data.frame(trait_name = colnames(bird.list$trait),
                         trait_type = c("Q", "N", "N", "N", "N"))

#Species trait summary
traits_summary <- sp.tr.summary(tr_cat = traits.cat, 
                                sp_tr = bird.list$trait, 
                                stop_if_NA = T)

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

#extract the space 
trait_space <- space_quality$"details_fspaces"$"sp_pc_coord"

#### calculate the functional indices ####
alpha_indices <- list()
FD_results <- data.frame()

# with the mFD package
for (i in 1:length(boot_L)){
    
  alpha_indices <- alpha.fd.multidim(sp_faxes_coord = trait_space[ , c("PC1", "PC2", "PC3", "PC4")],
                                          asb_sp_w = data.matrix(boot_L[[i]]),
                                          ind_vect = c("fdis", "feve", "fric", "fdiv"),
                                          scaling = TRUE,
                                          check_input = TRUE,
                                          details_returned = TRUE)
  
  FD_values <- alpha_indices$"functional_diversity_indices"
  
  FD_results_df <- FD_values %>% 
    as_tibble(rownames = "site") 
  
  FD_results <- rbind(FD_results, FD_results_df)

}
#the mFD package uses ape::pcoa() which automatically removes negative eigenvalues rather than applying a correction

colnames(FD_results)[2:6] <- c("Species_Richness", "FDis", "FEve", "FRic", "FDiv")

FD_results <- FD_results %>% 
  mutate(region = ifelse(site %in% c("FAM", "TUR", "COR"), "North", "South"))

plot_site_index <- function (index){
  ggplot(data = FD_results, aes(x = .data[[index]], 
                                y = factor(site, levels = rev(SOS_core_sites)), 
                                fill = region, color = region)) +
    geom_density_ridges() + 
    theme_classic() +
    xlab(index) +
    theme(axis.title.y = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
}

index_plots <- lapply(names(FD_results[2:6]), plot_site_index)

index_plots[[1]] + index_plots[[2]] + index_plots[[3]] + index_plots[[4]] + index_plots[[5]] + guide_area() + 
  plot_layout(ncol = 3) + plot_layout(guides = "collect")

