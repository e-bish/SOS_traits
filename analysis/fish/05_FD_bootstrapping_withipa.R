library(here)
library(tidyverse)
library(janitor)
library(rsample)
library(mFD)
library(FD)
library(ggridges)
library(patchwork)

set.seed(1993)

#load tidy fish data frame created in 02_tidy_data
load(here("data", "net_tidy.Rdata")) 
load(here("data", "fish.list.Rdata")) #object created in 03_create_matrices

SOS_core_sites <- c("FAM", "TUR", "COR", "SHR", "DOK", "EDG")

#### bootstrap tidy dataframe ####

expand_species <- net_tidy %>% 
  expand(nesting(year, month, site, ipa), ComName) %>% 
  filter(!is.na(ComName))

#create abundance matrix
mat_to_boot <- net_tidy %>% 
  filter(!is.na(ComName)) %>% 
  group_by(year, month, day, site, ipa, ComName) %>% 
  summarize(spp_sum = sum(species_count)) %>% #sum across depth stations on each day
  ungroup() %>% 
  select(!day) %>% #don't need this anymore
  full_join(expand_species) %>% 
  mutate(site = factor(site, levels = SOS_core_sites)) %>% 
  mutate(spp_sum = replace_na(spp_sum, 0)) %>% 
  mutate(ComName = replace(ComName, ComName == "Pacific Sandfish", "Pacific sandfish")) %>% #if it's capitolized it gets confused about the proper order
  arrange(ComName) 
#bootstrap & keep equal proportions for each month to maintain seasonal variability
#using the strata argument here instead of the group option, because the group option results in 
#some assessment datasets with zero rows (meaning that the analysis dataset is the exact same as the raw data)
#the strata argument works well when categorical variables are unbalanced but you want similar proportions
boot_obj <- bootstraps(mat_to_boot, strata = month, times = 999)

# for (i in seq_along(boot_obj$splits)) {
#   print(boot_obj$splits[[i]])
# }

#store the bootstrap dataframes in a list
boot_list <- list()
for (i in seq_along(boot_obj$splits)) {
  resample <- boot_obj$splits[[i]]
  boot_list[[i]] <- analysis(resample)
  
}

format_fish_L <- function(df) {
  
  fish_L <- df %>% 
    group_by(site, ipa, ComName) %>% 
    summarize(spp_total_avg = mean(spp_sum)) %>% #average across months/years
    ungroup() %>% 
    pivot_wider(names_from = ComName, values_from = spp_total_avg, values_fill = 0) %>% 
    mutate(site_ipa = paste(site, ipa, sep = "_"), .after = ipa, .keep = "unused") %>% 
    column_to_rownames(var = "site_ipa") %>% 
    clean_names() 
  
  return(fish_L)
}

#format the bootstrapped data frame into properly formatted matrices
boot_L <- lapply(boot_list, format_fish_L)

#remove species that aren't represented in assemblages 
remove_missing_spp <- function(df) {
  filtered_df <- df %>% 
    select_if(colSums(.) != 0)
  
  return(filtered_df)
}

boot_L_filtered <- lapply(boot_L, remove_missing_spp)

#### prep the functional trait space ####
traits.cat <- data.frame(trait_name = colnames(fish.list$trait),
                         trait_type = c("Q", "N", "N", "N", "N"))

#Species trait summary
traits_summary <- sp.tr.summary(tr_cat = traits.cat, 
                                sp_tr = fish.list$trait, 
                                stop_if_NA = T)

#create a list of trait matrices based on species represented in the bootstrapped L matrices
fish_Q_list <- list() 
dist_mat <- list()
space_quality <- list()
trait_space <- list()

for (i in 1:length(boot_L)) {
  
  fish_Q_list[[i]] <- fish.list$trait.t %>%
    rownames_to_column(var = "species") %>%
    filter(species %in% colnames(boot_L_filtered[[i]])) %>%
    column_to_rownames(var = "species")

  #create the trait space
  dist_mat[[i]] <- gowdis(fish_Q_list[[i]], ord = "podani")
  dist_mat[[i]] <- cailliez(dist_mat[[i]])

  #examine the quality of the potential functional spaces 
  space_quality[[i]] <- quality.fspaces(sp_dist = dist_mat[[i]],
                                        maxdim_pcoa = 10,
                                        deviation_weighting = "absolute",
                                        fdist_scaling = FALSE,
                                        fdendro = "ward.D2")
  
  #extract the space 
  trait_space[[i]] <- space_quality[[i]]$"details_fspaces"$"sp_pc_coord"
  
}

#look at the distribution of best space qualities
space_qual_df <- data.frame(matrix(ncol = 2))
names(space_qual_df) <- c("n_axes", "mad")

for (i in 1:length(space_quality)) {
  space_qual <- space_quality[[i]][1] %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "n_axes") 
    # slice(which.min(mad))
  
  space_qual_df <- rbind(space_qual_df, space_qual)
}

space_qual_df %>% 
  filter(!is.na(n_axes)) %>% 
  group_by(n_axes) %>% 
  summarize(n = n(), min = min(mad), max = max(mad), avg = mean(mad))
#not sure why these are so terrible now, maybe it's the correction?
#if not separating by ipa, 5 axes is more often the best trait space but 4 is close behind

space_qual_df %>% 
  filter(!is.na(n_axes)) %>% 
  ggplot() + 
  geom_histogram(aes(x = mad), fill = "lightblue", color = "skyblue") +
  xlab("Quality of the Trait Space") +
  theme_classic()

#### calculate the functional indices ####
# with the FD package
gowdist.list <- lapply(fish_Q_list, gowdis, ord = "podani")
FD_results <- list()
# CWM_results <- list() #gives an error for one of the bootstrapped matrices when you try to use this for getting the CWM

for (i in 1:length(boot_L)){

  fishFD <- dbFD(x = gowdist.list[[i]], #must be a distance object or df where character columns are factors
                 a = data.matrix(boot_L_filtered[[i]]),
                 corr = "cailliez", 
                 m = 5,
                 calc.FDiv = TRUE,
                 print.pco = FALSE)

  FD_values <- cbind(fishFD$nbsp, fishFD$FRic, fishFD$FEve, fishFD$FDiv, fishFD$FDis) #extract indices

  FD_results_df <- FD_values %>%
    as_tibble(rownames = "site_ipa")

  FD_results <- rbind(FD_results, FD_results_df)
  # CWM_results[[i]] <- fishFD$CWM
}

colnames(FD_results)[2:7] <- c("Species_Richness", "FRic", "FEve", "FDiv", "FDis")

FD_results <- FD_results %>%
  separate_wider_delim(site_ipa, delim = "_", names = c("site", "ipa"), cols_remove = FALSE) %>% 
  mutate(site = factor(site, levels = SOS_core_sites)) %>%
  mutate(region = ifelse(site %in% c("FAM", "TUR", "COR"), "North", "South"))

FD_boot_results <- FD_results
# save(FD_boot_results, file = "data/FD_boot_results.Rda")

#test for differences in calculations between packages
# all.equal(FD_results, FD_results_v2)

index_names <- c("Species Richness","F. Richness","F. Eveness","F. Divergence","F. Dispersion")

plot_site_index <- function (index){
  ggplot(data = FD_results, aes(x = .data[[index]], 
                                y = factor(site, levels = rev(SOS_core_sites)), 
                                fill = region, color = region)) +
    geom_density_ridges(alpha = 0.9) + 
    theme_classic() +
    theme(axis.title.y = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
}

index_plots <- lapply(names(FD_results[5:8]), plot_site_index)

for (i in 1:length(index_plots)) {
  index_plots[[i]] <- index_plots[[i]] + xlab(index_names[i+1])
}

#FEve and FDiv are scaled between 0,1
index_plots$"F. Evenness" <- index_plots$"F. Evenness" + xlim(0,1)
index_plots$"F. Divergence" <- index_plots$"F. Divergence" + xlim(0,1)

index_plots[[5]] <- ggplot(data = FD_results, 
                           aes(x = Species_Richness, 
                               y = factor(site, levels = rev(SOS_core_sites)), 
                               color = region)) +
  geom_density_ridges2(aes(fill = region), 
                       stat = "binline", 
                       binwidth = 1, 
                       scale = 0.95, 
                       alpha = 0.9) +
  theme_classic() +
  xlab("Species Richness") +
  theme(axis.title.y = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

index_plots[[5]]  + index_plots[[1]]  + guide_area() + index_plots[[2]] + index_plots[[3]] + index_plots[[4]] + 
  plot_layout(ncol = 3, guides = "collect")

# ggsave("docs/figures/fish_FDbootpatch.png")

#### test for differences ####

FD_boot_means <- FD_boot_results %>% 
  # group_by(site) %>% 
  summarize(across(where(is.numeric), mean))

#### taxonomic diversity ####

pivot_boot_L <- function(fish_L) {
  fish_L_long <- fish_L %>% 
    rownames_to_column("sample") %>% 
    pivot_longer(!c(sample), names_to = "species", values_to = "avg_n")
  
  return(fish_L_long)
}

boot_L_long <- lapply(boot_L, pivot_boot_L)

alpha_div <- list()

for (i in seq_along(boot_L_long)) {
  
  alpha_div[[i]] <- boot_L_long[[i]] %>% 
    separate_wider_delim(sample, delim = "_", names = c("site", "ipa"), cols_remove = TRUE) %>% 
    group_by(site) %>% 
    summarize(richness = specnumber(avg_n),
              shannon = diversity(avg_n, index = "shannon"),
              simpson = diversity(avg_n, index = "simpson"),
              invsimpson = diversity(avg_n, index = "invsimpson"),
              sum_avg_n = sum(avg_n)) %>% 
    ungroup() 
}

alpha_div_df <- alpha_div[[1]]

for (i in 2:length(alpha_div)) {
  alpha_div_df <- rbind(alpha_div_df, alpha_div[[i]])
}

alpha_div_df <- alpha_div_df %>% 
  mutate(region = ifelse(site %in% c("FAM", "TUR", "COR"), "North", "South"))

plot_site_TD_index <- function (index){
  ggplot(data = alpha_div_df, aes(x = .data[[index]], 
                                  y = factor(site, levels = rev(SOS_core_sites)), 
                                  fill = region, color = region)) +
    geom_density_ridges(alpha = 0.9) + 
    theme_classic() +
    theme(axis.title.y = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
}

td_index_plots <- lapply(names(alpha_div_df[2:6]), plot_site_TD_index)

td_index_names <- c("Species Richness", "Shannon", "Simpson", "InvSimpson", "Sum of Average Count")

for (i in 1:length(td_index_plots)) {
  td_index_plots[[i]] <- td_index_plots[[i]] + xlab(td_index_names[i])
}

td_index_plots[[1]]  + td_index_plots[[2]]  + guide_area() + td_index_plots[[3]] + td_index_plots[[4]] + td_index_plots[[5]] + 
  plot_layout(ncol = 3, guides = "collect")


index_plots[[5]]  + index_plots[[1]]  + guide_area() + index_plots[[2]] + index_plots[[3]] + index_plots[[4]] + 
  td_index_plots[[2]] +  td_index_plots[[3]] +  td_index_plots[[4]] +
  plot_layout(ncol = 3, guides = "collect")

ggsave("~/docs/figures/all_indices.png")


