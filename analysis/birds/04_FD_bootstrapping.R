library(here)
library(tidyverse)
library(rsample)
library(mFD)
library(ggridges)
library(patchwork)

set.seed(1993)

#load tidy bird data frame created in 01_tidy_data
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

#bootstrap & keep equal proportions for each month to maintain seasonal variability
boot_obj <- bootstraps(mat_to_boot, strata = month, times = 999)

#store the bootstrap dataframes in a list
boot_list <- list()
for (i in seq_along(boot_obj$splits)) {
  resample <- boot_obj$splits[[i]]
  boot_list[[i]] <- analysis(resample)
  
}

format_bird_L <- function(df) {
  
  bird_L <- df %>% 
    group_by(site, spp_code) %>% 
    summarize(spp_total_avg = mean(spp_sum)) %>% #average across months/years
    ungroup() %>% 
    pivot_wider(names_from = spp_code, values_from = spp_total_avg, values_fill = 0) %>% 
    column_to_rownames(var = "site")
  
  return(bird_L)
}

boot_L <- lapply(boot_list, format_bird_L)

#remove species that aren't represented in assemblages
remove_missing_spp <- function(df) {
  filtered_df <- df %>% 
    select_if(colSums(.) != 0)
  
  return(filtered_df)
}

boot_L_filtered <- lapply(boot_L, remove_missing_spp)

#### prep the functional trait space ####
traits.cat <- data.frame(trait_name = colnames(bird.list$trait),
                         trait_type = c("Q", "N", "N", "N", "N"))

#Species trait summary
traits_summary <- sp.tr.summary(tr_cat = traits.cat, 
                                sp_tr = bird.list$trait, 
                                stop_if_NA = T)

#create a list of trait matrices based on species represented in the bootstrapped L matrices
bird_Q_list <- list() 
dist_mat <- list()
space_quality <- list()
trait_space <- list()

for (i in 1:length(boot_L)) {
  bird_Q_list[[i]] <- bird.list$trait %>% 
    rownames_to_column(var = "species") %>% 
    filter(species %in% colnames(boot_L_filtered[[i]])) %>% 
    column_to_rownames(var = "species")

#create the trait space
dist_mat[[i]] <- funct.dist(sp_tr = bird_Q_list[[i]], 
                       tr_cat = traits.cat,
                       metric = "gower",
                       weight_type = "equal",
                       stop_if_NA = TRUE)

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
    rownames_to_column(var = "n_axes") %>% 
    slice(which.min(mad))
  
  space_qual_df <- rbind(space_qual_df, space_qual)
}

space_qual_df %>% 
  filter(!is.na(n_axes)) %>% 
  group_by(n_axes) %>% 
  summarize(n = n(), min = min(mad), max = max(mad), avg = mean(mad))
  #5 axes is almost always the best trait space but 4 is nearly as good

space_qual_df %>% 
  filter(!is.na(n_axes)) %>% 
  ggplot() + 
  geom_histogram(aes(x = mad), fill = "lightblue", color = "skyblue") +
  xlab("Quality of the Trait Space") +
  theme_classic()

#### calculate the functional indices ####
gowdist.list <- lapply(bird_Q_list, gowdis, ord = "podani")
FD_results <- data.frame()

# with the FD package
for (i in 1:length(boot_L)){
    
  birdFD <- dbFD(x = gowdist.list[[i]], #must be a distance object or df where character columns are factors
                 a = data.matrix(boot_L_filtered[[i]]),
                 corr = "none", #mFD package gives an explanation of why sqrt is misleading, and just removing the negative eigenvalues is preferred
                 m = 5,
                 calc.FDiv = TRUE,
                 print.pco = FALSE)
  
  FD_values <- cbind(birdFD$nbsp, birdFD$FRic, birdFD$FEve, birdFD$FDiv, birdFD$FDis) #extract indices
  
  FD_results_df <- FD_values %>%
    as_tibble(rownames = "site")
  
  FD_results <- rbind(FD_results, FD_results_df)

}

colnames(FD_results)[2:6] <- c("Species_Richness", "FRic", "FEve", "FDiv", "FDis")
index_names <- c("Species Richness", "F. Richness","F. Eveness","F. Divergence","F. Dispersion", )

FD_results <- FD_results %>% 
  mutate(site = factor(site, levels = SOS_core_sites)) %>% 
  mutate(region = ifelse(site %in% c("FAM", "TUR", "COR"), "North", "South"))

# FD_bird_boot_results <- FD_results
# save(FD_bird_boot_results, file = "data/FD_bird_boot_results.Rda")

plot_site_index <- function (index){
  ggplot(data = FD_results, aes(x = .data[[index]], 
                                y = factor(site, levels = rev(SOS_core_sites)), 
                                fill = region, color = region)) +
    geom_density_ridges(alpha = 0.9) + 
    theme_classic() +
    theme(axis.title.y = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
}

index_plots <- lapply(names(FD_results[3:6]), plot_site_index)

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

index_plots[[5]]  + index_plots[[1]]  + guide_area() + index_plots[[4]] + index_plots[[3]] + index_plots[[2]] + 
  plot_layout(ncol = 3, guides = "collect")

# ggsave("docs/figures/bird_FDbootpatch.png")

#### test for equal variances ####

kruskal.test(Species_Richness ~ site, data = FD_results) #different
kruskal.test(Species_Richness ~ region, data = FD_results) #different

kruskal.test(FRic ~ site, data = FD_results) #different
kruskal.test(FRic ~ region, data = FD_results) #different

kruskal.test(FEve ~ site, data = FD_results) #different
kruskal.test(FEve ~ region, data = FD_results) #different

kruskal.test(FDiv ~ site, data = FD_results) #different
kruskal.test(FDiv ~ region, data = FD_results) #different

kruskal.test(FDis ~ site, data = FD_results) #different
kruskal.test(FDis ~ region, data = FD_results) #different

#### test for significance ####
SR_model <- lm(Species_Richness ~ site, data = FD_results)

plot(SR_model)

confint(SR_model, level = 0.95)
