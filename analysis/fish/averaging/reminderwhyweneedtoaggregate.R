#load tidy fish data frame created in 02_tidy_data
load(here("data", "net_tidy.Rdata")) 

#### Create the L (abundance) matrix ####
# expand_species <- net_tidy %>% 
#   expand(nesting(year, month, site, ipa), ComName) %>% 
#   filter(!is.na(ComName))

fish_L <- net_tidy %>% #L is referring to the RLQ analysis
  filter(!is.na(ComName)) %>% 
  group_by(year, month, site, ipa, ComName) %>%
  summarize(spp_sum = sum(species_count)) %>% #sum across depths within a site
  ungroup() %>%
  # full_join(expand_species) %>% #add back in all of the events so we can capture 0s
  # arrange(year, month, site, ipa) %>% 
  # mutate(spp_sum =replace_na(spp_sum, 0)) %>% 
  # group_by(site, ipa, ComName) %>% 
  # summarize(spp_avg = mean(spp_sum)) %>% #average across sampling events at each site/ipa (unbalanced)
  pivot_wider(names_from = ComName, values_from = spp_sum, values_fill = 0) %>% 
  clean_names() %>% 
  ungroup() %>% 
  mutate(sample = paste(year, month, site, ipa, sep = "_"), .after = ipa) %>% 
  select(!1:4) %>% 
  column_to_rownames(var = "sample")

fish_Q <- fish.list$trait.t


#examine the quality of the potential functional spaces using the mFD package
dist_mat <- gowdis(fish_Q, ord = "classic")
space_quality <- quality.fspaces(sp_dist = dist_mat,
                                 maxdim_pcoa = 10,
                                 deviation_weighting = "absolute",
                                 fdist_scaling = FALSE,
                                 fdendro = "ward.D2")

round(space_quality$"quality_fspaces",3) #lowest value is the best (<.1 is good), meaning species pairs are accurately represented
#aka the distances in euclidean space are accurately reflecting the gowers distances

n_axes_to_retain <- 4

#plot the trait space with the first 4 axes
trait_space <- space_quality$"details_fspaces"$"sp_pc_coord"

alpha_indices <- alpha.fd.multidim(sp_faxes_coord = trait_space[ , 1:n_axes_to_retain],
                                   asb_sp_w = data.matrix(fish_L),
                                   ind_vect = c("fdis", "feve", "fric", "fdiv"),
                                   scaling = TRUE,
                                   check_input = TRUE,
                                   details_returned = TRUE)
