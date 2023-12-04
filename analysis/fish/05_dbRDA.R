
fish_dbrda.all <- dbrda(fish.list$abund ~ FRic + FEve + FDis + Rao,
                        data = FD_results,
                        distance = "bray",
                        sqrt.dist = TRUE,
                        add = FALSE,
                        dfun = vegdist,
                        metaMDSdist = FALSE,
                        na.action = na.exclude,
                        subset = NULL
)

dbrda_scores <- scores(fish_dbrda.all)
sites_scores <- as.data.frame(dbrda_scores[[1]])
biplot_scores <- as.data.frame(dbrda_scores[[2]])

dbRDA_ordiplot <- gg_ordiplot(ord = sites_scores, #for some reason the scale gets weird if you don't specify this
                              groups = FD_results$shoreline,
                              ellipse = TRUE,
                              hull = FALSE,
                              spiders = FALSE)
