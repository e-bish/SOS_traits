library(tidyverse)
library(here)
library(ade4)
library(mvabund)

#all
load(here("data", "fish.list.Rdata")) #object created in 02_tidy_data
load(here("data", "env_table.Rdata")) #object created in spatial analysis

#using presence/absence to summarize by site
abu <- fish.list$abund %>% 
  as_tibble(rownames = "sample") %>% 
  separate_wider_delim(sample, delim = "_", names = c("site", "year"), cols_remove = TRUE) %>% 
  select(!year) %>% 
  group_by(site) %>% 
  summarize_all(sum) %>% 
  column_to_rownames("site") %>% 
  decostand(method = "pa")

trait <- fish.list$trait.t

env <- env_table %>% 
  column_to_rownames("site")

#matrices by region
# abu.r <- fish.list$abund %>% 
#   as_tibble(rownames = "sample") %>% 
#   separate_wider_delim(sample, delim = "_", names = c("year", "month", "site"), cols_remove = FALSE) %>% 
#   select(!c(year, month, sample)) %>% 
#   mutate(site = ifelse(site %in% c("FAM", "TUR", "COR"), "north", "south")) %>% 
#   group_by(site) %>% 
#   summarize_all(sum) %>% 
#   column_to_rownames("site") %>% 
#   decostand(method = "pa")
# 
# env.r <- env_table %>% 
#   mutate(site = ifelse(site %in% c("FAM", "TUR", "COR"), "north", "south")) %>% 
#   group_by(site) %>% 
#   summarize_all(mean) %>% 
#   column_to_rownames("site")

#COA on abundance matrix L (site x species)
#aka what species correspond with what sites
coa.abu <- dudi.coa(abu, scannf = FALSE, nf=2)

#PCA on traits matrix Q (species x traits)
#aka what are the principle components of trait variation
pca.trait <- dudi.hillsmith(trait, scannf = FALSE, 
                      row.w = coa.abu$cw)

#PCA on environment matrix R (sites x environment)
#aka what are the principle components of environmental variation
#use dudi.hillsmith for cat and continuous variables or dudi.pca for just continuous
pca.env <- dudi.pca(env, scannf = FALSE, 
                    row.w = coa.abu$lw)

#Compute the RLQ Analysis
rlqF <- rlq(pca.env, coa.abu, pca.trait, 
            scannf = FALSE)

summary(rlqF)
#low correlation on the first axis (0.2645)
plot(rlqF)

#Plot traits score
t1 <- order(rlqF$c1[,1])
dotchart(rlqF$c1[t1,1], pch=16, 
         labels = names(trait)[t1])
abline(v=0, lty=2)

#species scores
# top 10 species with positive score
rlqF$mQ[order(rlqF$mQ[,1], decreasing = TRUE)[1:10],]

# top 10 species with negative score
rlqF$mQ[order(rlqF$mQ[,1])[1:10],]

#Plot environment score
e1 <- order(rlqF$l1[,1])
dotchart(rlqF$l1[e1,1], pch=16,
         labels = names(env)[e1])
abline(v=0, lty=2)

RLQ1 <- rlqF$lR[,1]

save(RLQ1, file = here("data", "RLQ1.Rdata"))

####### fourth corner analysis
#method adds a LASSO penalty which essentially reduces the number of predictors by 
#getting rid of interactions close to zero
fit <- traitglm(abu,env,trait,method="manyglm")
## why does lasso give an error?!
fit$fourth

plot(fit)

anova(fit, nBoot = 10)
# p = 0.727

a        = max( abs(fit$fourth.corner) )
colort   = colorRampPalette(c("blue","white","red")) 
plot.4th = levelplot(t(as.matrix(fit$fourth.corner)), xlab="Environmental Variables",
                     ylab="Species traits", col.regions=colort(100), at=seq(-a, a, length=100),
                     scales = list( x= list(rot = 45)))
print(plot.4th)

# another method
nrepet <- 49999
four.comb.adj <- p.adjust.4thcorner(four.comb,
                                    p.adjust.method.G = "fdr", 
                                    p.adjust.method.D = "fdr")
plot(four.comb.adj, alpha = 0.05, stat = "D2")
#nothing is significant lol


test_combine <- inner_join(FD_results, env_table, by = "site")
  
#redundancy analysis using euclidean distances between FD metrics
fish_rda <- rda(test_combine[5:9] ~ perc.ag + perc.natural + perc.developed + armor.500m + armor.10km,
                      data = test_combine,
                      sqrt.dist = TRUE,
                      add = FALSE,
                      dfun = vegdist,
                      metaMDSdist = FALSE,
                      na.action = na.fail,
                      subset = NULL)

summary(fish_rda)
plot(fish_rda)

rda_scores <- scores(fish_rda)
sites_scores <- as.data.frame(rda_scores[[1]])
biplot_scores <- as.data.frame(rda_scores[[2]])

RDA_ordiplot <- gg_ordiplot(ord = sites_scores, #for some reason the scale gets weird if you don't specify this
                              # groups = c("Species_Richness", "FRic", "FEve", "FDiv", "FDis"),
                              ellipse = TRUE,
                              hull = FALSE,
                              spiders = FALSE)


