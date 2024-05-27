library(tidyverse)
library(here)
library(ade4)
library(mvabund)

#all
load(here("data", "fish.list.Rdata")) #object created in 02_tidy_data
load(here("data", "env_table.Rdata")) #object created in spatial analysis

#just June 
# load(here("data", "fish.list.june.Rdata")) #object created in 02_tidy_data
# load(here("data", "env.data.june.Rdata")) #object created in 02_tidy_data

trait <- fish.list$trait %>% 
  as.data.frame() %>% 
  mutate(mean_length_mm = as.numeric(mean_length_mm)) %>% 
  mutate_if(is.character, as.factor)

abu <- sqrt(fish.list$abund) %>% #use the corrected abundance matrix??
  as_data_frame() #needs to be a df for the fourthcorner function
#correcting doesnt seem to hugely change the result

env <- env.data

#COA on abundance matrix L (site x species)
#aka what species correspond with what sites
### Should this be on the corrected abundance matrix? 
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
fit <- traitglm(abu,env,trait,method="glm1path")
## why does lasso give an error?!
fit$fourth

plot(fit)

anova(fit, nBoot = 10)
# p = 0.91

a        = max( abs(fit$fourth.corner) )
colort   = colorRampPalette(c("blue","white","red")) 
plot.4th = levelplot(t(as.matrix(fit$fourth.corner)), xlab="Environmental Variables",
                     ylab="Species traits", col.regions=colort(100), at=seq(-a, a, length=100),
                     scales = list( x= list(rot = 45)))
print(plot.4th)

# another method
nrepet <- 49999
four.comb <- fourthcorner(env, abu,
                                trait, modeltype = 6, p.adjust.method.G = "none",
                                p.adjust.method.D = "none", nrepet = nrepet)

plot(four.comb, alpha = 0.05, stat = "D2") 
# in D2 sthe association is measured between the quantitative variable and each category separately

#now adjust the p-values
four.comb.adj <- p.adjust.4thcorner(four.comb,
                                    p.adjust.method.G = "fdr", p.adjust.method.D = "fdr")
plot(four.comb.adj, alpha = 0.05, stat = "D2")
#nothing is significant lol