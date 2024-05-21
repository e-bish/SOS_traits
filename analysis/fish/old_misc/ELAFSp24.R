library(tidyverse)
library(here)
library(ggrepel)
library(FD)
library(ggordiplots)
library(patchwork)
library(vegan)
library(ggpubr)
library(ade4)
library(mvabund)

load(here("data", "fish.list.Rdata")) #object created in 02_tidy_data
load(here("data", "env_table.Rdata")) #object created in spatial analysis
here("analysis", "general_functions", "geb12299-sup-0002-si.r") %>% source()
here("analysis", "general_functions", "scree.r") %>% source()

#create species x species matrix
fish.gowdist <- gowdis(fish.list$trait)

#view pcoa
pcoa <- cmdscale(sqrt(fish.gowdist)) 
colnames(pcoa) <- c("pcoa1", "pcoa2")

pcoa %>%
  as_tibble(rownames = "species") %>%
  ggplot(aes(x = pcoa1, y = pcoa2)) +
  geom_point(size = 3) +
  labs(x="PCoA Axis 1 (27.53%)",
       y="PCoA Axis 2 (17.37%)") +
  theme_bw(base_size = 16)

pcoa %>%
  as_tibble(rownames = "species") %>%
  ggplot(aes(x = pcoa1, y = pcoa2)) +
  geom_text_repel(
    label=rownames(pcoa),
    max.time = 1,
    max.overlaps = Inf) +
  geom_point() +
  labs(x="PCoA Axis 1 (27.53%)",
       y="PCoA Axis 2 (17.37%)") +
  theme_bw(base_size = 16)

#### FD indices calculation ####################################################
fishFD <- dbFD(x = fish.list$trait, #must be a df where character columns are factors
               a = fish.list$abund,
               corr = "sqrt", #this is the default option to deal with negative eigenvalues
               m = 4, #seems to keep 4 axis no matter how many you specify?
               calc.FGR = TRUE, 
               clust.type = "ward.D2",
               calc.FDiv = TRUE, #won't return a value if there are only categorical traits
               print.pco = TRUE)
#it will prompt to see if you want to separate by groups (g) and ask how many (7)

important.indices <- cbind(fishFD$nbsp, fishFD$FRic, fishFD$FEve, fishFD$FDiv,
                           fishFD$FDis, fishFD$RaoQ) #extract indices

colnames(important.indices) <- c("NumbSpecies", "FRic", "FDiv", "FEve", "FDis", "Rao")

FD_results <- important.indices %>% 
  as_tibble(rownames = "sample") %>% 
  separate_wider_delim(sample, delim = "_", names = c("year", "month", "site"), cols_remove = FALSE) %>% 
  relocate(sample) %>% 
  mutate(month = factor(month, labels = c("Apr", "May", "Jun", "Jul", "Aug", "Sept"))) %>% 
  mutate(site = factor(site, levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG"))) %>% 
  replace(is.na(.), 0) #can't calculate indices with <3 species

FD_results %>% 
  filter(NumbSpecies >2) %>%
  group_by(site) %>% 
  summarize(across(where(is.numeric), mean))

DOK_df <- fish.list[[2]] %>% 
  as_tibble(rownames = "sample") %>% 
  separate_wider_delim(sample, delim = "_", names = c("year", "month", "site"), cols_remove = FALSE) %>% 
  filter(site == "DOK") %>% 
  select(-c(1:4)) %>% 
  select_if(colSums(.) != 0)

DOK_spp <- colnames(DOK_df)

pcoa %>% 
  as_tibble(rownames = "species") %>%
  filter(species %in% DOK_spp) %>% 
  ggplot(aes(x = pcoa1, y = pcoa2)) +
    geom_text_repel(
      label=DOK_spp,
      max.time = 1,
      max.overlaps = Inf) +
    geom_point() +
    labs(x="PCoA Axis 1 (27.53%)",
       y="PCoA Axis 2 (17.37%)") +
    theme_bw(base_size = 16)

COR_df <- fish.list[[2]] %>% 
  as_tibble(rownames = "sample") %>% 
  separate_wider_delim(sample, delim = "_", names = c("year", "month", "site"), cols_remove = FALSE) %>% 
  filter(site == "COR") %>% 
  select(-c(1:4)) %>% 
  select_if(colSums(.) != 0)

COR_spp <- colnames(COR_df)

pcoa %>% 
  as_tibble(rownames = "species") %>%
  filter(species %in% COR_spp) %>% 
  ggplot(aes(x = pcoa1, y = pcoa2)) +
  geom_text_repel(
    label=COR_spp,
    max.time = 1,
    max.overlaps = Inf) +
  geom_point() +
  labs(x="PCoA Axis 1 (27.53%)",
       y="PCoA Axis 2 (17.37%)") +
  theme_bw(base_size = 16)

plot_index <- function (index){
  
  ggplot(data = FD_results, aes(x = site, y = .data[[index]], fill = site)) +
    geom_boxplot(outlier.shape = NA) + 
    stat_compare_means(label = "p.signif", ref.group = ".all.", hide.ns = TRUE) + #using default kruskal-wallace test
    theme_classic() +
    ylab(index) +
    theme(axis.title.x = element_blank())
}

index_plots <- lapply(names(FD_results[5:8]), plot_index)

index_plots[[1]] + index_plots[[2]] + index_plots[[3]] + index_plots[[4]] + 
  plot_layout(ncol = 4, guides = 'collect')


adonis2(FD_results[5:8] ~ FD_results$site, strata = FD_results$year, method = "euc", permutations = 999)

kruskal.test(NumbSpecies ~ site, data = FD_results)

#### RLQ/Fourth corner ####

#try with just presence absence
abu <- fish.list$abund %>% 
  as_tibble(rownames = "sample") %>% 
  separate_wider_delim(sample, delim = "_", names = c("year", "month", "site"), cols_remove = FALSE) %>% 
  select(!c(year, month, sample)) %>% 
  group_by(site) %>% 
  summarize_all(sum) %>% 
  column_to_rownames("site") %>% 
  decostand(method = "pa")

trait <- fish.list$trait

env <- env_table %>% 
  column_to_rownames("site")

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


#method adds a LASSO penalty which essentially reduces the number of predictors by 
#getting rid of interactions close to zero
fit <- traitglm(abu,env,trait,family = "binomial", method="manyglm")
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

#other technique not working
nrepet <- 49999
four.comb <- fourthcorner(env, abu, trait, modeltype = 6, 
                          p.adjust.method.G = "fdr",
                          p.adjust.method.D = "fdr",
                          nrepet = nrepet)

plot(four.comb, alpha = 0.05, stat = "D2")

