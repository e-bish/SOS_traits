#adapted from Gomez et al. 2017  https://doi.org/10.1111/2041-210X.12856


library(rjags)

# Begin JAGS model specification to fit the model
# This section writes a text file that jags will read

#### The Beta n-mixture model for jags
sink("beta_nmixture.txt")
cat("model{

		# Priors
		# Priors for detectability
		a	~	dunif(0,10)
		b	~	dunif(0,10)
		
		# Priors for abundance
		for(i in 1:nspp){
			lambda[i]	~	dgamma(0.01,0.01)
		}
		
		# likelihood
		
		# Loop along clones		
		
		for(k in 1:K){
		
			# Estimate abundance for each speices
		
			for (i in 1:nspp){
		
				# Sample detection probability from distribution
		
				p[i,k] ~ dbeta(a,b)
				
				# Loop along sites
				for (j in 1:nsites){
					
					# Abundance porcess
					
					N[i,j,k] ~ dpois(lambda[i])
					
					for(t in 1:nvisits){
					
						# Detection process
						counts[j,t,i,k]~dbin(p[i,k],N[i,j,k])
				
					}
				}
			}
		}	
	}
")
sink()


# Define parameters required. 
# The following uses parallel computation to run the simulation and assumes that the computer has
# 20 cores to run parallel. To run the code appropriately find the number of cores in the computer and 
# adjust nreps and star.index accordingly. 

#### load count data here ####

nsites	<-	25
nvisits	<-	3
area.50m	<-	pi*0.50^2
alpha		<-	0.25*4.5
beta		<-	0.75*4.5
abundances	<-	c(1,2,3,4,5,7,10,15,25,40,55,65,75,85,100)
nreps		<-	20
K			<-	20
full.res	<-	list()
index	<-	1	# index is taken as an outargument from the batch/shell file used to run the program

start.index	<-	seq(0,500,by=20)
tmp.spp.counts.ls	<-	loP.counts[start.index[index]:(start.index[index]+nreps)]

# Estimating the abundance of the species using Beta mixture model
# This uses data cloning with 20 clons.

K		<-	20

spp.counts	<-	UIF.counts
spp.clons		<-	array(NA,dim=c(dim(spp.counts),K))
nspp			<-	dim(spp.counts)[3]

for (j in 1:K){
  
  spp.clons[,,,j]	<-	spp.counts
}

max.count		<-	apply(spp.clons,c(3,1,4),max)	

data_list		<-	list(counts=spp.clons,nspp=dim(spp.counts)[3],nsites=dim(spp.counts)[1],nvisits=dim(spp.counts)[2],K=K)	

inits			<-	function(){list(lambda=runif(nspp,0.1,10),a=2,b=2,N=max.count)}

params		<-	c("lambda","a","b")

jm			<-	jags.model(file="beta_nmixture.txt"
                   ,data=data_list,n.chains=2,n.adapt=1000,inits=inits)
out1			<-	coda.samples(model=jm,variable.names=params,n.iter=20000,thin=20)

hats			<-	summary(out1)$statistics[,1]

prop_sigma	<-	1.96*(sqrt(summary(out1)$statistics[,2]^2*20))
upper		<-	(hats+prop_sigma)*100/0.78
lower		<-	(hats-prop_sigma)*100/0.78

out1.sum		<-	summary(out1)$statistics

