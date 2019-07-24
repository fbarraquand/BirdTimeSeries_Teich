### FB 27/06/2019 -- Models of correlation between species, to test values of the Gross index / mean Pearson correlations
### In support of the response letter, not to include in the manuscript. 
### CP 24/07/2019 -- Modifying the code so that we can use Gross index

graphics.off()
rm(list=ls())

### Two case studies so far

library(MASS)
library(Matrix)
library(synchrony)
library(codyn)

source("../../../SCRIPTS/iaaft.R")
source("../../../SCRIPTS/test_synchrony_Gross.r")


nsamples = 35
anrands=0
amethod="shift"
type_correct="BH"
nspecies_list = 40
nrep=100

set.seed(42)

pdf("res_simulations.pdf")
par(mfrow=c(2,2))

### 1. 10 pairs with -0.8 correlation, unit variance, no correlation outside pairs
seq_sp=seq(2,nspecies_list,2)

rho = -0.8
rhomean_list=matrix(NA,ncol=length(seq_sp),nrow=nrep)
eta_list=matrix(NA,ncol=length(seq_sp),nrow=nrep)

sp=0

for(nspecies in seq_sp){
sp=sp+1

SigmaPair = matrix(c(1,rho,rho,1),2,2,byrow=TRUE)

Sigma = bdiag(rep(list(SigmaPair),nspecies/2))
Sigma
mu = rep(0,nspecies)

for(r in 1:nrep){
x = mvrnorm(n = nsamples, mu, Sigma)

#cor(x)

#mean(cor(x)) #CP These two  are different, am wondering why.
rhomean_list[r,sp]=meancorr(x)$obs
## Why the hell can't I get the right result the first time?

#Let's turn the matrix into a data.frame
name_species=paste("Sp",1:nspecies,sep="")
sp_data_frame=c()
abundance=c()
dates=c()
for(i in 1:nspecies){
	sp_data_frame=c(sp_data_frame,rep(name_species[i],nsamples))
	abundance=c(abundance,x[,i])
	dates=c(dates,1:nsamples)
}
data_x=data.frame(dates,sp_data_frame,abundance)

#synchrony(df = data_x, metric = "Gross") ## argh, needs to be reformated in a dataframe 
## (easier to recompute from scratch based on custom function?) CP: Yes, if we wanted to use only the Gross index, but maybe we also want to check if the test works well
## check out source("../../../SCRIPTS/test_synchrony_Gross.r")

#Trying Gross
synch_x=community_sync_Gross(data_x,nrands=anrands,method=amethod)
eta_list[r,sp]=synch_x$obs
}
}

par(mar=c(.5,4,3.5,1))
eta_mean=apply(eta_list,2,mean)
eta_sd=apply(eta_list,2,sd)
eta_min=eta_mean-eta_sd
eta_max=eta_mean+eta_sd
plot( seq_sp,eta_mean,t="p",xlab="",ylab="",pch=16,xaxt="n",ylim=c(-0.85,0.5),main="eta")
 arrows(seq_sp,eta_min,seq_sp,eta_max,angle=90,length=0.1,code=3)
mtext("Paircorr",side=2,las=3,line=3,font=2)

par(mar=c(.5,4,3.5,1))
rho_mean=apply(rhomean_list,2,mean)
rho_sd=apply(rhomean_list,2,sd)
rho_min=rho_mean-rho_sd
rho_max=rho_mean+rho_sd
plot( seq_sp,rho_mean,t="o",xlab="",ylab="",pch=16,xaxt="n",ylim=c(-0.85,0.5),main="rhomean")
arrows(seq_sp,rho_min,seq_sp,rho_max,angle=90,length=0.1,code=3)

### 2. Modified Beta distributed rho entries, with a skew towards positive (skew towards negative, Beta(2,4) does not work)

rhomean_list=matrix(NA,ncol=length(seq_sp),nrow=nrep)
eta_list=matrix(NA,ncol=length(seq_sp),nrow=nrep)

sp=0

for(nspecies in seq_sp){
sp=sp+1
for(r in 1:nrep){
IsPosDef = FALSE
iter=0

################ FB first version
#CP Let's use brute force algorithm
if(1==0){
while(IsPosDef==FALSE&iter<5000){
	iter=iter+1
  #SigmaBis = matrix(-1+2*rbeta(n=nspecies*nspecies, 5, 2),nspecies,nspecies)
  #SigmaBis = matrix(-1+2*rbeta(n=nspecies*nspecies, 5, 2),nspecies,nspecies)
  #SigmaBis = matrix(-1+2*rbeta(n=nspecies*nspecies, 3, 3),nspecies,nspecies)
  SigmaBis = matrix(-1+2*rbeta(n=nspecies*nspecies, 15, 15),nspecies,nspecies)
  diag(SigmaBis) = rep(1,nspecies)
  #SigmaBis[lower.tri(SigmaBis)] = SigmaBis[upper.tri(SigmaBis)]
  SigmaBis[upper.tri(SigmaBis)] = t(SigmaBis)[upper.tri(SigmaBis)]
#  SigmaBis
#  isSymmetric(SigmaBis)
#  det(SigmaBis-diag(nspecies))
  IsPosDef = min(Re(eigen(SigmaBis)$values))>10^(-15) #Was 0.05 before
  #IsPosDef = rcond(SigmaBis)>0.05 #0.05# 0.01
}
}
############ Up to 18 species with max_iter=500, and min(eigenvalue>=10^-15

##FB & CP 2nd version
#Trying # https://www.r-bloggers.com/fixing-non-positive-definite-correlation-matrices-using-r-2/
if(1==0){
  SigmaBis = matrix(-1+2*rbeta(n=nspecies*nspecies, 15, 15),nspecies,nspecies)
  diag(SigmaBis) = rep(1,nspecies)
  SigmaBis[upper.tri(SigmaBis)] = t(SigmaBis)[upper.tri(SigmaBis)]
  IsPosDef = min(Re(eigen(SigmaBis)$values))>0.05

  newSigmaBis=SigmaBis

while(IsPosDef==FALSE&iter<5000){
	iter=iter+1

    # replace -ve eigen values with small +ve number
    newEig <- eigen(newSigmaBis)
    newEig2 <- ifelse(newEig$values < 0, 0, newEig$values)

    # create modified matrix eqn 5 from Brissette et al 2007, inv = transp for
    # eig vectors
    newSigmaBis <- newEig$vectors %*% diag(newEig2) %*% t(newEig$vectors)

    # normalize modified matrix eqn 6 from Brissette et al 2007
    newSigmaBis <- newSigmaBis/sqrt(diag(newSigmaBis) %*% t(diag(newSigmaBis)))

	  IsPosDef = min(Re(eigen(newSigmaBis)$values))>10^(-10) #was 0.05
}
}
#Not better: up to 10 species when using 5000 iter and min(eigenvalue)>10^-10
#Funny: here, 10^(-15) does not lead to a positive SigmaBis according to mvrnorm

  SigmaBis = matrix(-1+2*rbeta(n=nspecies*nspecies, 15, 15),nspecies,nspecies)
  diag(SigmaBis) = rep(1,nspecies)
  SigmaBis[upper.tri(SigmaBis)] = t(SigmaBis)[upper.tri(SigmaBis)]
  IsPosDef = min(Re(eigen(SigmaBis)$values))>0.05

if(!IsPosDef){
SigmaBis <- as.matrix(nearPD(SigmaBis, corr=T, do2eigen=T,posd.tol=10^(-10))$mat)
          IsPosDef = min(Re(eigen(SigmaBis)$values))>10^(-10) #was 0.05
}

if(IsPosDef){
#hist(SigmaBis)
### NV
#rcond(Hilbert(nspecies))
#rcond(SigmaBis)
mu = rep(0,nspecies)
xbis = mvrnorm(n = nsamples, mu, SigmaBis) # difficult to get a sufficiently pos-def matrix in there!
rhomean_list[r,sp]=meancorr(xbis)$obs


# Some ideas there that may help quickening this process, maybe to explore
# https://www.r-bloggers.com/fixing-non-positive-definite-correlation-matrices-using-r-2/


#Let's turn the matrix into a data.frame
name_species=paste("Sp",1:nspecies,sep="")
sp_data_frame=c()
abundance=c()
dates=c()
for(i in 1:nspecies){
        sp_data_frame=c(sp_data_frame,rep(name_species[i],nsamples))
        abundance=c(abundance,xbis[,i])
        dates=c(dates,1:nsamples)
}
data_xbis=data.frame(dates,sp_data_frame,abundance)

#Trying Gross
synch_x=community_sync_Gross(data_xbis,nrands=anrands,method=amethod)
eta_list[r,sp]=synch_x$obs
} #end if(IsPosDef)
} #end r
}

par(mar=c(4,4,1,1))
eta_mean=apply(eta_list,2,mean)
eta_sd=apply(eta_list,2,sd)
eta_min=eta_mean-eta_sd
eta_max=eta_mean+eta_sd
plot( seq_sp,eta_mean,t="p",xlab="nspecies",ylab="",pch=16,ylim=c(-0.85,0.5),main="eta")
 arrows(seq_sp,eta_min,seq_sp,eta_max,angle=90,length=0.1,code=3)
mtext("Beta",side=2,las=3,line=3,font=2)

par(mar=c(4,4,1,1))
rho_mean=apply(rhomean_list,2,mean)
rho_sd=apply(rhomean_list,2,sd)
rho_min=rho_mean-rho_sd
rho_max=rho_mean+rho_sd
plot( seq_sp,rho_mean,t="o",xlab="nspecies",ylab="",pch=16,ylim=c(-0.85,0.5),main="rhomean")
 arrows(seq_sp,rho_min,seq_sp,rho_max,angle=90,length=0.1,code=3)

dev.off()
