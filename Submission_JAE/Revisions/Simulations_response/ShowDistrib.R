#### CP 24/07/2019
#### Checking the difference between the distributions we want in the correlation matrix (sometimes not positive definite), and the correlation we finally have, both with nearPD and with the observed x
graphics.off()
rm(list=ls())

library(Matrix)
library(synchrony)

set.seed(42)

seq_nsp=c(10,20,30,40)
nsamples = 35

pdf("compare_distrib.pdf")
par(mfrow=c(4,3))
for(nspecies in seq_nsp){
#Quasi-normal distribution
  SigmaBis = matrix(-1+2*rbeta(n=nspecies*nspecies, 15, 15),nspecies,nspecies)
  diag(SigmaBis) = rep(1,nspecies)
  SigmaBis[upper.tri(SigmaBis)] = t(SigmaBis)[upper.tri(SigmaBis)]
   if(nspecies==40){
	ax="Input"
   }else{
	ax=""
   }
  hist(SigmaBis,breaks=10,main="",ylab=paste(nspecies,"spp"),xlab=ax,xlim=c(-1,1))

   SigmaBis <- as.matrix(nearPD(SigmaBis, corr=T, do2eigen=T,posd.tol=10^(-10))$mat)
   if(nspecies==10){
	amain="Quasi-normal"
   }else{
	amain=""
   }
   if(nspecies==40){
	ax="near pos. def."
   }else{
	ax=""
   }
   hist(SigmaBis,breaks=10,main=amain,ylab="",xlab=ax,xlim=c(-1,1))

   mu = rep(0,nspecies)
   xbis = mvrnorm(n = nsamples, mu, SigmaBis) # difficult to get a sufficiently pos-def matrix in there!
   if(nspecies==40){
	ax="Output"
   }else{
	ax=""
   }
   hist(cor(xbis),breaks=10,main="",ylab="",xlab=ax,xlim=c(-1,1))
}
#Compensation
for(nspecies in seq_nsp){
  SigmaBis = matrix(-1+2*rbeta(n=nspecies*nspecies, 2, 4),nspecies,nspecies)
  diag(SigmaBis) = rep(1,nspecies)
  SigmaBis[upper.tri(SigmaBis)] = t(SigmaBis)[upper.tri(SigmaBis)]
   if(nspecies==40){
	ax="Input"
   }else{
	ax=""
   }
  hist(SigmaBis,breaks=10,main="",ylab=paste(nspecies,"spp"),xlab=ax,xlim=c(-1,1))

   SigmaBis <- as.matrix(nearPD(SigmaBis, corr=T, do2eigen=T,posd.tol=10^(-10))$mat)
   if(nspecies==10){
        amain="Compensation"
   }else{
        amain=""
   }
   if(nspecies==40){
	ax="near pos. def."
   }else{
	ax=""
   }
   hist(SigmaBis,breaks=10,main=amain,ylab="",xlab=ax,xlim=c(-1,1))

   mu = rep(0,nspecies)
   xbis = mvrnorm(n = nsamples, mu, SigmaBis) # difficult to get a sufficiently pos-def matrix in there!
   if(nspecies==40){
	ax="Output"
   }else{
	ax=""
   }
   hist(cor(xbis),breaks=10,main="",ylab="",xlab=ax,xlim=c(-1,1))
}

#Synchrony
for(nspecies in seq_nsp){
SigmaBis = matrix(-1+2*rbeta(n=nspecies*nspecies, 4, 2),nspecies,nspecies)
  diag(SigmaBis) = rep(1,nspecies)
  SigmaBis[upper.tri(SigmaBis)] = t(SigmaBis)[upper.tri(SigmaBis)]
   if(nspecies==40){
	ax="Input"
   }else{
	ax=""
   }
  hist(SigmaBis,breaks=10,main="",ylab=paste(nspecies,"spp"),xlab=ax,xlim=c(-1,1))

   SigmaBis <- as.matrix(nearPD(SigmaBis, corr=T, do2eigen=T,posd.tol=10^(-10))$mat)
   if(nspecies==10){
        amain="Synchrony"
   }else{
        amain=""
   }
   if(nspecies==40){
	ax="near pos. def."
   }else{
	ax=""
   }
   hist(SigmaBis,breaks=10,main=amain,ylab="",xlab=ax,xlim=c(-1,1))

   mu = rep(0,nspecies)
   xbis = mvrnorm(n = nsamples, mu, SigmaBis) # difficult to get a sufficiently pos-def matrix in there!
   if(nspecies==40){
	ax="Output"
   }else{
	ax=""
   }
   hist(cor(xbis),breaks=10,main="",ylab="",xlab=ax,xlim=c(-1,1))
}

dev.off()
