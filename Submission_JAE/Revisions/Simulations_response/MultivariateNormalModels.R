### FB 27/06/2019 -- Models of correlation between species, to test values of the Gross index / mean Pearson correlations
### In support of the response letter, not to include in the manuscript. 

### Two case studies so far

library(MASS)
library(Matrix)
library(synchrony)
library(codyn)


nsamples = 35

### 1. 10 pairs with -0.8 correlation, unit variance, no correlation outside pairs

rho = -0.8
SigmaPair = matrix(c(1,rho,rho,1),2,2,byrow=TRUE)

Sigma = bdiag(rep(list(SigmaPair),10))
Sigma
mu = rep(0,20)
x = mvrnorm(n = nsamples, mu, Sigma)

cor(x)
mean(cor(x))
meancorr(x)
## Why the hell can't I get the right result the first time?
synchrony(df = x, metric = "Gross") ## argh, needs to be reformated in a dataframe 
## (easier to recompute from scratch based on custom function?)
## check out source("../../../SCRIPTS/test_synchrony_Gross.r")

### 2. Modified Beta distributed rho entries, with a skew towards positive (skew towards negative, Beta(2,4) does not work)
nspecies = 20
IsPosDef = FALSE
while(IsPosDef==FALSE){
  #SigmaBis = matrix(-1+2*rbeta(n=nspecies*nspecies, 5, 2),nspecies,nspecies)
  #SigmaBis = matrix(-1+2*rbeta(n=nspecies*nspecies, 5, 2),nspecies,nspecies)
  #SigmaBis = matrix(-1+2*rbeta(n=nspecies*nspecies, 3, 3),nspecies,nspecies)
  SigmaBis = matrix(-1+2*rbeta(n=nspecies*nspecies, 15, 15),nspecies,nspecies)
  diag(SigmaBis) = rep(1,nspecies)
  #SigmaBis[lower.tri(SigmaBis)] = SigmaBis[upper.tri(SigmaBis)]
  SigmaBis[upper.tri(SigmaBis)] = t(SigmaBis)[upper.tri(SigmaBis)]
  SigmaBis
  isSymmetric(SigmaBis)
  det(SigmaBis-diag(nspecies))
  IsPosDef = min(Re(eigen(SigmaBis)$values))>0.05
  #IsPosDef = rcond(SigmaBis)>0.05 #0.05# 0.01
  print(IsPosDef)
}
hist(SigmaBis)
### NV
rcond(Hilbert(nspecies))
rcond(SigmaBis)
xbis = mvrnorm(n = nsamples, mu, SigmaBis) # difficult to get a sufficiently pos-def matrix in there!
# Some ideas there that may help quickening this process, maybe to explore
# https://www.r-bloggers.com/fixing-non-positive-definite-correlation-matrices-using-r-2/
meancorr(xbis)