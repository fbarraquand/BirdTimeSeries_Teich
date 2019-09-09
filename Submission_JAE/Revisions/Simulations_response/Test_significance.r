#CP 25/07/2019 This script is meant to test the significance of the eta value for different distributions

graphics.off()
rm(list=ls())

library(Matrix)
library(synchrony)
library(MASS)

source("../../../SCRIPTS/iaaft.R")
source("../../../SCRIPTS/test_synchrony_Gross.r")

nsamples = 35
anrands=100
amethod="shift"
nspecies_list = c(2,6,10,20,30,40)
thresh=0.1

alpha=c(15,2,4)
beta=c(15,4,2)
type_dist=c("Pair","Quasi-normal","Compensation","Synchrony")

set.seed(42)

pdf("test_for_significance.pdf")
par(mfrow=c(2,2),mar=c(4,4,3,2))

#####Pairwise
rho = -0.8
sp=0
plot(0,0,t="n",ylim=c(-1,1),xlim=c(1,length(nspecies_list)),xlab="",ylab=expression(eta),main=type_dist[1],xaxt="n")
for(nspecies in nspecies_list){
sp=sp+1

SigmaPair = matrix(c(1,rho,rho,1),2,2,byrow=TRUE)

Sigma = bdiag(rep(list(SigmaPair),nspecies/2))
mu = rep(0,nspecies)

x = mvrnorm(n = nsamples, mu, Sigma)
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
if(nspecies==2){
synch_x=community_sync_Gross(data_x,nrands=anrands,method="iaaft")
}else{
synch_x=community_sync_Gross(data_x,nrands=anrands,method=amethod)
}
obs=synch_x$obs
p_s=synch_x$pval
points(sp,obs,pch=16,cex=1)
boxplot(synch_x$rands[1:anrands],at=sp,add=T,boxwex=0.25,range=0,yaxt="n",xaxt="n")
if (p_s<=thresh){
points(sp,as.numeric(obs)+0.1,pch='*',col="red",cex=2)
                }
}


####Beta


for(t_sigma in 1:length(alpha)){
sp=0
if(t_sigma>1){
	xl="Nb species"
}else{
	xl=""
}
if(t_sigma==2){
	yl=expression(eta)
}else{
	yl=""
}
plot(0,0,t="n",ylim=c(-1,1),xlim=c(1,length(nspecies_list)),xlab=xl,ylab=yl,main=type_dist[t_sigma+1],xaxt="n")
if(t_sigma>1){
	axis(1,at=1:length(nspecies_list),nspecies_list)
}
for(nspecies in nspecies_list){
sp=sp+1
IsPosDef = FALSE
  SigmaBis = matrix(-1+2*rbeta(n=nspecies*nspecies, alpha[t_sigma], beta[t_sigma]),nspecies,nspecies)
  diag(SigmaBis) = rep(1,nspecies)
  SigmaBis[upper.tri(SigmaBis)] = t(SigmaBis)[upper.tri(SigmaBis)]
  IsPosDef = min(Re(eigen(SigmaBis)$values))>0.05

if(!IsPosDef){
SigmaBis <- as.matrix(nearPD(SigmaBis, corr=T, do2eigen=T,posd.tol=10^(-10))$mat)
          IsPosDef = min(Re(eigen(SigmaBis)$values))>10^(-10) #was 0.05
}

if(IsPosDef){
mu = rep(0,nspecies)
x = mvrnorm(n = nsamples, mu, SigmaBis)
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
if(nspecies==2){
synch_x=community_sync_Gross(data_x,nrands=anrands,method="iaaft")
}else{
synch_x=community_sync_Gross(data_x,nrands=anrands,method=amethod)
}
obs=synch_x$obs
p_s=synch_x$pval
points(sp,obs,pch=16,cex=1)
boxplot(synch_x$rands[1:anrands],at=sp,add=T,boxwex=0.25,range=0,yaxt="n",xaxt="n")
if (p_s<=thresh){
points(sp,as.numeric(obs)+0.1,pch='*',col="red",cex=2)
                }
}
}
}
dev.off()
