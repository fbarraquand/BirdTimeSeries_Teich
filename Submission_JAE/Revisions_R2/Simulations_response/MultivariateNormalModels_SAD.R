### FB 27/06/2019 -- Models of correlation between species, to test values of the Gross index / mean Pearson correlations
### In support of the response letter, not to include in the manuscript. 
### CP 24/07/2019 -- Modifying the code so that we can use Gross index
### CP 02/07/2020 -- Using a different SAD for birds, as well as monthly abundances

graphics.off()
rm(list=ls())

### Two case studies so far

library(MASS)
library(Matrix)
library(synchrony)
library(codyn)

source("../../../SCRIPTS/iaaft.R")
source("../../../SCRIPTS/test_synchrony_Gross.r")


nsamples = 35 *12 #We want to look at monthly data
anrands=0
amethod="shift"
type_correct="BH"
nspecies_list = 40 #Starting out with a random number of species just to check that it works
nrep=2 #Should be 100
mu_min_coeff=0.1
mu_max_coeff=0.9
m=3.26
c=0.5

alpha=c(15,2,4)
beta=c(15,4,2)
type_dist=c("Pair","Quasi-normal","Compensation","Synchrony")

set.seed(42)


#pdf("res_simulations_different_distributions_SAD.pdf")
par(mfrow=c(4,2))

### 1. 10 pairs with -0.8 correlation, unit variance, no correlation outside pairs
#seq_sp=seq(2,nspecies_list,2)


seq_sp=c(4,40) ###REMOVE

#if(1==0){
rho = -0.8
rhomean_list=matrix(NA,ncol=length(seq_sp),nrow=nrep)
eta_list=matrix(NA,ncol=length(seq_sp),nrow=nrep)

sp=0

for(nspecies in seq_sp){
sp=sp+1
#Define mu
log_mu=rnorm(nspecies,m)
mu=exp(log_mu)


for(r in 1:nrep){
x=matrix(NA,nsamples,nspecies)
for(t in 1:nsamples){
#Compute mu_t
mu_t=mu_min_coeff*mu+mu*(mu_max_coeff-mu_min_coeff)*sin(2*pi*t/12)
sigma_i=c*mu_t

SigmaPair = matrix(rho,nrow=nspecies,ncol=nspecies)
diag(SigmaPair)=1

for(i in 1:nspecies){
	for(j in 1:nspecies){
		if(j==i){
			SigmaPair[i,j]=1
		}else{
			SigmaPair[i,j]=SigmaPair[i,j]*sigma_i[i]*sigma_i[j]
		}
	}
}

IsPosDef = min(Re(eigen(SigmaPair)$values))>0.05

if(!IsPosDef){
SigmaPair <- as.matrix(nearPD(SigmaPair, corr=T, do2eigen=T,posd.tol=10^(-5))$mat)
          IsPosDef = min(Re(eigen(SigmaPair)$values))>10^(-10) #was 0.05
}


x[t,] = mvrnorm(n = 1, mu_t, SigmaPair)

}
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

data_tot=cbind(rep(r,nrow(data_x)),data_x)
colnames(data_tot)=c("Rep","Time","Species","Abundance")
if(nspecies==4|nspecies==40){
write.table(data_tot,paste("MockData_SAD_",nspecies,"sp_pair.csv",sep=""),sep=";",dec=".",row.names=F,append=T)
}
} #end loop on nrep
}#end loop on nspecies_seq

par(mar=c(3,5,2.5,1))
eta_mean=apply(eta_list,2,mean)
eta_sd=apply(eta_list,2,sd)
eta_min=eta_mean-eta_sd
eta_max=eta_mean+eta_sd
#plot( seq_sp,eta_mean,t="p",xlab="",ylab="",pch=16,xaxt="n",ylim=c(-0.85,0.65),main="eta")
plot( seq_sp,eta_mean,t="p",xlab="",ylab="",pch=16,xaxt="n",ylim=c(min(eta_min),max(eta_max)),main=expression(eta))
 arrows(seq_sp,eta_min,seq_sp,eta_max,angle=90,length=0.1,code=3)
mtext("Paircorr",side=2,las=3,line=3,font=2)
abline(h=0.0,col="red")

par(mar=c(3,5,2.5,1))
rho_mean=apply(rhomean_list,2,mean)
rho_sd=apply(rhomean_list,2,sd)
rho_min=rho_mean-rho_sd
rho_max=rho_mean+rho_sd
#plot( seq_sp,rho_mean,t="o",xlab="",ylab="",pch=16,xaxt="n",ylim=c(-0.85,0.65),main="rhomean")
plot( seq_sp,rho_mean,t="o",xlab="",ylab="",pch=16,xaxt="n",ylim=c(min(rho_min),max(rho_max)),main=expression(bar(rho)))
arrows(seq_sp,rho_min,seq_sp,rho_max,angle=90,length=0.1,code=3)
abline(h=0.0,col="red")


#}#end 1==0
### 2. Modified Beta distributed rho entries, with a skew towards positive (skew towards negative, Beta(2,4) does not work)

rhomean_list=matrix(NA,ncol=length(seq_sp),nrow=nrep)
eta_list=matrix(NA,ncol=length(seq_sp),nrow=nrep)


for(t_sigma in 1:length(alpha)){
sp=0
print(paste('alpha',alpha[t_sigma],"beta",beta[t_sigma]))
for(nspecies in seq_sp){
print(nspecies)
sp=sp+1
#Define mu
log_mu=rnorm(nspecies,m)
mu=exp(log_mu)


for(r in 1:nrep){
x=matrix(NA,nsamples,nspecies)
for(t in 1:nsamples){
#Compute mu_t
mu_t=mu_min_coeff*mu+mu*(mu_max_coeff-mu_min_coeff)*sin(2*pi*t/12)
sigma_i=c*mu_t

IsPosDef = FALSE
iter=0

  SigmaPair = matrix(-1+2*rbeta(n=nspecies*nspecies, alpha[t_sigma], beta[t_sigma]),nspecies,nspecies)
  diag(SigmaPair) = rep(1,nspecies)

for(i in 1:nspecies){
        for(j in 1:nspecies){
                if(j==i){
                        SigmaPair[i,j]=1
                }else{
                        SigmaPair[i,j]=SigmaPair[i,j]*sigma_i[i]*sigma_i[j]
                }
        }
}
  SigmaPair[upper.tri(SigmaPair)] = t(SigmaPair)[upper.tri(SigmaPair)]


IsPosDef = min(Re(eigen(SigmaPair)$values))>0.05

if(!IsPosDef){
SigmaPair <- as.matrix(nearPD(SigmaPair, corr=T, do2eigen=T,posd.tol=10^(-5))$mat)
          IsPosDef = min(Re(eigen(SigmaPair)$values))>10^(-10) #was 0.05
}


x[t,] = mvrnorm(n = 1, mu_t, SigmaPair)

}
#cor(x)

#mean(cor(x)) #CP These two  are different, am wondering why.
rhomean_list[r,sp]=meancorr(x)$obs

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

#Trying Gross
synch_x=community_sync_Gross(data_x,nrands=anrands,method=amethod)
eta_list[r,sp]=synch_x$obs

data_tot=cbind(rep(r,nrow(data_x)),data_x)
colnames(data_tot)=c("Rep","Time","Species","Abundance")
if(nspecies==4|nspecies==40){
write.table(data_tot,paste("MockData_SAD_",nspecies,"sp_alpha",alpha[t_sigma],"_beta",beta[t_sigma],".csv",sep=""),sep=";",dec=".",row.names=F,append=T)
}
} #end r
} #end nspecies

eta_mean=apply(eta_list,2,mean)
eta_sd=apply(eta_list,2,sd)
eta_min=eta_mean-eta_sd
eta_max=eta_mean+eta_sd
if(t_sigma==length(alpha)){
par(mar=c(3,5,2.5,1))
ax=""
#plot( seq_sp,eta_mean,t="p",xlab=ax,ylab="",pch=16,ylim=c(-0.85,0.65))
plot( seq_sp,eta_mean,t="p",xlab=ax,ylab="",pch=16,ylim=c(min(eta_min),max(eta_max)))
mtext("Nb species",1,2,cex=0.75)
}else{
par(mar=c(3,5,2.5,1))
ax=""
#plot( seq_sp,eta_mean,t="p",xlab=ax,ylab="",pch=16,ylim=c(-0.85,0.65),xaxt="n")
plot( seq_sp,eta_mean,t="p",xlab=ax,ylab="",pch=16,ylim=c(min(eta_min),max(eta_max)),xaxt="n")
}
 arrows(seq_sp,eta_min,seq_sp,eta_max,angle=90,length=0.1,code=3)
mtext(type_dist[t_sigma+1],side=2,las=3,line=3,font=2)
abline(h=0.0,col="red")

rho_mean=apply(rhomean_list,2,mean)
rho_sd=apply(rhomean_list,2,sd)
rho_min=rho_mean-rho_sd
rho_max=rho_mean+rho_sd
if(t_sigma==length(alpha)){
#plot( seq_sp,rho_mean,t="o",xlab=ax,ylab="",pch=16,ylim=c(-0.85,0.65))
plot( seq_sp,rho_mean,t="o",xlab=ax,ylab="",pch=16,ylim=c(min(rho_min),max(rho_max)))
mtext("Nb species",1,2,cex=0.75)
}else{
#plot( seq_sp,rho_mean,t="o",xlab=ax,ylab="",pch=16,ylim=c(-0.85,0.65),xaxt="n")
plot( seq_sp,rho_mean,t="o",xlab=ax,ylab="",pch=16,ylim=c(min(rho_min),max(rho_max)),xaxt="n")
}
 arrows(seq_sp,rho_min,seq_sp,rho_max,angle=90,length=0.1,code=3)
abline(h=0.0,col="red")
}

dev.off()
