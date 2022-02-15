graphics.off()
rm(list=ls())

### Two case studies so far

library(MASS)
library(Matrix)
library(synchrony)
library(codyn)

source("../../SCRIPTS/iaaft.R")
source("../../SCRIPTS/test_synchrony_Gross.r")
source("../../SCRIPTS/image_mvcwt_for_colormaps.r") 

nsamples = 35 *12 #We want to look at monthly data
anrands=1000
amethod="shift"
type_correct="BH"
nrep=1 #Should be 100
mu_min_coeff=0.1
mu_max_coeff=0.9
m=3.26
c=0.5
rho = -0.8
doyouload=T
alpha=c(2)
beta=c(4)
type_dist=c("Compensation")
seq_sp=c(4)
norm=F
set.seed(42)




### 2. Modified Beta distributed rho entries, with a skew towards positive (skew towards negative, Beta(2,4) does not work)
rhomean_list=matrix(NA,ncol=length(seq_sp),nrow=nrep)
eta_list=matrix(NA,ncol=length(seq_sp),nrow=nrep)


for(t_sigma in 1:length(alpha)){
#for(t_sigma in 2:2){
sp=0
print(paste('alpha',alpha[t_sigma],"beta",beta[t_sigma]))
for(nspecies in seq_sp){

print(nspecies)
sp=sp+1
#Define mu
log_mu=rnorm(nspecies,m)
mu=exp(log_mu)
id_most_abundant=which(mu==max(mu))
coef_multiple=rep(1,length(mu))
if(nspecies==40){
coef_multiple[id_most_abundant]=100
}else if(nspecies==4){
coef_multiple[id_most_abundant]=7
}else{
}

for(r in 1:nrep){
x=matrix(NA,nsamples,nspecies)
for(t in 1:nsamples){
#Compute mu_t

mu_t=mu_min_coeff*mu+coef_multiple*mu*(mu_max_coeff-mu_min_coeff)*(1+0.5*sin(2*pi*t/12+pi)) #+pi so that it decreases at the beginning of the year
#mu_t=mu_min_coeff*mu+mu*(mu_max_coeff-mu_min_coeff)*sin(2*pi*t/12)
sigma_i=c*mu_t

IsPosDef = FALSE
iter=0
while(!IsPosDef){
  SigmaPair = matrix(-1+2*rbeta(n=nspecies*nspecies, alpha[t_sigma], beta[t_sigma]),nspecies,nspecies)

for(i in 1:nspecies){
        for(j in 1:nspecies){
                        SigmaPair[i,j]=SigmaPair[i,j]*sigma_i[i]*sigma_i[j]
        }
}
  SigmaPair[upper.tri(SigmaPair)] = t(SigmaPair)[upper.tri(SigmaPair)]


IsPosDef = min(Re(eigen(SigmaPair)$values))>0.05

if(!IsPosDef){
try(SigmaPair <- as.matrix(nearPD(SigmaPair, corr=F, do2eigen=T,posd.tol=10^(-5))$mat),silent=T)
          IsPosDef = min(Re(eigen(SigmaPair)$values))>10^(-10) #was 0.05
}
}
x[t,] = mvrnorm(n = 1, mu_t, SigmaPair)
}
#cor(x)

l=apply(x,2,function(x) diff(range(x)))
print("Ratio amplitude x/sum(others)")
l_bis=l[-id_most_abundant]
print(l[id_most_abundant]/(sum(l_bis)))
#mean(cor(x)) #CP These two  are different, am wondering why.
rhomean_list[r,sp]=meancorr(x)$obs

#pdf("MockData_SAD_timeseries.pdf",width=20,height=6)
plot(1:nsamples,x[,2],col="grey",t="o",pch=16,ylim=range(c(x)),xlab="Time",ylab="Abundance")
lines(x[,1],col="black",t="o",pch=16)
lines(x[,3],col="red",t="o",pch=16)
lines(x[,4],col="blue",t="o",pch=16)
dev.off()

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

data_tot=cbind(rep(r,nrow(data_x)),data_x)
colnames(data_tot)=c("Rep","Time","Species","Abundance")
#write.table(data_tot,paste("../../../Teich_resultsLFS/simulated_timeseries_very_skewed_SAD/MockData_SAD_",nspecies,"sp_alpha",alpha[t_sigma],"_beta",beta[t_sigma],".csv",sep=""),sep=";",dec=".",row.names=F)

tab=cbind(mu,SigmaPair)
colnames(tab)=c("mu",paste("Sp",1:nspecies))
#write.table(tab,paste("../../../Teich_resultsLFS/simulated_timeseries_very_skewed_SAD/MuSigma_SAD_",nspecies,"sp_alpha",alpha[t_sigma],"_beta",beta[t_sigma],".csv",sep=""),sep=";",dec=".",row.names=F,append=T)
} #end r
} #end nspecies
} #end t_sigma

end_of_file_seq=paste(nspecies,"sp_alpha",alpha[t_sigma],"_beta",beta[t_sigma],sep="")
explain=c("compensation")

for(e in 1:length(end_of_file_seq)){
#for(e in 1:1){
end_of_file=end_of_file_seq[e]
#############################
tab_bm=read.table(paste("../../../Teich_resultsLFS/simulated_timeseries_very_skewed_SAD/MockData_SAD_",end_of_file,".csv",sep=""),sep=";",dec=".",header=T)

rep=1

tab_bm=subset(tab_bm,Rep==rep)

tab_bm$Species=as.character(tab_bm$Species)
tab_bm$Time=as.numeric(as.character(tab_bm$Time))
tab_bm$Abundance=as.numeric(as.character(tab_bm$Abundance))

dates_tmp=unique(tab_bm$Time)

#Convert 1,2,3 to month and year

dates=rep(NA,length(dates_tmp))
yy=0
for(d in 1:length(dates)){
        if(dates_tmp[d]%%12==1){
                yy=yy+1
        }
        mm=dates_tmp[d]%%12
        if(mm==0){mm=12}
        tmp=as.Date(paste(yy,mm,"15",sep="-"))
        dates[d]=tmp
}

sp=unique(tab_bm$Species)

tab_species=matrix(0,nrow=length(dates),ncol=length(sp))
colnames(tab_species)=sp
for(id in 1:length(dates)){
        for (s in sp){
                id_d=which(tab_bm$Time==dates_tmp[id]&tab_bm$Species==s)
                if(length(id_d)>0){
                        tab_species[id,s]=tab_bm$Abundance[id_d]
                }
        }
}

for(normalize in norm){
if(normalize){
end_nor="scaled"
}else{
end_nor="NOTscaled"
}

tab_species_tmp=tab_species
for(s in sp){
if(normalize){
        tab_species_tmp[,s]=scale(tab_species[,s])
}

}

x=(dates-dates[1])/365.25

#This function computes the Morlet wavelet transform for each bird species separately
print(Sys.time())
seq_x=seq(1,365.25*35,length.out=420)/365.25

mm=mvcwt(seq_x,tab_species_tmp,min.scale=mean(diff(seq_x))*3,max.scale=10.0,nscales=100,loc=seq_x) ####WARNING. THIS ONLY WORKS BECAUSE x IS ARTIFICIAL)

print(paste(Sys.time(),"after mvcwt"))

#This function computes the wavelet ratio of the whole community (see Keitt's paper in 2008)
if(!doyouload){

ref_wmr=wmr(mm)
ref_val=ref_wmr$z[,,1]

tab_values_iaaft=array(NA,dim=c(length(mm$x),length(mm$y),anrands+1))
tab_values_iaaft[,,anrands+1]=ref_val
prog.bar = txtProgressBar(min = 0, max = anrands+1,style = 3)
for(i in 1:anrands){
        setTxtProgressBar(prog.bar, i)
        tab_tmp=tab_species_tmp
        for(s in 1:ncol(tab_species)){
        tab_tmp[,s]=iaaft_surrogate(tab_species_tmp[,s])
        }
	mmtmp=mvcwt(x,tab_tmp,min.scale=mean(diff(seq_x))*3,max.scale=10.0,nscales=100,loc=seq_x)
        wmr_tmp=wmr(mmtmp)
        tab_values_iaaft[,,i]=wmr_tmp$z[,,1]
}

tab_pval=array(NA,dim=c(length(mm$x),length(mm$y),1))
for(i in 1:length(mm$x)){
        for(j in 1:length(mm$y)){
#                tab_pval[i,j,1]= 2*min(sum(tab_values_iaaft[i,j,] >= ref_val[i,j]),sum(tab_values_iaaft[i,j,] < ref_val[i,j]))/(anrands+1)
                tab_pval[i,j,1]= sum(tab_values_iaaft[i,j,] <= ref_val[i,j])/(anrands+1)
                if(tab_pval[i,j,1]>1){stop()}

        }
}
ref_wmr$z.boot=tab_pval

if(length(ref_wmr$x)>length(ref_wmr$y)){
	yy=c(ref_wmr$y,rep(NA,length(ref_wmr$x)-length(ref_wmr$y)))
	xx=ref_wmr$x
}else{
	xx=c(ref_wmr$x,rep(NA,length(ref_wmr$y)-length(ref_wmr$x)))
	yy=ref_wmr$y
}
tab_xy=cbind(xx,yy)
colnames(tab_xy)=c("x","y")
write.table(tab_xy,paste("../../../Teich_resultsLFS/simulated_timeseries_very_skewed_SAD/tab_xy_mr_simulated_data_",end_nor,"_with",anrands,"_",end_of_file,"_IAAFT.csv",sep=""),sep=";",dec=".",col.names=T,row.names=F)

tab_z=ref_wmr$z
write.table(as.matrix(tab_z[,,1]),paste("../../../Teich_resultsLFS/simulated_timeseries_very_skewed_SAD/tab_z_mr_simulated_data_",end_nor,"_with",anrands,"_",end_of_file,"_IAAFT.csv",sep=""),sep=";",dec=".",col.names=F,row.names=F)

tab_z.boot=ref_wmr$z.boot
write.table(as.matrix(tab_z.boot[,,1]),paste("../../../Teich_resultsLFS/simulated_timeseries_very_skewed_SAD/tab_zboot_mr_simulated_data_",end_nor,"_with",anrands,"_",end_of_file,"_IAAFT.csv",sep=""),sep=";",dec=".",col.names=F,row.names=F)

}else{
ref_wmr=wmr(mm)

tmp_xy=read.csv(paste("../../../Teich_resultsLFS/simulated_timeseries_very_skewed_SAD/tab_xy_mr_simulated_data_",end_nor,"_with",anrands,"_",end_of_file,"_IAAFT.csv",sep=""),header=T,sep=";",dec=".")
ref_wmr$x=tmp_xy[!is.na(tmp_xy[,"x"]),"x"]
ref_wmr$y=tmp_xy[!is.na(tmp_xy[,"y"]),"y"]


tmp_z=as.matrix(read.csv(paste("../../../Teich_resultsLFS/simulated_timeseries_very_skewed_SAD/tab_z_mr_simulated_data_",end_nor,"_with",anrands,"_",end_of_file,"_IAAFT.csv",sep=""),header=F,sep=";",dec="."))
tmp_array_z=array(0,dim=c(dim(tmp_z),1))
tmp_array_z[,,1]=tmp_z
ref_wmr$z=tmp_array_z

tmp_z.boot=as.matrix(read.csv(paste("../../../Teich_resultsLFS/simulated_timeseries_very_skewed_SAD/tab_zboot_mr_simulated_data_",end_nor,"_with",anrands,"_",end_of_file,"_IAAFT.csv",sep=""),header=F,sep=";",dec="."))
tmp_array_z.boot=array(0,dim=c(dim(tmp_z.boot),1))
tmp_array_z.boot[,,1]=tmp_z.boot
ref_wmr$z.boot=tmp_array_z.boot
}


pdf(paste("Skewed_SAD_",nspecies,"sp_compensation_IAAFT.pdf",sep=""),width=7,height=4)
layout(matrix(c(1,1,2),nrow=1,ncol=3,byrow=T),widths=c(6,2))
print(paste(Sys.time(),"before image"))
par(mar=c(5,5,2,3))
image_mvcwt_for_colormaps(ref_wmr,reset.par=F,cex.axis=4,z.fun="Mod",adj="None")
mtext("Scale (years)", side=2, line=-1.75,at=0.55, outer = TRUE,cex=1.1)
mtext("Years",side=1, line=-2, at=0.475,outer = TRUE,cex=1.1)
print(paste(Sys.time(),"after image"))

#abline(v=2006,lwd=3,col="black") #This is supposed to change in 2006 with water management
print("After wavelet")
print(Sys.time())
dev.off()

}
}
