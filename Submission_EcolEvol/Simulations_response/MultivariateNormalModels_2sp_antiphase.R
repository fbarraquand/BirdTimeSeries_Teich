graphics.off()
rm(list=ls())

### Two case studies so far

library(MASS)
library(Matrix)
library(synchrony)
library(codyn)
library(mvcwt)

source("../SCRIPTS/iaaft.R")
source("../SCRIPTS/test_synchrony_Gross.r")
source("../SCRIPTS/image_mvcwt_for_colormaps.r") 
library("RColorBrewer")

nsamples = 35 *12 #We want to look at monthly data
anrands=1000
amethod="shift"
type_correct="BH"
nrep=1 #Should be 100
mu_min_coeff=0.1
mu_max_coeff=0.9
mu=c(200,500)
c=0.5
rho = 1
doyouload=T
type_dist=c("Compensation")
nspecies=c(2)
norm=F
set.seed(42)



### 2. Modified Beta distributed rho entries, with a skew towards positive (skew towards negative, Beta(2,4) does not work)

for(r in 1:nrep){
x=matrix(NA,nsamples,nspecies)
for(t in 1:nsamples){
#Compute mu_t
mu_t=mu_min_coeff*mu+mu*(mu_max_coeff-mu_min_coeff)*(1+0.5*sin(2*pi*t/12+c(0,pi))) #We want to antiphase populations
sigma_i=c*mu_t
sigma_i[2]=4*sigma_i[1]
#mu_t=mu_min_coeff*mu+mu*(mu_max_coeff-mu_min_coeff)*sin(2*pi*t/12)
SigmaPair = matrix(c(sigma_i[1]^2,0,0,sigma_i[2]^2),2,2,byrow=TRUE)

x[t,] = mvrnorm(n = 1, mu_t, SigmaPair)
}
#cor(x)

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
write.table(data_tot,paste("../../../Teich_resultsLFS/simulated_timeseries_very_skewed_SAD/MockData_SAD_2sp_antiphase.csv",sep=""),sep=";",dec=".",row.names=F)

tab=cbind(mu,SigmaPair)
colnames(tab)=c("mu",paste("Sp",1:nspecies))
write.table(tab,paste("../../../Teich_resultsLFS/simulated_timeseries_very_skewed_SAD/MuSigma_SAD_2sp_antiphase.csv",sep=""),sep=";",dec=".",row.names=F,append=T)
} #end r

pdf("MockData_cormorantheronegret_timeseries.pdf",width=20,height=6)
plot(x[,2],col="grey",ylim=range(x[,2]),t="o",pch=16)
lines(x[,1],col="black",t="o",pch=16)
dev.off()

end_of_file_seq="2sp_antiphase"
explain=c("antiphase")

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


pdf(paste("Skewed_SAD_",end_of_file,"_IAAFT.pdf",sep=""),width=7,height=3)
layout(matrix(c(1,1,2),nrow=1,ncol=3,byrow=T),widths=c(6,2))
print(paste(Sys.time(),"before image"))
par(mar=c(3,5,2,3))
image_mvcwt_for_colormaps(ref_wmr,reset.par=F,cex.axis=4,z.fun="Mod",adj="None")
#mtext("b)",side=2,line=-2,at=0.48,cex=1.5,outer=T,las=1)
print(paste(Sys.time(),"after image"))

#abline(v=2006,lwd=3,col="black") #This is supposed to change in 2006 with water management
print("After wavelet")
print(Sys.time())
dev.off()

}
}
