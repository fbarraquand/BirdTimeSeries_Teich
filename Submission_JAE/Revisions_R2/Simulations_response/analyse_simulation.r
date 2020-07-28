### 2018 CPicoche: This script computes the wavelet-based synchrony index (Keitt 2008) for the wader community
### 2019/07/03: CP Relaunched with BH-correction 
### 2019/07/25: CP Relaunched adding Calidris 
### 2019/09/04: CP Modified to take into account biomasses and increased a lot the number of species to ignore
### 2020/07/01: CP normalize values

graphics.off()
rm(list=ls())
library('mvcwt')
#source("../../../SCRIPTS/image_mvcwt.r") #Add to change the image function to have a nice Color Bar
source("image_mvcwt_for_colormaps.r") #Add to change the image function to have a nice Color Bar
library("lubridate")
library("RColorBrewer")

set.seed(42)
norm=c(FALSE,TRUE)

#end_of_file_seq=c("4sp_pair","40sp_pair","4sp_alpha15_beta15","40sp_alpha15_beta15","4sp_alpha2_beta4","40sp_alpha2_beta4","4sp_alpha4_beta2","40sp_alpha4_beta2")
#explain=c("rho=-0.8","rho=-0.8","quasi-normal","quasi-normal","compensation","compensation","synchrony","synchrony")
end_of_file_seq=c("4sp_alpha2_beta4","4sp_alpha4_beta2")
explain=c("compensation","synchrony")

for(e in 1:length(end_of_file_seq)){
#for(e in 1:1){
end_of_file=end_of_file_seq[e]
#############################
tab_bm=read.table(paste("MockData_SAD_",end_of_file,".csv",sep=""),sep=";",dec=".",header=T)

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
mm=mvcwt(x,tab_species_tmp,min.scale=0.2,max.scale=10.0)

#This function computes the wavelet ratio of the whole community (see Keitt's paper in 2008)
mr = wmr.boot(mm, smoothing = 1,reps=50)

#png('OUT/Figure3.png',width=800)
pdf(paste("wavelet_simu_",end_nor,"_",end_of_file,"TESTforcontourlines.pdf",sep=""))
  image_mvcwt_for_colormaps(mr,reset.par=F,cex.axis=4,z.fun="Mod",main=explain[e])

#abline(v=2006,lwd=3,col="black") #This is supposed to change in 2006 with water management
print(Sys.time())
dev.off()
}
}
