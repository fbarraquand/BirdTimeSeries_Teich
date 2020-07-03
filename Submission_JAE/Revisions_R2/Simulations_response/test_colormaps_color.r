####  CP 03/07/2020 Just trying new colormaps for the wavelet stuff

source("image_mvcwt_for_colormaps.r") #Add to change the image function to have a nice Color Bar
graphics.off()
if(1==0){

rm(list=ls())

library('mvcwt')
library("lubridate")
library("RColorBrewer")


#############################
tab_bm=read.table("MockData_SAD_4sp_alpha2_beta4.csv",sep=";",dec=".",header=T)

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

x=(dates-dates[1])/365.25

#This function computes the Morlet wavelet transform for each bird species separately
print(Sys.time())
mm=mvcwt(x,tab_species,min.scale=0.2,max.scale=10.0)

#This function computes the wavelet ratio of the whole community (see Keitt's paper in 2008)
mr = wmr.boot(mm, smoothing = 1,reps=2)

}
#png('OUT/Figure3.png',width=800)
pdf(paste("wavelet_simu_4sp_alpha15_beta15_test_colors.pdf",sep=""),width=10,height=15)
layout(matrix(c(1,2,3,4,5,6),nrow=3,ncol=2,byrow=T),widths=c(10,2))

  image_mvcwt_for_colormaps(mr,reset.par=F,cex.axis=4,z.fun="Mod",col.palette="Spectral")
  image_mvcwt_for_colormaps(mr,reset.par=F,cex.axis=4,z.fun="Mod",col.palette="Viridis")
  image_mvcwt_for_colormaps(mr,reset.par=F,cex.axis=4,z.fun="Mod",col.palette="Inferno",inv=F)

#abline(v=2006,lwd=3,col="black") #This is supposed to change in 2006 with water management
print(Sys.time())
dev.off()

pdf(paste("wavelet_simu_4sp_alpha15_beta15_test_colors_2.pdf",sep=""),width=10,height=15)
layout(matrix(c(1,2,3,4,5,6),nrow=3,ncol=2,byrow=T),widths=c(10,2))

  image_mvcwt_for_colormaps(mr,reset.par=F,cex.axis=4,z.fun="Mod",col.palette="Magma")
  image_mvcwt_for_colormaps(mr,reset.par=F,cex.axis=4,z.fun="Mod",col.palette="Plasma")
  image_mvcwt_for_colormaps(mr,reset.par=F,cex.axis=4,z.fun="Mod",col.palette="Oranges",inv=T)

#abline(v=2006,lwd=3,col="black") #This is supposed to change in 2006 with water management
print(Sys.time())
dev.off()

pdf(paste("wavelet_simu_4sp_alpha15_beta15_test_colors_3.pdf",sep=""),width=10,height=15)
layout(matrix(c(1,2,3,4,5,6),nrow=3,ncol=2,byrow=T),widths=c(10,2))

  image_mvcwt_for_colormaps(mr,reset.par=F,cex.axis=4,z.fun="Mod",col.palette="PuBu")
  image_mvcwt_for_colormaps(mr,reset.par=F,cex.axis=4,z.fun="Mod",col.palette="YlGnBu")
  image_mvcwt_for_colormaps(mr,reset.par=F,cex.axis=4,z.fun="Mod",col.palette="BrBG")

#abline(v=2006,lwd=3,col="black") #This is supposed to change in 2006 with water management
print(Sys.time())
dev.off()

