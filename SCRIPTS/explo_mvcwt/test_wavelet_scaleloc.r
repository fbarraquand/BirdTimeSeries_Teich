# Script by CPicoche

rm(list=ls())
graphics.off()
library('mvcwt')
source("wmr_boot.r")
source("../image_mvcwt_for_colormaps.r")
library("RColorBrewer")

set.seed(42)
log_b=F
normalize=F
end_bio="abundances"

db=read.csv(paste("IN/summed_",end_bio,"_v2_wtoutrarespecies.csv",sep=""),sep=";",header=T)
db$Date=as.Date(db$Date)

db=subset(db,Nom_latin %in% c("Cormorant","HeronEgret"))
db_tmp=db
names(db_tmp)=c("sp_data_frame","dates","abundance")
db_av=subset(db_tmp,dates<as.Date("2007-01-01"))
db_ap=subset(db_tmp,dates>=as.Date("2007-01-01"))
db_tmp$dates=as.numeric(db_tmp$dates)
db_av$dates=as.numeric(db_av$dates)
db_ap$dates=as.numeric(db_ap$dates)

dates=unique(db$Date)
tab=matrix(0,nrow=length(dates),ncol=2)
colnames(tab)=c("Cormorant","HeronEgret")
rownames(tab)=dates
for(id in 1:length(dates)){
        for (s in c('Cormorant',"HeronEgret")){
                id_d=which(db$Date==dates[id]&db$Nom_latin==s)
                if(length(id_d)>0){
                        if(log_b){
                                tab[id,s]=log(db$Nombre[id_d]+1)
                        }else{
                                tab[id,s]=db$Nombre[id_d]
                        }
                }
        }
}

if(normalize){
        tab[,"Cormorant"]=scale(tab[,'Cormorant'])
        tab[,"HeronEgret"]=scale(tab[,'HeronEgret'])
}

x=(dates-dates[1])/365.25
year_min=1981
if(1==0){
#This function computes the Morlet wavelet transform for each bird species separately
print(paste("Beginning reference wavelet",Sys.time()))
mm_ref=mvcwt(x,tab,min.scale=0.2,max.scale=10.0,nscales=423,loc=regularize(x))
wmr_ref=wmr.boot(mm_ref, smoothing = 1,reps=100)
print(paste("End reference wavelet",Sys.time()))

#mm_ns=mvcwt(x,tab,min.scale=0.2,max.scale=10.0,nscales=212,loc=regularize(x)) #Half the default nscale
#wmr_ns=wmr(mm_ns)

#mm_loc=mvcwt(x,tab,min.scale=0.2,max.scale=10.0,nscales=423,loc=regularize(x,nsteps=length(x)/2))
#wmr_loc=wmr(mm_loc)

print(paste("Beginning small dim wavelet",Sys.time()))
mm_locns=mvcwt(x,tab,min.scale=0.2,max.scale=10.0,nscales=212,loc=regularize(x,nsteps=length(x)/2))
wmr_locns=wmr.boot(mm_locns,smoothing=1,reps=100)
print(paste("End smalldim wavelet",Sys.time()))

pdf("reference_image.pdf")
layout(matrix(c(1,2),nrow=1,ncol=2,byrow=T),widths=c(10,2))
image_mvcwt_for_colormaps(wmr_ref)
dev.off()
pdf("halfscales_image.pdf")
layout(matrix(c(1,2),nrow=1,ncol=2,byrow=T),widths=c(10,2))
image_mvcwt_for_colormaps(wmr_ns)
dev.off()
pdf("halfloc_image.pdf")
layout(matrix(c(1,2),nrow=1,ncol=2,byrow=T),widths=c(10,2))
image_mvcwt_for_colormaps(wmr_loc)
dev.off()
pdf("halfloc_halfscales_image.pdf")
layout(matrix(c(1,2),nrow=1,ncol=2,byrow=T),widths=c(10,2))
image_mvcwt_for_colormaps(wmr_locns)
dev.off()


print(paste("Beginning 0.25 dim wavelet",Sys.time()))
mm_locns=mvcwt(x,tab,min.scale=0.2,max.scale=10.0,nscales=ceiling(423)/4,loc=regularize(x,nsteps=ceiling(length(x)/4)))
wmr_locns4=wmr.boot(mm_locns,smoothing=1,reps=100)
print(paste("End smalldim wavelet",Sys.time()))

pdf("quarterloc_quarterscales_image.pdf")
layout(matrix(c(1,2),nrow=1,ncol=2,byrow=T),widths=c(10,2))
image_mvcwt_for_colormaps(wmr_locns4)
dev.off()

pdf("quarterloc_quarterscales_image_no_adjustment.pdf")
layout(matrix(c(1,2),nrow=1,ncol=2,byrow=T),widths=c(10,2))
image_mvcwt_for_colormaps(wmr_locns4,adj="None")
dev.off()
}

mm_locns3=mvcwt(x,tab,min.scale=0.2,max.scale=10.0,nscales=100,loc=regularize(x,nsteps=ceiling(length(x)/2)))
wmr_locns3=wmr.boot(mm_locns3,smoothing=1,reps=100)
print(paste("End smalldim wavelet",Sys.time()))

pdf("smallscale_image.pdf")
layout(matrix(c(1,2),nrow=1,ncol=2,byrow=T),widths=c(10,2))
image_mvcwt_for_colormaps(wmr_locns3,adj="None")
dev.off()
