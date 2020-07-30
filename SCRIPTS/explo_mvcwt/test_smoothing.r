# Script by CPicoche

rm(list=ls())
graphics.off()
source("../test_synchrony_Gross.r")
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
#This function computes the Morlet wavelet transform for each bird species separately
mm=mvcwt(x,tab,min.scale=0.2,max.scale=10.0)

print(paste("Smooth2",Sys.time()))
distrib2=wmr.boot(mm, smoothing = 2,reps=100)
pdf("smooth2.pdf")
layout(matrix(c(1,2),nrow=1,ncol=2,byrow=T),widths=c(10,2))
image_mvcwt_for_colormaps(distrib2,reset.par=F,cex.axis=4,z.fun="Mod")
dev.off()
stop()


print(paste("Smooth 1",Sys.time()))
distrib1=wmr.boot(mm, smoothing = 1,reps=100)
print(paste("Smooth 5",Sys.time()))
distrib5=wmr.boot(mm, smoothing = 5,reps=100)
print(paste("Smooth 10",Sys.time()))
distrib10=wmr.boot(mm, smoothing = 10,reps=100)
print(paste("Smooth 100",Sys.time()))
distrib100=wmr.boot(mm, smoothing = 100,reps=100)

pdf("smooth1.pdf")
layout(matrix(c(1,2),nrow=1,ncol=2,byrow=T),widths=c(10,2))
image_mvcwt_for_colormaps(distrib1,reset.par=F,cex.axis=4,z.fun="Mod")
dev.off()

pdf("smooth5.pdf")
layout(matrix(c(1,2),nrow=1,ncol=2,byrow=T),widths=c(10,2))
image_mvcwt_for_colormaps(distrib5,reset.par=F,cex.axis=4,z.fun="Mod")
dev.off()

pdf("smooth10.pdf")
layout(matrix(c(1,2),nrow=1,ncol=2,byrow=T),widths=c(10,2))
image_mvcwt_for_colormaps(distrib10,reset.par=F,cex.axis=4,z.fun="Mod")
dev.off()

pdf("smooth100.pdf")
layout(matrix(c(1,2),nrow=1,ncol=2,byrow=T),widths=c(10,2))
image_mvcwt_for_colormaps(distrib100,reset.par=F,cex.axis=4,z.fun="Mod")
dev.off()
