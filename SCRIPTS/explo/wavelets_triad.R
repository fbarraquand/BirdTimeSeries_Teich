### 2018 This code allows us to plot the wavelet figure for Cormorant/Egret/Heron

graphics.off()
rm(list=ls())
library('mvcwt')
source("SCRIPTS/image_mvcwt.r") #Add to change the image function to have a nice Color Bar
library("RColorBrewer")

db=read.csv("IN/summed_abundances.csv",sep=";",header=T)
db$Date=as.Date(db$Date)

db=subset(db,Nom_latin %in% c("Cormorant","HeronEgret"))
db_tmp=db
names(db_tmp)=c("sp_data_frame","dates","abundance")
db_av=subset(db_tmp,dates<as.Date("2007-01-01"))
db_ap=subset(db_tmp,dates>=as.Date("2007-01-01"))
db_tmp$dates=as.numeric(db_tmp$dates)
db_av$dates=as.numeric(db_av$dates)
db_ap$dates=as.numeric(db_ap$dates)

#synch_all=community_sync_Gross(db_tmp,nrands=100) #-0.09735, p=0.74
#synch_av=community_sync_Gross(db_av,nrands=100) #-0.2028886, p=0.34
#synch_ap=community_sync_Gross(db_ap,nrands=100) #0.2758426, p=0.11


dates=unique(db$Date)
tab=matrix(0,nrow=length(dates),ncol=2)
colnames(tab)=c("Cormorant","HeronEgret")
rownames(tab)=dates
for(id in 1:length(dates)){
        for (s in c('Cormorant',"HeronEgret")){
                id_d=which(db$Date==dates[id]&db$Nom_latin==s)
                if(length(id_d)>0){
                        tab[id,s]=db$Nombre[id_d]
                }
        }
}

x=(dates-dates[1])/365.25
year_min=1981
#This function computes the Morlet wavelet transform for each bird species separately
print(Sys.time())
mm=mvcwt(x,tab,min.scale=0.2,max.scale=10.0)

#This function computes the wavelet ratio of the whole community (see Keitt's paper in 2008)
mr = wmr.boot(mm, smoothing = 1,reps=100)
mr$x=mr$x+year_min #Change the dates to be "human-readable"

png('Submission_JAE/Revisions/wavelet_triad_BH.png',width=800) #Was BY before, image_mvcwt has been changed on 2019/07/03

image_mvcwt(mr,reset.par=F,cex.axis=4,z.fun="Mod")

#abline(v=2006,lwd=3,col="black") #This is supposed to change in 2006 with water management
print(Sys.time())
dev.off()

