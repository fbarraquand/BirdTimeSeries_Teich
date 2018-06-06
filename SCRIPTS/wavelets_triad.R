#graphics.off()
rm(list=ls())
library('mvcwt')
source("SCRIPTS/image_mvcwt.r") #Add to change the image function to have a nice Color Bar
library("RColorBrewer")

db=read.csv("IN/summed_abundances.csv",sep=";",header=T)
db$Date=as.Date(db$Date)

db=subset(db,Nom_latin %in% c("Cormorant","HeronEgret"))

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

png('OUT/wavelet_Driad.png',width=800)
image_mvcwt(mr,reset.par=F,cex.axis=4,z.fun="Mod")

#abline(v=2006,lwd=3,col="black") #This is supposed to change in 2006 with water management
print(Sys.time())
dev.off()

