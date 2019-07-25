### 2018 CPicoche: This script computes the wavelet-based synchrony index (Keitt 2008) for the wader community
### 2019/07/03: CP Relaunched with BH-correction 
### 2019/07/25: CP Relaunched adding Calidris 

graphics.off()
rm(list=ls())
library('mvcwt')
source("SCRIPTS/image_mvcwt.r") #Add to change the image function to have a nice Color Bar
library("lubridate")
library("RColorBrewer")

##################################################################################
# -- importation des données du Teich
#DBt<-read.csv(file="/home/cpicoche/Documents/Data_to_treat/TRANSFERT_LIMICOLES/IN/DBWithMonthlyPhotoTeich.csv",header=TRUE,sep=",",dec=".")
DBt<-read.csv(file="/home/cpicoche/Documents/Birds/BirdTimeSeries_Teich/IN/DBWithMonthlyPhotoTeich_completed.csv",header=TRUE,sep=",",dec=".")
DBt = subset(DBt,(DBt$Protocol==1 | DBt$Protocol==2) & (DBt$Lieu_dit=="USN00-Réserve ornithologique (générique)" | DBt$Lieu_dit=="USN01-Artigues-Réserve ornithologique")  )#& !(DBt$Annee %in% c(2007,2008)))
DBt$Date=as.Date(as.character(DBt$Date))
minAnnee = as.numeric(format(min(DBt$Date), format = "%Y"))
maxAnnee = as.numeric(format(max(DBt$Date), format = "%Y"))

DBt<-DBt[!colnames(DBt) %in% c("X","Ref")] #edit du 14/04/2017 : vu qu'il y a 3 doublons
DBt = unique.matrix(DBt) # élimination de cette manière, mais il faudrait aller regarder dans la photo mensuelle 
# les doublons sont présents dans le fichier d'origine, 
# pour les virer avec unique.matrix je suis obligée de sortir les deux premières colonnes

#-----------------------------------------------------------------------------------------
#Combien d'espèces? 
limicoles = sort(c("Recurvirostra avosetta","Limosa limosa","Limosa lapponica","Calidris temminckii","Calidris canutus",
              "Calidris alba","Calidris alpina","Calidris minuta","Calidris maritima" ,"Gallinago gallinago",
              "Tringa flavipes","Tringa nebularia","Tringa erythropus","Tringa ochropus","Tringa totanus",
              "Tringa glareola","Actitis hypoleucos","Philomachus pugnax","Numenius arquata","Numenius phaeopus",
              "Himantopus himantopus","Charadrius hiaticula","Charadrius alexandrinus","Haematopus ostralegus",
              "Burrhinus oedicnemus","Charadrius dubius","Phalaropus lobatus","Pluvialis squatarola",
              "Pluvialis apricaria","Arenaria interpres","Vanellus vanellus","Calidris ferruginea"))
sp_to_ignore=c("Anas discors","Anas americana","Calidris melanotos","Calidris pusilla","Calidris ruficollis", "Calidris fuscicollis", "Calidris himantopus", "Burhinus oedicnemus","Phalaropus lobatus","Charadrius alexandrinus","Haematopus ostralegus","Calidris maritima")
limicoles=limicoles[!(limicoles %in% sp_to_ignore)]

#############################


#Warning : I'm cutting the beginning of the series here, because there are no waders before 1980
year_min=1980
DBt=subset(DBt,year(Date)>year_min)

dates=sort(unique(DBt$Date))

tab_limicoles=matrix(0,nrow=length(dates),ncol=length(limicoles))
colnames(tab_limicoles)=limicoles

for(id in 1:length(dates)){
	for (s in limicoles){
		id_d=which(DBt$Date==dates[id]&DBt$Nom_latin==s)
		if(length(id_d)>0){
			tab_limicoles[id,s]=DBt$Nombre[id_d]
		}
	}
}

x=(dates-dates[1])/365.25

#This function computes the Morlet wavelet transform for each bird species separately
print(Sys.time())
mm=mvcwt(x,tab_limicoles,min.scale=0.2,max.scale=10.0)

#This function computes the wavelet ratio of the whole community (see Keitt's paper in 2008)
mr = wmr.boot(mm, smoothing = 1,reps=100)
mr$x=mr$x+year_min #Change the dates to be "human-readable"

#png('OUT/Figure3.png',width=800)
png("Submission_JAE/Revisions/Figure3_BH.png",width=800)
image_mvcwt(mr,reset.par=F,cex.axis=4,z.fun="Mod")

#abline(v=2006,lwd=3,col="black") #This is supposed to change in 2006 with water management
print(Sys.time())
dev.off()
