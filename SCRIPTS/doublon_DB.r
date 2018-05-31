#analyse des données après la création de la photo mensuelle.
rm(list=ls())
graphics.off()
library("lattice")
library ("RColorBrewer") #pour la génération automatique de couleur
library("corrplot")
library("lubridate")


##################################################################################
# -- importation des données du Teich
#DBt<-read.csv(file="/home/caluome/Documents/DATA/DATA/le_Teich/DBWithMonthlyPhotoTeich.csv",header=TRUE,sep=",",dec=".")
DBt<-read.csv(file="./IN/DBWithMonthlyPhotoTeich_completed.csv",header=TRUE,sep=",",dec=".")

DBt$Date=as.Date(DBt$Date)
DBt = subset(DBt,(DBt$Protocol==1 | DBt$Protocol==2) & (DBt$Lieu_dit=="USN00-Réserve ornithologique (générique)" | DBt$Lieu_dit=="USN01-Artigues-Réserve ornithologique")) #Before 2007, not interesting for what I want

#DBt = subset(DBt,(((DBt$Protocol==1 | DBt$Protocol==2) & DBt$Lieu_dit=="USN00-Réserve ornithologique (générique)") | ((DBt$Protocol==1 | DBt$Protocol==2)  & DBt$Lieu_dit=="USN01-Artigues-Réserve ornithologique")))

species=sort(unique(DBt$Nom_latin))
for (s in species){
	print(s)
	for (a in 1973:2016){
		for (m in 1:12){
			plou=subset(DBt, Nom_latin==s & Mois==m & Annee==a)
			if (length(unique(plou$Jour))>1){
				print(paste(a,m,unique(plou$Jour)))	
			}
		}
	}
}
