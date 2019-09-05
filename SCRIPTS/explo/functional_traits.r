rm(list=ls())
graphics.off()

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
SpeciesL = as.character(unique(DBt$Nom_latin)) #-> 280
SpeciesF = as.character(unique(DBt$Nom_espece)) #même chose mais avec les noms français

#minimun d'observation pour garder une espèce
vec_n_obs_oiseaux_t=rep(0,length(SpeciesL))
for (i in 1:length(SpeciesL)){
  vec_n_obs_oiseaux_t[i]=sum(as.character(DBt$Nom_latin)==SpeciesL[i],na.rm = TRUE)
}
oiseaux_Frequents_t=SpeciesL[vec_n_obs_oiseaux_t>75] #60/279


tmp<-read.csv(file="/home/cpicoche/Documents/Birds/BirdTimeSeries_Teich/IN/BirdFuncDat.csv",header=TRUE,sep=";",dec=".")
#tmp<-read.table(file="/home/cpicoche/Documents/Birds/BirdTimeSeries_Teich/IN/BirdFuncDat.txt",header=TRUE,sep="\t",dec=".") #Hamish Wilman, Jonathan Belmaker, Jennifer Simpson, Carolina de la Rosa, Marcelo M. Rivadeneira, and Walter Jetz. 2014. EltonTraits 1.0: Species-level foraging attributes of the world's birds and mammals. Ecology 95:2027. http://dx.doi.org/10.1890/13-1917.1

sok=0
stot=0
for (s in SpeciesL){
	stot=stot+1
	if(s %in% tmp$Scientific){
		sok=sok+1
	}
}
print(paste("For all birds",sok/stot))

sok=0
stot=0
for (s in oiseaux_Frequents_t){
        stot=stot+1
        if(s %in% tmp$Scientific){
                sok=sok+1
        }
}
print(paste("For frequent birds",sok/stot))

limicoles = sort(c("Recurvirostra avosetta","Limosa limosa","Limosa lapponica","Calidris temminckii","Calidris canutus",
              "Calidris alba","Calidris alpina","Calidris minuta","Calidris maritima" ,"Gallinago gallinago",
              "Tringa flavipes","Tringa nebularia","Tringa erythropus","Tringa ochropus","Tringa totanus",
              "Tringa glareola","Actitis hypoleucos","Philomachus pugnax","Numenius arquata","Numenius phaeopus",
              "Himantopus himantopus","Charadrius hiaticula","Charadrius alexandrinus","Haematopus ostralegus",
              "Charadrius dubius","Phalaropus lobatus","Pluvialis squatarola",
              "Pluvialis apricaria","Arenaria interpres","Vanellus vanellus"))

sok=0
stot=0
for (s in limicoles){
        stot=stot+1
        if(s %in% tmp$Scientific){
                sok=sok+1
        }
}
print(paste("For waders",sok/stot))

