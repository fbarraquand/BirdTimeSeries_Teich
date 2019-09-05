rm(list=ls())
graphics.off()

DBt<-read.csv(file="/home/cpicoche/Documents/Birds/BirdTimeSeries_Teich/IN/DBWithMonthlyPhotoTeich_completed.csv",header=TRUE,sep=",",dec=".")
DBt = subset(DBt,(DBt$Protocol==1 | DBt$Protocol==2) & (DBt$Lieu_dit=="USN00-Réserve ornithologique (générique)" | DBt$Lieu_dit=="USN01-Artigues-Réserve ornithologique")  )

#No Bartramia, Coenocorypha
a=DBt$Nombre[grep("Limnodromus scolopaceus",DBt$Nom_latin)]
a=a[a>0]
print(length(a))
print(mean(a))
a=DBt$Nombre[grep("Lymnocryptes minimus",DBt$Nom_latin)]
a=a[a>0]
print(length(a))
print(mean(a))
a=DBt$Nombre[grep("Scolopax",DBt$Nom_latin)]
a=a[a>0]
print(length(a))
print(mean(a))
a=DBt$Nombre[grep("Xenus",DBt$Nom_latin)]
a=a[a>0]
print(length(a))
print(mean(a))
a=DBt$Nombre[grep("Actitis hypoleucos",DBt$Nom_latin)]
a=a[a>0]
print(length(a))
print(mean(a))
a=DBt$Nombre[grep("Arenaria interpres",DBt$Nom_latin)]
a=a[a>0]
print(length(a))
print(mean(a))

tab_tmp=read.csv(file="./IN/Initial_files/data_ROT20160324.csv",header=TRUE,sep="\t")
sp=tab_tmp$Nom_latin
fam=tab_tmp$Famille

limicoles = sort(c("Recurvirostra avosetta","Limosa limosa","Limosa lapponica","Calidris temminckii","Calidris canutus",
              "Calidris alba","Calidris alpina","Calidris minuta","Calidris maritima" ,"Gallinago gallinago",
              "Tringa flavipes","Tringa nebularia","Tringa erythropus","Tringa ochropus","Tringa totanus",
              "Tringa glareola","Actitis hypoleucos","Philomachus pugnax","Numenius arquata","Numenius phaeopus",
              "Himantopus himantopus","Charadrius hiaticula","Charadrius alexandrinus","Haematopus ostralegus",
              "Charadrius dubius","Phalaropus lobatus","Pluvialis squatarola",
              "Pluvialis apricaria","Arenaria interpres","Vanellus vanellus"))

#Removed Burrhinus oedicnemus
fam_limicoles=rep(NA,length(limicoles))
for (f in 1:length(limicoles)){
	fam_limicoles[f]=as.character(unique(fam[sp==limicoles[f]]))
	if(fam_limicoles[f]!="Scolopacidae"){
		print(limicoles[f])
		print(fam_limicoles[f])
	}
}
#
