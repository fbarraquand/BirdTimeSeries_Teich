#CP 2018 : To avoid problems with synchrony values. 

#We want the synchrony values inside the groups of Anas, Calidris, waders and frequent birds
#I'm relaunching this to (a) have 'readable' tables and values of synchronys ignoring the values before 1981 and b) have all birds counted at the same time

#CP 25/07/2019 Added "Calidris ferruginea" to waders
#CP 25/07/2019 Used the normalisation_timeseries.r script, just added the biomasses

rm(list=ls())
graphics.off()

source("SCRIPTS/outputtimeseries_seasonal.R")
sp_to_ignore=c("Anas discors","Anas americana","Calidris melanotos","Calidris pusilla","Calidris ruficollis", "Calidris fuscicollis", "Calidris himantopus", "Burhinus oedicnemus","Phalaropus lobatus","Charadrius alexandrinus","Haematopus ostralegus","Calidris maritima","Aythya nyroca","Bucephala clangula","Melanitta nigra","Mergus serrator","Clangula hyemalis","Alopochen aegyptiaca", "Aix galericulata","Cygnus atratus","Tadorna ferruginea","Branta leucopsis","Anser fabalis","Anser albifrons","Cygnus cygnus","Mergus merganser","Anser brachyrhynchus")

tab_bm=read.table("IN/Information_Trait_Oiseaux_20190903_allimportantbirds_meanmass.csv",sep=";",dec=".",header=T)


DBt<-read.csv(file="/home/cpicoche/Documents/Birds/BirdTimeSeries_Teich/IN/DBWithMonthlyPhotoTeich_completed.csv",header=TRUE,sep=",",dec=".")
DBt = subset(DBt,(DBt$Protocol==1 | DBt$Protocol==2) & (DBt$Lieu_dit=="USN00-Réserve ornithologique (générique)" | DBt$Lieu_dit=="USN01-Artigues-Réserve ornithologique") & (DBt$Annee>1980))
DBt<-DBt[!colnames(DBt) %in% c("X","Ref")] #edit du 14/04/2017 : there are 3 duplicates
DBt = unique.matrix(DBt) # delete with this line,  but I must check the monthly picture
sp=unique(DBt$Nom_latin)

print('Anas')
sp_anas=sp[grep("Anas",sp)]
sp_anas=setdiff(sp_anas,sp_to_ignore)
tab=subset(DBt,(DBt$Nom_latin %in% sp_anas))
tab_bis=tab
for(i in 1:nrow(tab_bis)){
  tab_bis$Nombre[i]=tab_bis$Nombre[i]*tab_bm$Mean_mass[as.character(tab_bm$Species)==as.character(tab_bis$Nom_latin[i])]
}
tab_anas=season2_average(tab_bis)

print('Calidris')
sp_calidris=sp[grep("Calidris",sp)]
sp_calidris=setdiff(sp_calidris,sp_to_ignore)
tab=subset(DBt,(DBt$Nom_latin %in% sp_calidris))
tab_bis=tab
for(i in 1:nrow(tab_bis)){
  tab_bis$Nombre[i]=tab_bis$Nombre[i]*tab_bm$Mean_mass[as.character(tab_bm$Species)==as.character(tab_bis$Nom_latin[i])]
}
tab_calidris=season2_average(tab_bis)

print('Waders')
limicoles = c("Recurvirostra avosetta","Limosa limosa","Limosa lapponica","Calidris temminckii","Calidris canutus",
              "Calidris alba","Calidris alpina","Calidris minuta","Calidris maritima" ,"Gallinago gallinago",
              "Tringa nebularia","Tringa erythropus","Tringa ochropus","Tringa totanus",
              "Tringa glareola","Actitis hypoleucos","Philomachus pugnax","Numenius arquata","Numenius phaeopus",
              "Himantopus himantopus","Charadrius hiaticula","Charadrius alexandrinus","Haematopus ostralegus",
              "Burhinus oedicnemus","Charadrius dubius","Phalaropus lobatus","Pluvialis squatarola",
              "Pluvialis apricaria","Arenaria interpres","Vanellus vanellus","Calidris ferruginea")
#Tringa flavipes never appers with Protocol>0
limicoles=setdiff(limicoles,sp_to_ignore)
tab=subset(DBt,(DBt$Nom_latin %in% limicoles))
tab_bis=tab
for(i in 1:nrow(tab_bis)){
  tab_bis$Nombre[i]=tab_bis$Nombre[i]*tab_bm$Mean_mass[as.character(tab_bm$Species)==as.character(tab_bis$Nom_latin[i])]
}
tab_waders=season2_average(tab_bis)

print('Frequent birds')
vec_n_obs_oiseaux_t=rep(0,length(sp))
for (i in 1:length(sp)){
  vec_n_obs_oiseaux_t[i]=sum(as.character(DBt$Nom_latin)==sp[i],na.rm = TRUE)
}
oiseaux_Frequents_t=sp[vec_n_obs_oiseaux_t>75]
oiseaux_Frequents_t=setdiff(oiseaux_Frequents_t,sp_to_ignore)
tab=subset(DBt,(DBt$Nom_latin %in% oiseaux_Frequents_t))
tab_bis=tab
for(i in 1:nrow(tab_bis)){
  tab_bis$Nombre[i]=tab_bis$Nombre[i]*tab_bm$Mean_mass[as.character(tab_bm$Species)==as.character(tab_bis$Nom_latin[i])]
}
tab_freq=season2_average(tab_bis)

dd=as.character(1981:2015)

dates=rep(dd,length(sp_anas))
sp_all=c()
abundance_cold=c()
abundance_warm=c()
for (i in 1:length(sp_anas)){
  ss=as.character(sp_anas[i])
  sp_all=c(sp_all,rep(ss,length(dd)))
  abundance_cold=c(abundance_cold,tab_anas[ss,'Cold',dd])
  abundance_warm=c(abundance_warm,tab_anas[ss,'Warm',dd])
}
data_frame_anas_cold=data.frame(dates,sp_all,abundance_cold)
data_frame_anas_warm=data.frame(dates,sp_all,abundance_warm)
write.csv(data_frame_anas_cold,"IN/coldseason_anas_detailed_biomasses.txt")
write.csv(data_frame_anas_warm,"IN/warmseason_anas_detailed_biomasses.txt")

dates=rep(dd,length(sp_calidris))
sp_all=c()
abundance_cold=c()
abundance_warm=c()
for (i in 1:length(sp_calidris)){
  ss=as.character(sp_calidris[i])
  sp_all=c(sp_all,rep(ss,length(dd)))
  abundance_cold=c(abundance_cold,tab_calidris[ss,'Cold',dd])
  abundance_warm=c(abundance_warm,tab_calidris[ss,'Warm',dd])
}
data_frame_calidris_cold=data.frame(dates,sp_all,abundance_cold)
data_frame_calidris_warm=data.frame(dates,sp_all,abundance_warm)
write.csv(data_frame_calidris_cold,"IN/coldseason_calidris_detailed_biomasses.txt")
write.csv(data_frame_calidris_warm,"IN/warmseason_calidris_detailed_biomasses.txt")

dates=rep(dd,length(limicoles))
sp_all=c()
abundance_cold=c()
abundance_warm=c()
for (i in 1:length(limicoles)){
  ss=as.character(limicoles[i])
  sp_all=c(sp_all,rep(ss,length(dd)))
  abundance_cold=c(abundance_cold,tab_waders[ss,'Cold',dd])
  abundance_warm=c(abundance_warm,tab_waders[ss,'Warm',dd])
}
data_frame_waders_cold=data.frame(dates,sp_all,abundance_cold)
data_frame_waders_warm=data.frame(dates,sp_all,abundance_warm)
write.csv(data_frame_waders_cold,"IN/coldseason_waders_detailed_biomasses.txt")
write.csv(data_frame_waders_warm,"IN/warmseason_waders_detailed_biomasses.txt")

dates=rep(dd,length(oiseaux_Frequents_t))
sp_all=c()
abundance_cold=c()
abundance_warm=c()
for (i in 1:length(oiseaux_Frequents_t)){
  ss=as.character(oiseaux_Frequents_t[i])
  sp_all=c(sp_all,rep(ss,length(dd)))
  abundance_cold=c(abundance_cold,tab_freq[ss,'Cold',dd])
  abundance_warm=c(abundance_warm,tab_freq[ss,'Warm',dd])
}
data_frame_freq_cold=data.frame(dates,sp_all,abundance_cold)
data_frame_freq_warm=data.frame(dates,sp_all,abundance_warm)
write.csv(data_frame_freq_cold,"IN/coldseason_freq_detailed_biomasses.txt")
write.csv(data_frame_freq_warm,"IN/warmseason_freq_detailed_biomasses.txt")

print('Ducks')
tab_tmp=read.csv(file="./IN/Initial_files/data_ROT20160324.csv",header=TRUE,sep="\t")
sp_all=tab_tmp$Nom_latin
fam=tab_tmp$Famille
sp_duck=c(as.character(unique(sp_all[fam=="Anatidae"])),"Fulica atra")
sp_duck=intersect(sp_duck,DBt$Nom_latin)
sp_duck=setdiff(sp_duck,sp_to_ignore)
tab=subset(DBt,(DBt$Nom_latin %in% sp_duck))
tab_bis=tab
for(i in 1:nrow(tab_bis)){
  tab_bis$Nombre[i]=tab_bis$Nombre[i]*tab_bm$Mean_mass[as.character(tab_bm$Species)==as.character(tab_bis$Nom_latin[i])]
}
tab_duck=season2_average(tab_bis)

dates=rep(dd,length(sp_duck))
sp_all=c()
abundance_cold=c()
abundance_warm=c()
for (i in 1:length(sp_duck)){
  ss=as.character(sp_duck[i])
  sp_all=c(sp_all,rep(ss,length(dd)))
  abundance_cold=c(abundance_cold,tab_duck[ss,'Cold',dd])
  abundance_warm=c(abundance_warm,tab_duck[ss,'Warm',dd])
}
data_frame_duck_cold=data.frame(dates,sp_all,abundance_cold)
data_frame_duck_warm=data.frame(dates,sp_all,abundance_warm)
write.csv(data_frame_duck_cold,"IN/coldseason_duck_detailed_biomasses.txt")
write.csv(data_frame_duck_warm,"IN/warmseason_duck_detailed_biomasses.txt")



