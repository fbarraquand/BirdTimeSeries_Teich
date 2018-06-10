#CP 2018 : To avoid problems with synchrony values. 

#We want the synchrony values inside the groups of Anas, Calidris, waders and frequent birds
#I'm relaunching this to (a) have 'readable' tables and values of synchronys ignoring the values before 1981

#UNFINISHED : need to put evtg in a dataframe

rm(list=ls())
graphics.off()

source("SCRIPTS/outputtimeseries_seasonal.R")


DBt<-read.csv(file="/home/cpicoche/Documents/Birds/BirdTimeSeries_Teich/IN/DBWithMonthlyPhotoTeich_completed.csv",header=TRUE,sep=",",dec=".")
DBt = subset(DBt,(DBt$Protocol==1 | DBt$Protocol==2) & (DBt$Lieu_dit=="USN00-Réserve ornithologique (générique)" | DBt$Lieu_dit=="USN01-Artigues-Réserve ornithologique") & (DBt$Annee>1980))
DBt<-DBt[!colnames(DBt) %in% c("X","Ref")] #edit du 14/04/2017 : there are 3 duplicates
DBt = unique.matrix(DBt) # delete with this line,  but I must check the monthly picture
sp=unique(DBt$Nom_latin)

print('Anas')
sp_anas=sp[grep("Anas",sp)]
tmp=subset(DBt,DBt$Nom_latin %in% sp_anas)
tab_anas=season2_average(tmp)


print('Calidris')
sp_calidris=sp[grep("Calidris",sp)]
tmp=subset(DBt,DBt$Nom_latin %in% sp_calidris)
tab_calidris=season2_average(tmp)

print('Waders')
limicoles = c("Recurvirostra avosetta","Limosa limosa","Limosa lapponica","Calidris temminckii","Calidris canutus",
              "Calidris alba","Calidris alpina","Calidris minuta","Calidris maritima" ,"Gallinago gallinago",
              "Tringa flavipes","Tringa nebularia","Tringa erythropus","Tringa ochropus","Tringa totanus",
              "Tringa glareola","Actitis hypoleucos","Philomachus pugnax","Numenius arquata","Numenius phaeopus",
              "Himantopus himantopus","Charadrius hiaticula","Charadrius alexandrinus","Haematopus ostralegus",
              "Burhinus oedicnemus","Charadrius dubius","Phalaropus lobatus","Pluvialis squatarola",
              "Pluvialis apricaria","Arenaria interpres","Vanellus vanellus")
tmp=subset(DBt,DBt$Nom_latin %in% limicoles)
tab_waders=season2_average(tmp)

print('Frequent birds')
vec_n_obs_oiseaux_t=rep(0,length(sp_all))
for (i in 1:length(sp_all)){
  vec_n_obs_oiseaux_t[i]=sum(as.character(DBt$Nom_latin)==sp_all[i],na.rm = TRUE)
}
oiseaux_Frequents_t=sp_all[vec_n_obs_oiseaux_t>75]
tmp=subset(DBt,DBt$Nom_latin %in% oiseaux_Frequents)
tab_waders=season2_average(tmp)

