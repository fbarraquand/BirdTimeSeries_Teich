##CP July 2018
##CP 25/07/19 There was a pb in wader def. It should also contain Calidris ferruginea
##CP 03/09/19 Now using the script "summed_timeseries" but weighting the numbers with the mass of each species

rm(list=ls())
graphics.off()

source("SCRIPTS/outputtimeseries_seasonal.R")

sp_to_ignore=c("Anas discors","Anas americana","Calidris melanotos","Calidris pusilla","Calidris ruficollis", "Calidris fuscicollis", "Calidris himantopus", "Burhinus oedicnemus","Phalaropus lobatus","Charadrius alexandrinus","Haematopus ostralegus","Calidris maritima","Aythya nyroca","Bucephala clangula","Melanitta nigra","Mergus serrator","Clangula hyemalis","Alopochen aegyptiaca", "Aix galericulata","Cygnus atratus","Tadorna ferruginea","Branta leucopsis","Anser fabalis","Anser albifrons","Cygnus cygnus","Mergus merganser","Anser brachyrhynchus")

DBt<-read.csv(file="/home/cpicoche/Documents/Birds/BirdTimeSeries_Teich/IN/DBWithMonthlyPhotoTeich_completed.csv",header=TRUE,sep=",",dec=".")
DBt = subset(DBt,(DBt$Protocol==1 | DBt$Protocol==2) & (DBt$Lieu_dit=="USN00-Réserve ornithologique (générique)" | DBt$Lieu_dit=="USN01-Artigues-Réserve ornithologique") & (DBt$Annee>1980))
DBt<-DBt[!colnames(DBt) %in% c("X","Ref")] #edit du 14/04/2017 : there are 3 duplicates
DBt = unique.matrix(DBt) # delete with this line,  but I must check the monthly picture
sp=unique(DBt$Nom_latin)

tab_bm=read.table("IN/Information_Trait_Oiseaux_20190903_allimportantbirds_meanmass.csv",sep=";",dec=".",header=T)

print('Anas')
sp_anas=sp[grep("Anas",sp)]
tab=subset(DBt,((DBt$Nom_latin %in% sp_anas)&!(DBt$Nom_latin %in% sp_to_ignore)))
tab_bis=tab
for(i in 1:nrow(tab_bis)){
  tab_bis$Nombre[i]=tab_bis$Nombre[i]*tab_bm$Mean_mass[as.character(tab_bm$Species)==as.character(tab_bis$Nom_latin[i])]
}
anas_sum=sum_of_species(tab_bis,sp_anas,"Anas")
anas_season=season2_average(anas_sum)

print('Calidris')
sp_calidris=sp[grep("Calidris",sp)]
tab=subset(DBt,((DBt$Nom_latin %in% sp_calidris)&!(DBt$Nom_latin %in% sp_to_ignore)))
tab_bis=tab
for(i in 1:nrow(tab_bis)){
  tab_bis$Nombre[i]=tab_bis$Nombre[i]*tab_bm$Mean_mass[as.character(tab_bm$Species)==as.character(tab_bis$Nom_latin[i])]
}
calidris_sum=sum_of_species(tab_bis,sp_calidris,"Calidris")
calidris_season=season2_average(calidris_sum)

print('Waders')
limicoles = c("Recurvirostra avosetta","Limosa limosa","Limosa lapponica","Calidiris temminckii","Calidris canutus",
              "Calidris alba","Calidris alpina","Calidris minuta","Calidris maritima" ,"Gallinago gallinago",
              "Tringa nebularia","Tringa erythropus","Tringa ochropus","Tringa totanus",
              "Tringa glareola","Actitis hypoleucos","Philomachus pugnax","Numenius arquata","Numenius phaeopus",
              "Himantopus himantopus","Charadrius hiaticula","Charadrius alexandrinus","Haematopus ostralegus",
              "Burhinus oedicnemus","Charadrius dubius","Phalaropus lobatus","Pluvialis squatarola",
              "Pluvialis apricaria","Arenaria interpres","Vanellus vanellus","Calidris ferruginea")
tab=subset(DBt,((DBt$Nom_latin %in% limicoles)&!(DBt$Nom_latin %in% sp_to_ignore)))
tab_bis=tab
for(i in 1:nrow(tab_bis)){
  tab_bis$Nombre[i]=tab_bis$Nombre[i]*tab_bm$Mean_mass[as.character(tab_bm$Species)==as.character(tab_bis$Nom_latin[i])]
}
wader_sum=sum_of_species(tab_bis,limicoles,"Waders")
wader_season=season2_average(wader_sum)

print('Ducks')
tab_tmp=read.csv(file="./IN/Initial_files/data_ROT20160324.csv",header=TRUE,sep="\t")
sp_all=tab_tmp$Nom_latin
fam=tab_tmp$Famille

sp_duck=c(as.character(unique(sp_all[fam=="Anatidae"])),"Fulica atra")
tab=subset(DBt,((DBt$Nom_latin %in% sp_duck)&!(DBt$Nom_latin %in% sp_to_ignore)))
tab_bis=tab
for(i in 1:nrow(tab_bis)){
  tab_bis$Nombre[i]=tab_bis$Nombre[i]*tab_bm$Mean_mass[as.character(tab_bm$Species)==as.character(tab_bis$Nom_latin[i])]
}
duck_sum=sum_of_species(tab_bis,sp_duck,"Ducks")
duck_season=season2_average(duck_sum)

### We don't want to look at things that are not waders, it's not worh ti
if(1==0){
print('Other birds than waders') ##This will be the longest one
sp_other=setdiff(sp_all,limicoles)
other_than_wader_sum=sum_of_species(DBt,sp_other,"Not-waders")
other_than_wader_season=season2_average(other_than_wader_sum)
}

print('Frequent birds')
vec_n_obs_oiseaux_t=rep(0,length(sp_all))
for (i in 1:length(sp_all)){
  vec_n_obs_oiseaux_t[i]=sum(as.character(DBt$Nom_latin)==sp_all[i],na.rm = TRUE)
}
oiseaux_Frequents_t=sp_all[vec_n_obs_oiseaux_t>75]
oiseaux_Frequents_t=oiseaux_Frequents_t[!(oiseaux_Frequents_t %in% limicoles)]
tab=subset(DBt,((DBt$Nom_latin %in% oiseaux_Frequents_t)&!(DBt$Nom_latin %in% sp_to_ignore)))
tab_bis=tab
for(i in 1:nrow(tab_bis)){
  tab_bis$Nombre[i]=tab_bis$Nombre[i]*tab_bm$Mean_mass[as.character(tab_bm$Species)==as.character(tab_bis$Nom_latin[i])]
}
freq_sum=sum_of_species(tab_bis,oiseaux_Frequents_t,'Freq')
freq_season=season2_average(freq_sum)

#Small trick for Cormorant (we don't do a sum but it gives the format I want)
print('Cormorant')
tab=subset(DBt,DBt$Nom_latin=="Phalacrocorax carbo")
tab_bis=tab
for(i in 1:nrow(tab_bis)){
  tab_bis$Nombre[i]=tab_bis$Nombre[i]*tab_bm$Mean_mass[as.character(tab_bm$Species)==as.character(tab_bis$Nom_latin[i])]
}
cormorant_sum=sum_of_species(tab_bis,c("Phalacrocorax carbo"),'Cormorant')
cormorant_season=season2_average(cormorant_sum)

print('Finally Heron+Egret')
tab=subset(DBt,DBt$Nom_latin %in% c("Ardea cinerea","Egretta garzetta"))
tab_bis=tab
for(i in 1:nrow(tab_bis)){
  tab_bis$Nombre[i]=tab_bis$Nombre[i]*tab_bm$Mean_mass[as.character(tab_bm$Species)==as.character(tab_bis$Nom_latin[i])]
}
heregr_sum=sum_of_species(DBt,c("Ardea cinerea","Egretta garzetta"),'HeronEgret')
heregr_season=season2_average(heregr_sum)

#### Hereafter, I remove all mentions of the "Not wader" category that we won't use anymore, and that would take too much time to retrieve all info on biomasses
all_sum=rbind(anas_sum,calidris_sum,wader_sum,duck_sum,freq_sum,cormorant_sum,heregr_sum)
write.table(all_sum,'IN/summed_abundances_biomasses.csv',sep=";")

table_seasonal_cold=matrix(NA,nrow=length(anas_season[,"Cold",]),ncol=7)
rownames(table_seasonal_cold)=dimnames(anas_season)[[3]]
colnames(table_seasonal_cold)=c('Anas','Calidris','Waders','Ducks','Freq','Cormorant','HeronEgret')
table_seasonal_cold[,'Anas']=anas_season[,'Cold',]
table_seasonal_cold[,'Calidris']=calidris_season[,'Cold',]
table_seasonal_cold[,'Waders']=wader_season[,'Cold',]
table_seasonal_cold[,'Ducks']=duck_season[,'Cold',]
table_seasonal_cold[,'Freq']=freq_season[,'Cold',]
table_seasonal_cold[,'Cormorant']=cormorant_season[,'Cold',]
table_seasonal_cold[,'HeronEgret']=heregr_season[,'Cold',]
write.table(table_seasonal_cold,"IN/coldseason_abundances_summed_biomasses.csv",sep=";")

table_seasonal_warm=matrix(NA,nrow=length(anas_season[,"Warm",]),ncol=7)
rownames(table_seasonal_warm)=dimnames(anas_season)[[3]]
colnames(table_seasonal_warm)=c('Anas','Calidris','Waders','Ducks','Freq','Cormorant','HeronEgret')
table_seasonal_warm[,'Anas']=anas_season[,'Warm',]
table_seasonal_warm[,'Calidris']=calidris_season[,'Warm',]
table_seasonal_warm[,'Waders']=wader_season[,'Warm',]
table_seasonal_warm[,'Ducks']=duck_season[,'Warm',]
table_seasonal_warm[,'Freq']=freq_season[,'Warm',]
table_seasonal_warm[,'Cormorant']=cormorant_season[,'Warm',]
table_seasonal_warm[,'HeronEgret']=heregr_season[,'Warm',]
write.table(table_seasonal_warm,"IN/warmseason_abundances_summed_biomasses.csv",sep=";")

#Produce data frames and not table
building_line='data_frame_cold=rbind('
for (nem in colnames(table_seasonal_cold)){
  building_line=paste(building_line,'cbind(rep(\"',nem,'\",',dim(table_seasonal_cold)[1],'),table_seasonal_cold[,\"',nem,'\"]),',sep="")
}
building_line=substr(building_line,1,nchar(building_line)-1)
building_line=paste(building_line,')',sep="")
eval(parse(text=building_line))
data_frame_cold=cbind(rownames(data_frame_cold),data_frame_cold)
write.table(data_frame_cold,"IN/coldseason_abundances_asdataframe_summed_biomasses.csv",sep=";",col.names=c("Date","Species","Abundance"),row.names=F)

#Produce data frames and not table
building_line='data_frame_warm=rbind('
for (nem in colnames(table_seasonal_warm)){
  building_line=paste(building_line,'cbind(rep(\"',nem,'\",',dim(table_seasonal_warm)[1],'),table_seasonal_warm[,\"',nem,'\"]),',sep="")
}
building_line=substr(building_line,1,nchar(building_line)-1)
building_line=paste(building_line,')',sep="")
eval(parse(text=building_line))
data_frame_warm=cbind(rownames(data_frame_warm),data_frame_warm)
write.table(data_frame_warm,"IN/warmseason_abundances_asdataframe_summed_biomasses.csv",sep=";",col.names=c("Date","Species","Abundance"),row.names=F)

