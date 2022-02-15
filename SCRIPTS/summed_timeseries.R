##CP July 2018
##CP 25/07/19 There was a pb in wader def. It should also contain Calidris ferruginea
##CP 04/09/2019: Added Tringa, to compare them to Calidris and removed the "not wader" category that we don't use anymore
#CP 05/09/2019 Added "Philomachus pugnax" to calidris; changed DBt$Nom_latin to character to avoid any problem with factors, add a way to take into account biomasses (to avoid writing another script); removed sp that we will ignore later on. The index v2 in the output files corresponds to this version.

rm(list=ls())
graphics.off()

source("SCRIPTS/outputtimeseries_seasonal.R")
biomass=F #Do we weight the abundances with the biomasses? 
if(biomass){
end_bio="biomasses"
tab_bm=read.table("IN/Information_Trait_Oiseaux_20190903_allimportantbirds_meanmass.csv",sep=";",dec=".",header=T)
}else{
end_bio="abundances"
}
ignore_sp=FALSE #Do we remove the rare species from the beginning of the analysis?
if(ignore_sp){
end_ignore="wtoutrarespecies"
sp_to_ignore=c("Anas discors","Anas americana","Calidris melanotos","Calidris pusilla","Calidris ruficollis", "Calidris fuscicollis", "Calidris himantopus", "Burhinus oedicnemus","Phalaropus lobatus","Charadrius alexandrinus","Haematopus ostralegus","Calidris maritima","Aythya nyroca","Bucephala clangula","Melanitta nigra","Mergus serrator","Clangula hyemalis","Alopochen aegyptiaca", "Aix galericulata","Cygnus atratus","Tadorna ferruginea","Branta leucopsis","Anser fabalis","Anser albifrons","Cygnus cygnus","Mergus merganser","Anser brachyrhynchus")
}else{
end_ignore="wrarespecies"
sp_to_ignore=c()
}


DBt<-read.csv(file="IN/DBWithMonthlyPhotoTeich_completed.csv",header=TRUE,sep=",",dec=".")
DBt = subset(DBt,(DBt$Protocol==1 | DBt$Protocol==2) & (DBt$Lieu_dit=="USN00-Réserve ornithologique (générique)" | DBt$Lieu_dit=="USN01-Artigues-Réserve ornithologique") & (DBt$Annee>1980))
DBt<-DBt[!colnames(DBt) %in% c("X","Ref")] #edit du 14/04/2017 : there are 3 duplicates
DBt = unique.matrix(DBt) # delete with this line,  but I must check the monthly picture
DBt$Nom_latin=as.character(DBt$Nom_latin)
sp=unique(DBt$Nom_latin)

print('Anas')
sp_anas=sp[grep("Anas",sp)]
sp_anas=setdiff(sp_anas,sp_to_ignore)
tmp=subset(DBt,DBt$Nom_latin %in% sp_anas)
if(biomass){
tmp_bis=tmp
for(i in 1:nrow(tmp_bis)){
  tmp_bis$Nombre[i]=tmp_bis$Nombre[i]*tab_bm$Mean_mass[as.character(tab_bm$Species)==as.character(tmp_bis$Nom_latin[i])]
}
tmp=tmp_bis
}
anas_sum=sum_of_species(tmp,sp_anas,"Anas")
anas_season=season2_average(anas_sum)

print('Calidris')
sp_calidris=c(sp[grep("Calidris",sp)],"Philomachus pugnax")
sp_calidris=setdiff(sp_calidris,sp_to_ignore)
tmp=subset(DBt,DBt$Nom_latin %in% sp_calidris)
if(biomass){
tmp_bis=tmp
for(i in 1:nrow(tmp_bis)){
  tmp_bis$Nombre[i]=tmp_bis$Nombre[i]*tab_bm$Mean_mass[as.character(tab_bm$Species)==as.character(tmp_bis$Nom_latin[i])]
}
tmp=tmp_bis
}
calidris_sum=sum_of_species(tmp,sp_calidris,"Calidris")
calidris_season=season2_average(calidris_sum)

print('Tringa')
sp_tringa=sp[grep("Tringa",sp)]
sp_tringa=setdiff(sp_tringa,sp_to_ignore)
tmp=subset(DBt,DBt$Nom_latin %in% sp_tringa)
if(biomass){
tmp_bis=tmp
for(i in 1:nrow(tmp_bis)){
  tmp_bis$Nombre[i]=tmp_bis$Nombre[i]*tab_bm$Mean_mass[as.character(tab_bm$Species)==as.character(tmp_bis$Nom_latin[i])]
}
tmp=tmp_bis
}
tringa_sum=sum_of_species(tmp,sp_tringa,"Tringa")
tringa_season=season2_average(tringa_sum)

print('Waders')
limicoles = c("Recurvirostra avosetta","Limosa limosa","Limosa lapponica","Calidris temminckii","Calidris canutus",
              "Calidris alba","Calidris alpina","Calidris minuta","Calidris maritima" ,"Gallinago gallinago",
              "Tringa nebularia","Tringa erythropus","Tringa ochropus","Tringa totanus",
              "Tringa glareola","Actitis hypoleucos","Philomachus pugnax","Numenius arquata","Numenius phaeopus",
              "Himantopus himantopus","Charadrius hiaticula","Charadrius alexandrinus","Haematopus ostralegus",
              "Burhinus oedicnemus","Charadrius dubius","Phalaropus lobatus","Pluvialis squatarola",
              "Pluvialis apricaria","Arenaria interpres","Vanellus vanellus","Calidris ferruginea")
limicoles=setdiff(limicoles,sp_to_ignore)
tmp=subset(DBt,DBt$Nom_latin %in% limicoles)
if(biomass){
tmp_bis=tmp
for(i in 1:nrow(tmp_bis)){
  tmp_bis$Nombre[i]=tmp_bis$Nombre[i]*tab_bm$Mean_mass[as.character(tab_bm$Species)==as.character(tmp_bis$Nom_latin[i])]
}
tmp=tmp_bis
}
wader_sum=sum_of_species(DBt,limicoles,"Waders")
wader_season=season2_average(wader_sum)

print('Ducks')
tab_tmp=read.csv(file="./IN/Initial_files/data_ROT20160324.csv",header=TRUE,sep="\t")
sp_all=tab_tmp$Nom_latin
fam=tab_tmp$Famille
sp_duck=c(as.character(unique(sp_all[fam=="Anatidae"])),"Fulica atra")
sp_duck=setdiff(sp_duck,sp_to_ignore)
tmp=subset(DBt,DBt$Nom_latin %in% sp_duck)
if(biomass){
tmp_bis=tmp
for(i in 1:nrow(tmp_bis)){
  tmp_bis$Nombre[i]=tmp_bis$Nombre[i]*tab_bm$Mean_mass[as.character(tab_bm$Species)==as.character(tmp_bis$Nom_latin[i])]
}
tmp=tmp_bis
}
duck_sum=sum_of_species(tmp,sp_duck,"Ducks")
duck_season=season2_average(duck_sum)

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
oiseaux_Frequents_t=setdiff(oiseaux_Frequents_t,sp_to_ignore)
tmp=subset(DBt,DBt$Nom_latin %in% oiseaux_Frequents_t)
if(biomass){
tmp_bis=tmp
for(i in 1:nrow(tmp_bis)){
  tmp_bis$Nombre[i]=tmp_bis$Nombre[i]*tab_bm$Mean_mass[as.character(tab_bm$Species)==as.character(tmp_bis$Nom_latin[i])]
}
tmp=tmp_bis
}
freq_sum=sum_of_species(tmp,oiseaux_Frequents_t,'Freq')
freq_season=season2_average(freq_sum)

#Small trick for Cormorant (we don't do a sum but it gives the format I want)
print('Cormorant')
tmp=subset(DBt,DBt$Nom_latin %in% c("Phalacrocorax carbo"))
if(biomass){
tmp_bis=tmp
for(i in 1:nrow(tmp_bis)){
  tmp_bis$Nombre[i]=tmp_bis$Nombre[i]*tab_bm$Mean_mass[as.character(tab_bm$Species)==as.character(tmp_bis$Nom_latin[i])]
}
tmp=tmp_bis
}
cormorant_sum=sum_of_species(tmp,c("Phalacrocorax carbo"),'Cormorant')
cormorant_season=season2_average(cormorant_sum)

print('Finally Heron+Egret')
sp_her=c("Ardea cinerea","Egretta garzetta")
tmp=subset(DBt,DBt$Nom_latin %in% sp_her)
if(biomass){
tmp_bis=tmp
for(i in 1:nrow(tmp_bis)){
  tmp_bis$Nombre[i]=tmp_bis$Nombre[i]*tab_bm$Mean_mass[as.character(tab_bm$Species)==as.character(tmp_bis$Nom_latin[i])]
}
tmp=tmp_bis
}
heregr_sum=sum_of_species(tmp,c("Ardea cinerea","Egretta garzetta"),'HeronEgret')
heregr_season=season2_average(heregr_sum)

#all_sum=rbind(anas_sum,calidris_sum,wader_sum,duck_sum,other_than_wader_sum,freq_sum,cormorant_sum,heregr_sum)
all_sum=rbind(anas_sum,calidris_sum,wader_sum,duck_sum,tringa_sum,freq_sum,cormorant_sum,heregr_sum)
write.table(all_sum,paste("IN/summed_",end_bio,"_v2_",end_ignore,".csv",sep=""),sep=";")


table_seasonal_cold=matrix(NA,nrow=length(anas_season[,"Cold",]),ncol=8)
rownames(table_seasonal_cold)=dimnames(anas_season)[[3]]
#colnames(table_seasonal_cold)=c('Anas','Calidris','Waders','Ducks','Not-waders','Freq','Cormorant','HeronEgret')
colnames(table_seasonal_cold)=c('Anas','Calidris','Waders','Ducks','Tringa','Freq','Cormorant','HeronEgret')
table_seasonal_cold[,'Anas']=anas_season[,'Cold',]
table_seasonal_cold[,'Calidris']=calidris_season[,'Cold',]
table_seasonal_cold[,'Waders']=wader_season[,'Cold',]
table_seasonal_cold[,'Ducks']=duck_season[,'Cold',]
#table_seasonal_cold[,'Not-waders']=other_than_wader_season[,'Cold',]
table_seasonal_cold[,'Tringa']=tringa_season[,'Cold',]
table_seasonal_cold[,'Freq']=freq_season[,'Cold',]
table_seasonal_cold[,'Cormorant']=cormorant_season[,'Cold',]
table_seasonal_cold[,'HeronEgret']=heregr_season[,'Cold',]
write.table(table_seasonal_cold,paste("IN/coldseason_",end_bio,"_summed_v2_",end_ignore,".csv",sep=""),sep=";")

table_seasonal_warm=matrix(NA,nrow=length(anas_season[,"Warm",]),ncol=8)
rownames(table_seasonal_warm)=dimnames(anas_season)[[3]]
#colnames(table_seasonal_warm)=c('Anas','Calidris','Waders','Ducks','Not-waders','Freq','Cormorant','HeronEgret')
colnames(table_seasonal_warm)=c('Anas','Calidris','Waders','Ducks','Tringa','Freq','Cormorant','HeronEgret')
table_seasonal_warm[,'Anas']=anas_season[,'Warm',]
table_seasonal_warm[,'Calidris']=calidris_season[,'Warm',]
table_seasonal_warm[,'Waders']=wader_season[,'Warm',]
table_seasonal_warm[,'Ducks']=duck_season[,'Warm',]
#table_seasonal_warm[,'Not-waders']=other_than_wader_season[,'Warm',]
table_seasonal_warm[,'Tringa']=tringa_season[,'Warm',]
table_seasonal_warm[,'Freq']=freq_season[,'Warm',]
table_seasonal_warm[,'Cormorant']=cormorant_season[,'Warm',]
table_seasonal_warm[,'HeronEgret']=heregr_season[,'Warm',]
write.table(table_seasonal_warm,paste("IN/warmseason_",end_bio,"_summed_v2_",end_ignore,".csv",sep=""),sep=";")

#Produce data frames and not table
building_line='data_frame_cold=rbind('
for (nem in colnames(table_seasonal_cold)){
	building_line=paste(building_line,'cbind(rep(\"',nem,'\",',dim(table_seasonal_cold)[1],'),table_seasonal_cold[,\"',nem,'\"]),',sep="")
}
building_line=substr(building_line,1,nchar(building_line)-1)
building_line=paste(building_line,')',sep="")
eval(parse(text=building_line))
data_frame_cold=cbind(rownames(data_frame_cold),data_frame_cold)
write.table(data_frame_cold,paste("IN/coldseason_",end_bio,"_asdataframe_summed_v2_",end_ignore,".csv",sep=""),sep=";",col.names=c("Date","Species","Abundance"),row.names=F)

#Produce data frames and not table
building_line='data_frame_warm=rbind('
for (nem in colnames(table_seasonal_warm)){
        building_line=paste(building_line,'cbind(rep(\"',nem,'\",',dim(table_seasonal_warm)[1],'),table_seasonal_warm[,\"',nem,'\"]),',sep="")
}
building_line=substr(building_line,1,nchar(building_line)-1)
building_line=paste(building_line,')',sep="")
eval(parse(text=building_line))
data_frame_warm=cbind(rownames(data_frame_warm),data_frame_warm)
write.table(data_frame_warm,paste("IN/warmseason_",end_bio,"_asdataframe_summed_v2_",end_ignore,".csv",sep=""),sep=";",col.names=c("Date","Species","Abundance"),row.names=F)

