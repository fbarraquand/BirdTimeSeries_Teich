#CP 2018 : To avoid problems with synchrony values. 

#We want the synchrony values inside the groups of Anas, Calidris, waders and frequent birds
#I'm relaunching this to (a) have 'readable' tables and values of synchronys ignoring the values before 1981 and b) have all birds counted at the same time

#CP 25/07/2019 Added "Calidris ferruginea" to waders
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
ignore_sp=F #Do we remove the rare species from the beginning of the analysis?
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
tab_anas=season2_average(tmp)

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
tab_calidris=season2_average(tmp)

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
tmp=subset(DBt,DBt$Nom_latin %in% limicoles)
if(biomass){
tmp_bis=tmp
for(i in 1:nrow(tmp_bis)){
  tmp_bis$Nombre[i]=tmp_bis$Nombre[i]*tab_bm$Mean_mass[as.character(tab_bm$Species)==as.character(tmp_bis$Nom_latin[i])]
}
tmp=tmp_bis
}
tab_waders=season2_average(tmp)

print('Frequent birds')
vec_n_obs_oiseaux_t=rep(0,length(sp))
for (i in 1:length(sp)){
  vec_n_obs_oiseaux_t[i]=sum(as.character(DBt$Nom_latin)==sp[i],na.rm = TRUE)
}
oiseaux_Frequents_t=sp[vec_n_obs_oiseaux_t>75]
oiseaux_Frequents_t=setdiff(oiseaux_Frequents_t,sp_to_ignore)
tmp=subset(DBt,DBt$Nom_latin %in% oiseaux_Frequents_t)
if(biomass){
tmp_bis=tmp
for(i in 1:nrow(tmp_bis)){
  tmp_bis$Nombre[i]=tmp_bis$Nombre[i]*tab_bm$Mean_mass[as.character(tab_bm$Species)==as.character(tmp_bis$Nom_latin[i])]
}
tmp=tmp_bis
}
tab_freq=season2_average(tmp)

dd=as.character(1981:2015)
	
dates=rep(dd,length(sp_anas))
sp_all=c()
abundance_cold=c()
abundance_warm=c()
pdf(paste("OUT/closer_look_anas_v2_",end_bio,"_",end_ignore,".pdf",sep=""))
par(mfrow=c(3,3))
for (i in 1:length(sp_anas)){
	ss=as.character(sp_anas[i])
	sp_all=c(sp_all,rep(ss,length(dd)))
	abundance_cold=c(abundance_cold,tab_anas[ss,'Cold',dd])
	abundance_warm=c(abundance_warm,tab_anas[ss,'Warm',dd])
	plot(1981:2015,tab_anas[ss,'Cold',dd],ylim=c(0,max(tab_anas[ss,,],na.rm=T)),xlim=c(1981,2015),col="blue",t="o",pch=16,ylab=ss,xlab="")
	lines(1981:2015,tab_anas[ss,'Warm',dd],col="red",t="o",pch=16)
}
data_frame_anas_cold=data.frame(dates,sp_all,abundance_cold)
data_frame_anas_warm=data.frame(dates,sp_all,abundance_warm)
write.csv(data_frame_anas_cold,paste("IN/coldseason_anas_detailed_v2_",end_bio,"_",end_ignore,".txt",sep=""))
write.csv(data_frame_anas_warm,paste("IN/warmseason_anas_detailed_v2_",end_bio,"_",end_ignore,".txt",sep=""))
dev.off()

dates=rep(dd,length(sp_calidris))
sp_all=c()
abundance_cold=c()
abundance_warm=c()
#pdf("OUT/closer_look_calidris.pdf")
pdf(paste("OUT/closer_look_calidris_v2_",end_bio,"_",end_ignore,".pdf",sep=""))
par(mfrow=c(4,3))
for (i in 1:length(sp_calidris)){
        ss=as.character(sp_calidris[i])
        sp_all=c(sp_all,rep(ss,length(dd)))
        abundance_cold=c(abundance_cold,tab_calidris[ss,'Cold',dd])
        abundance_warm=c(abundance_warm,tab_calidris[ss,'Warm',dd])
        plot(1981:2015,tab_calidris[ss,'Cold',dd],ylim=c(0,max(tab_calidris[ss,,],na.rm=T)),xlim=c(1981,2015),col="blue",t="o",pch=16,ylab=ss,xlab="")
        lines(1981:2015,tab_calidris[ss,'Warm',dd],col="red",t="o",pch=16)
}
data_frame_calidris_cold=data.frame(dates,sp_all,abundance_cold)
data_frame_calidris_warm=data.frame(dates,sp_all,abundance_warm)
write.csv(data_frame_calidris_cold,paste("IN/coldseason_calidris_detailed_v2_",end_bio,"_",end_ignore,".txt",sep=""))
write.csv(data_frame_calidris_warm,paste("IN/warmseason_calidris_detailed_v2_",end_bio,"_",end_ignore,".txt",sep=""))
dev.off()

dates=rep(dd,length(limicoles))
sp_all=c()
abundance_cold=c()
abundance_warm=c()
#pdf("OUT/closer_look_waders.pdf")
pdf(paste("OUT/closer_look_waders_v2_",end_bio,"_",end_ignore,".pdf",sep=""))
par(mfrow=c(4,4))
for (i in 1:length(limicoles)){
        ss=as.character(limicoles[i])
        sp_all=c(sp_all,rep(ss,length(dd)))
        abundance_cold=c(abundance_cold,tab_waders[ss,'Cold',dd])
        abundance_warm=c(abundance_warm,tab_waders[ss,'Warm',dd])
	if(i==17){
	par(mfrow=c(4,4))
	}
        plot(1981:2015,tab_waders[ss,'Cold',dd],ylim=c(0,max(tab_waders[ss,,],na.rm=T)),xlim=c(1981,2015),col="blue",t="o",pch=16,ylab=ss,xlab="")
        lines(1981:2015,tab_waders[ss,'Warm',dd],col="red",t="o",pch=16)
}
data_frame_waders_cold=data.frame(dates,sp_all,abundance_cold)
data_frame_waders_warm=data.frame(dates,sp_all,abundance_warm)
write.csv(data_frame_waders_cold,paste("IN/coldseason_waders_detailed_v2_",end_bio,"_",end_ignore,".txt",sep=""))
write.csv(data_frame_waders_warm,paste("IN/warmseason_waders_detailed_v2_",end_bio,"_",end_ignore,".txt",sep=""))
dev.off()

dates=rep(dd,length(oiseaux_Frequents_t))
sp_all=c()
abundance_cold=c()
abundance_warm=c()
pdf(paste("OUT/closer_look_freq_v2_",end_bio,"_",end_ignore,".pdf",sep=""))
for (i in 1:length(oiseaux_Frequents_t)){
	if(i%%12==1){
	par(mfrow=c(4,3))
	}
        ss=as.character(oiseaux_Frequents_t[i])
        sp_all=c(sp_all,rep(ss,length(dd)))
        abundance_cold=c(abundance_cold,tab_freq[ss,'Cold',dd])
        abundance_warm=c(abundance_warm,tab_freq[ss,'Warm',dd])
        plot(1981:2015,tab_freq[ss,'Cold',dd],ylim=c(0,max(tab_freq[ss,,],na.rm=T)),xlim=c(1981,2015),col="blue",t="o",pch=16,ylab=ss,xlab="")
        lines(1981:2015,tab_freq[ss,'Warm',dd],col="red",t="o",pch=16)
}
data_frame_freq_cold=data.frame(dates,sp_all,abundance_cold)
data_frame_freq_warm=data.frame(dates,sp_all,abundance_warm)
write.csv(data_frame_freq_cold,paste("IN/coldseason_freq_detailed_v2_",end_bio,"_",end_ignore,".txt",sep=""))
write.csv(data_frame_freq_warm,paste("IN/warmseason_freq_detailed_v2_",end_bio,"_",end_ignore,".txt",sep=""))
dev.off()

print('Ducks')
tab_tmp=read.csv(file="./IN/Initial_files/data_ROT20160324.csv",header=TRUE,sep="\t")
sp_all=tab_tmp$Nom_latin
fam=tab_tmp$Famille
sp_duck=c(as.character(unique(sp_all[fam=="Anatidae"])),"Fulica atra")
sp_duck=intersect(sp_duck,DBt$Nom_latin)
sp_duck=setdiff(sp_duck,sp_to_ignore)
tmp=subset(DBt,DBt$Nom_latin %in% sp_duck)
if(biomass){
tmp_bis=tmp
for(i in 1:nrow(tmp_bis)){
  tmp_bis$Nombre[i]=tmp_bis$Nombre[i]*tab_bm$Mean_mass[as.character(tab_bm$Species)==as.character(tmp_bis$Nom_latin[i])]
}
tmp=tmp_bis
}
tab_duck=season2_average(tmp)

dates=rep(dd,length(sp_duck))
sp_all=c()
abundance_cold=c()
abundance_warm=c()
pdf(paste("OUT/closer_look_ducks_v2_",end_bio,"_",end_ignore,".pdf",sep=""))
for (i in 1:length(sp_duck)){
        if(i%%12==1){
        par(mfrow=c(4,3))
        }
        ss=as.character(sp_duck[i])
        sp_all=c(sp_all,rep(ss,length(dd)))
        abundance_cold=c(abundance_cold,tab_duck[ss,'Cold',dd])
        abundance_warm=c(abundance_warm,tab_duck[ss,'Warm',dd])
        plot(1981:2015,tab_duck[ss,'Cold',dd],ylim=c(0,max(tab_duck[ss,,],na.rm=T)),xlim=c(1981,2015),col="blue",t="o",pch=16,ylab=ss,xlab="")
        lines(1981:2015,tab_duck[ss,'Warm',dd],col="red",t="o",pch=16)
}
data_frame_duck_cold=data.frame(dates,sp_all,abundance_cold)
data_frame_duck_warm=data.frame(dates,sp_all,abundance_warm)
write.csv(data_frame_duck_cold,paste("IN/coldseason_ducks_detailed_v2_",end_bio,"_",end_ignore,".txt",sep=""))
write.csv(data_frame_duck_warm,paste("IN/warmseason_ducks_detailed_v2_",end_bio,"_",end_ignore,".txt",sep=""))
dev.off()
