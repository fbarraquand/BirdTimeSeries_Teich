#CP 2018 : To avoid problems with synchrony values. 

#We want the synchrony values inside the groups of Anas, Calidris, waders and frequent birds
#I'm relaunching this to (a) have 'readable' tables and values of synchronys ignoring the values before 1981 and b) have all birds counted at the same time

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
              "Tringa nebularia","Tringa erythropus","Tringa ochropus","Tringa totanus",
              "Tringa glareola","Actitis hypoleucos","Philomachus pugnax","Numenius arquata","Numenius phaeopus",
              "Himantopus himantopus","Charadrius hiaticula","Charadrius alexandrinus","Haematopus ostralegus",
              "Burhinus oedicnemus","Charadrius dubius","Phalaropus lobatus","Pluvialis squatarola",
              "Pluvialis apricaria","Arenaria interpres","Vanellus vanellus")
#Tringa flavipes never appers with Protocol>0
tmp=subset(DBt,DBt$Nom_latin %in% limicoles)
tab_waders=season2_average(tmp)

print('Frequent birds')
vec_n_obs_oiseaux_t=rep(0,length(sp))
for (i in 1:length(sp)){
  vec_n_obs_oiseaux_t[i]=sum(as.character(DBt$Nom_latin)==sp[i],na.rm = TRUE)
}
oiseaux_Frequents_t=sp[vec_n_obs_oiseaux_t>75]
tmp=subset(DBt,DBt$Nom_latin %in% oiseaux_Frequents_t)
tab_freq=season2_average(tmp)

dd=as.character(1981:2015)
	
dates=rep(dd,length(sp_anas))
sp_all=c()
abundance_cold=c()
abundance_warm=c()
pdf("OUT/closer_look_anas.pdf")
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
write.csv(data_frame_anas_cold,"IN/coldseason_anas_detailed.txt")
write.csv(data_frame_anas_warm,"IN/warmseason_anas_detailed.txt")
dev.off()

dates=rep(dd,length(sp_calidris))
sp_all=c()
abundance_cold=c()
abundance_warm=c()
pdf("OUT/closer_look_calidris.pdf")
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
write.csv(data_frame_calidris_cold,"IN/coldseason_calidris_detailed.txt")
write.csv(data_frame_calidris_warm,"IN/warmseason_calidris_detailed.txt")
dev.off()

dates=rep(dd,length(limicoles))
sp_all=c()
abundance_cold=c()
abundance_warm=c()
pdf("OUT/closer_look_waders.pdf")
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
write.csv(data_frame_waders_cold,"IN/coldseason_waders_detailed.txt")
write.csv(data_frame_waders_warm,"IN/warmseason_waders_detailed.txt")
dev.off()

dates=rep(dd,length(oiseaux_Frequents_t))
sp_all=c()
abundance_cold=c()
abundance_warm=c()
pdf("OUT/closer_look_freq.pdf")
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
write.csv(data_frame_freq_cold,"IN/coldseason_freq_detailed.txt")
write.csv(data_frame_freq_warm,"IN/warmseason_freq_detailed.txt")
dev.off()

print('Ducks')
tab_tmp=read.csv(file="./IN/Initial_files/data_ROT20160324.csv",header=TRUE,sep="\t")
sp_all=tab_tmp$Nom_latin
fam=tab_tmp$Famille
sp_duck=c(as.character(unique(sp_all[fam=="Anatidae"])),"Fulica atra")
sp_duck=intersect(sp_duck,DBt$Nom_latin)
tmp=subset(DBt,DBt$Nom_latin %in% sp_duck)
tab_duck=season2_average(tmp)

dates=rep(dd,length(sp_duck))
sp_all=c()
abundance_cold=c()
abundance_warm=c()
pdf("OUT/closer_look_ducks.pdf")
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
write.csv(data_frame_freq_cold,"IN/coldseason_duck_detailed.txt")
write.csv(data_frame_freq_warm,"IN/warmseason_duck_detailed.txt")
dev.off()



