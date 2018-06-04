#CP 2018 output txt files for sum of species (anas, calidris, waders, ducks, all birds other than waders, héron+aigrette), then compute seasonal average for cold and warm season (anas, calidris, waders, ducks, frequent birds)

rm(list=ls())
graphics.off()

sum_of_species=function(tab, sp,nom){
	tab=subset(tab,tab$Nom_latin %in% sp, c("Nombre","Date"))
	tab_sum=aggregate(tab$Nombre,list(tab$Date),sum)
	tab_sum=cbind(rep(nom,length(tab_sum[[1]])),tab_sum)
        names(tab_sum)=c("Nom_latin","Date","Nombre")
	return(tab_sum)
}

season2_average=function(tab){
	require('lubridate')
	Hivernage = c(11,12,1,2)
	Nichage   = c(5,6,7,8)
	yy=unique(year(tab$Date))
	sp=unique(tab$Nom_latin)
	array_mean=array(0,dim=c(length(sp),2,length(yy)),dimnames=list(sp,c("Cold","Warm"),as.character(yy)))
	
	tab_yy=array(NA,dim=c(length(yy),2),dimnames=list(as.character(yy),c('Cold','Warm')))
	for (y in yy){
		length_warm=length(Nichage)
		length_cold=length(Hivernage)
		if(y==min(yy)){
			list_month=month(tab$Date[year(tab$Date)==y])
			if(min(list_month)>min(Nichage)&max(list_month)<max(Nichage)){
				length_warm=length(intersect(Nichage,list_month))
			}else if(min(list_month)>max(Nichage)){
				length_warm=NA
			}
			length_cold=length(intersect(c(1,2),list_month))
		}else if (y==max(yy)){
			list_month=month(tab$Date[year(tab$Date)==y])
			length_warm=length(intersect(Nichage,list_month))
			list_month_previous=month(tab$Date[year(tab$Date)==(y-1)])
			length_cold=length(intersect(c(11,12),list_month_previous))+length(intersect(c(1,2),list_month))
		}
		tab_yy[as.character(y),]=c(length_cold,length_warm)
	}
	for (dd in 1:length(tab$Date)){
        	d=tab$Date[dd]
	        y=as.character(year(as.Date(d)))
        	mois=as.integer(month(as.Date(d)))
		for (s in sp){
			 if(length(tab$Nombre[tab$Date==d&tab$Nom_latin==s])>0){
				if(mois>=min(Nichage)&mois<=max(Nichage)){
			               array_mean[s,"Warm",y]=array_mean[s,"Warm",y]+1/tab_yy[y,'Warm']*tab$Nombre[tab$Date==d&tab$Nom_latin==s]
        			}else if(mois>=Hivernage[1]){
					if (y<max(yy)){ #Because we can do nothing with 11/2016
	                			array_mean[s,"Cold",as.character(as.integer(y)+1)]=array_mean[s,"Cold",as.character(as.integer(y)+1)]+1/tab_yy[as.character(as.integer(y)+1),'Cold']*tab$Nombre[tab$Date==d&tab$Nom_latin==s]
					}
        			}else if(mois>0&mois<=Hivernage[4]){
			        	        array_mean[s,"Cold",y]=array_mean[s,"Cold",y]+1/tab_yy[y,'Cold']*tab$Nombre[tab$Date==d&tab$Nom_latin==s]
        			}

			}
		}
	}
	for(y in yy){
		if(tab_yy[as.character(y),'Cold']==0){
			array_mean[,'Cold',as.character(y)]=NA
		}
		if(tab_yy[as.character(y),'Warm']==0){
			array_mean[,'Warm',as.character(y)]=NA
		}
	}
	return(array_mean)
}



DBt<-read.csv(file="/home/cpicoche/Documents/Birds/BirdTimeSeries_Teich/IN/DBWithMonthlyPhotoTeich_completed.csv",header=TRUE,sep=",",dec=".")
DBt = subset(DBt,(DBt$Protocol==1 | DBt$Protocol==2) & (DBt$Lieu_dit=="USN00-Réserve ornithologique (générique)" | DBt$Lieu_dit=="USN01-Artigues-Réserve ornithologique") & (DBt$Annee>1980))
sp=unique(DBt$Nom_latin)

print('Anas')
sp_anas=sp[grep("Anas",sp)]
anas_sum=sum_of_species(DBt,sp_anas,"Anas")
anas_season=season2_average(anas_sum)

print('Calidris')
sp_calidris=sp[grep("Calidris",sp)]
calidris_sum=sum_of_species(DBt,sp_calidris,"Calidris")
calidris_season=season2_average(calidris_sum)

print('Waders')
limicoles = c("Recurvirostra avosetta","Limosa limosa","Limosa lapponica","Calidris temminckii","Calidris canutus",
              "Calidris alba","Calidris alpina","Calidris minuta","Calidris maritima" ,"Gallinago gallinago",
              "Tringa flavipes","Tringa nebularia","Tringa erythropus","Tringa ochropus","Tringa totanus",
              "Tringa glareola","Actitis hypoleucos","Philomachus pugnax","Numenius arquata","Numenius phaeopus",
              "Himantopus himantopus","Charadrius hiaticula","Charadrius alexandrinus","Haematopus ostralegus",
              "Burhinus oedicnemus","Charadrius dubius","Phalaropus lobatus","Pluvialis squatarola",
              "Pluvialis apricaria","Arenaria interpres","Vanellus vanellus")
wader_sum=sum_of_species(DBt,limicoles,"Waders")
wader_season=season2_average(wader_sum)

print('Ducks')
tab_tmp=read.csv(file="./IN/Initial_files/data_ROT20160324.csv",header=TRUE,sep="\t")
sp_all=tab_tmp$Nom_latin
fam=tab_tmp$Famille

sp_duck=c(as.character(unique(sp_all[fam=="Anatidae"])),"Fulica atra")
duck_sum=sum_of_species(DBt,sp_duck,"Ducks")
duck_season=season2_average(duck_sum)

print('Other birds than waders') ##This will be the longest one
sp_other=setdiff(sp_all,limicoles)
other_than_wader_sum=sum_of_species(DBt,sp_other,"Not-waders")
other_than_wader_season=season2_average(other_than_wader_sum)

print('Frequent birds')
vec_n_obs_oiseaux_t=rep(0,length(sp_all))
for (i in 1:length(sp_all)){
  vec_n_obs_oiseaux_t[i]=sum(as.character(DBt$Nom_latin)==sp_all[i],na.rm = TRUE)
}
oiseaux_Frequents_t=sp_all[vec_n_obs_oiseaux_t>75]

freq_sum=sum_of_species(DBt,oiseaux_Frequents_t,'Freq')
freq_season=season2_average(freq_sum)

#Small trick for Cormorant (we don't do a sum but it gives the format I want)
print('Cormorant')
cormorant_sum=sum_of_species(DBt,c("Phalacrocorax carbo"),'Cormorant')
cormorant_season=season2_average(cormorant_sum)

print('Finally Heron+Egret')
heregr_sum=sum_of_species(DBt,c("Ardea cinerea","Egretta garzetta"),'HeronEgret')
heregr_season=season2_average(heregr_sum)

all_sum=rbind(anas_sum,calidris_sum,wader_sum,duck_sum,other_than_wader_sum,freq_sum,cormorant_sum,heregr_sum)
write.table(all_sum,'IN/summed_abundances.csv',sep=";")

table_seasonal_cold=matrix(NA,nrow=length(anas_season[,"Cold",]),ncol=8)
rownames(table_seasonal_cold)=dimnames(anas_season)[[3]]
colnames(table_seasonal_cold)=c('Anas','Calidris','Waders','Ducks','Not-waders','Freq','Cormorant','HeronEgret')
table_seasonal_cold[,'Anas']=anas_season[,'Cold',]
table_seasonal_cold[,'Calidris']=calidris_season[,'Cold',]
table_seasonal_cold[,'Waders']=wader_season[,'Cold',]
table_seasonal_cold[,'Ducks']=duck_season[,'Cold',]
table_seasonal_cold[,'Not-waders']=other_than_wader_season[,'Cold',]
table_seasonal_cold[,'Freq']=freq_season[,'Cold',]
table_seasonal_cold[,'Cormorant']=cormorant_season[,'Cold',]
table_seasonal_cold[,'HeronEgret']=heregr_season[,'Cold',]
write.table(table_seasonal_cold,"IN/coldseason_abundances.csv",sep=";")

table_seasonal_warm=matrix(NA,nrow=length(anas_season[,"Warm",]),ncol=8)
rownames(table_seasonal_warm)=dimnames(anas_season)[[3]]
colnames(table_seasonal_warm)=c('Anas','Calidris','Waders','Ducks','Not-waders','Freq','Cormorant','HeronEgret')
table_seasonal_warm[,'Anas']=anas_season[,'Warm',]
table_seasonal_warm[,'Calidris']=calidris_season[,'Warm',]
table_seasonal_warm[,'Waders']=wader_season[,'Warm',]
table_seasonal_warm[,'Ducks']=duck_season[,'Warm',]
table_seasonal_warm[,'Not-waders']=other_than_wader_season[,'Warm',]
table_seasonal_warm[,'Freq']=freq_season[,'Warm',]
table_seasonal_warm[,'Cormorant']=cormorant_season[,'Warm',]
table_seasonal_warm[,'HeronEgret']=heregr_season[,'Warm',]
write.table(table_seasonal_warm,"IN/warmseason_abundances.csv",sep=";")

#Produce data frames and not table
building_line='data_frame_cold=rbind('
for (nem in colnames(table_seasonal_cold)){
	building_line=paste(building_line,'cbind(rep(\"',nem,'\",',dim(table_seasonal_cold)[1],'),table_seasonal_cold[,\"',nem,'\"]),',sep="")
}
building_line=substr(building_line,1,nchar(building_line)-1)
building_line=paste(building_line,')',sep="")
eval(parse(text=building_line))
data_frame_cold=cbind(rownames(data_frame_cold),data_frame_cold)
write.table(data_frame_cold,"IN/coldseason_abundances_asdataframe.csv",sep=";",col.names=c("Date","Species","Abundance"),row.names=F)

#Produce data frames and not table
building_line='data_frame_warm=rbind('
for (nem in colnames(table_seasonal_warm)){
        building_line=paste(building_line,'cbind(rep(\"',nem,'\",',dim(table_seasonal_warm)[1],'),table_seasonal_warm[,\"',nem,'\"]),',sep="")
}
building_line=substr(building_line,1,nchar(building_line)-1)
building_line=paste(building_line,')',sep="")
eval(parse(text=building_line))
data_frame_warm=cbind(rownames(data_frame_warm),data_frame_warm)
write.table(data_frame_warm,"IN/warmseason_abundances_asdataframe.csv",sep=";",col.names=c("Date","Species","Abundance"),row.names=F)


