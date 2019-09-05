#CP 30/05/18 Ducks
rm(list=ls())
graphics.off()
library("lubridate")

#Define Ducks : All Anatidae + Fulica atra
tab_tmp=read.csv(file="./IN/Initial_files/data_ROT20160324.csv",header=TRUE,sep="\t")
sp=tab_tmp$Nom_latin
fam=tab_tmp$Famille

sp=c(as.character(unique(sp[fam=="Anatidae"])),"Fulica atra") #Other Rallidae do not dive to feed, I considered this feature as a definition of ducks

Hivernage = c(11,12,1,2)
Nichage   = c(5,6,7,8)

#Upload our DB
DBt<-read.csv(file="/home/cpicoche/Documents/Birds/BirdTimeSeries_Teich/IN/DBWithMonthlyPhotoTeich_completed.csv",header=TRUE,sep=",",dec=".")
DBt = subset(DBt,(DBt$Protocol==1 | DBt$Protocol==2) & (DBt$Lieu_dit=="USN00-Réserve ornithologique (générique)" | DBt$Lieu_dit=="USN01-Artigues-Réserve ornithologique")  )
DBt$Date=as.Date(as.character(DBt$Date))

DBt<-DBt[!colnames(DBt) %in% c("X","Ref")]
DBt = unique.matrix(DBt)

yy=unique(year(DBt$Date))
DBt$Nombre=as.numeric(DBt$Nombre)

array_mean_hivnich=array(0,dim=c(length(sp),2,length(yy)),dimnames=list(sp,c("Wintering","Breeding"),as.character(yy)))
for (dd in 1:length(DBt$Date)){
        d=DBt$Date[dd]
        y=as.character(year(as.Date(d)))
        mois=as.integer(month(as.Date(d)))
        tmp_DBT=subset(DBt,Date==d)
        sp1=intersect(sp,tmp_DBT$Nom_latin)
        for (s in sp1){
        if(length(DBt$Nombre[DBt$Date==d&DBt$Nom_latin==s])>0){
        if(mois>4&mois<9){
                array_mean_hivnich[s,"Breeding",y]=array_mean_hivnich[s,"Breeding",y]+1/length(Nichage)*DBt$Nombre[DBt$Date==d&DBt$Nom_latin==s]
        }else if(mois>10){
                array_mean_hivnich[s,"Wintering",as.character(as.integer(y)+1)]=array_mean_hivnich[s,"Wintering",as.character(as.integer(y)+1)]+1/length(Hivernage)*DBt$Nombre[DBt$Date==d&DBt$Nom_latin==s]
        }else if(mois>0&mois<3){
                array_mean_hivnich[s,"Wintering",y]=array_mean_hivnich[s,"Wintering",y]+1/length(Hivernage)*DBt$Nombre[DBt$Date==d&DBt$Nom_latin==s]
        }
        }
}
}

array_tmp_hivnich=array_mean_hivnich
sp_data_frame=c()
abundance_hivernage=c()
abundance_nichage=c()
dates=c()
yy1=yy[yy>2006]
for (i in 1:length(sp)){
        sp_data_frame=c(sp_data_frame,rep(sp[i],length(yy1)))
        abundance_hivernage=c(abundance_hivernage,array_tmp_hivnich[sp[i],"Wintering",which(yy>2006)])
        abundance_nichage=c(abundance_nichage,array_tmp_hivnich[sp[i],"Breeding",which(yy>2006)])
        dates=c(dates,paste(yy1,"-06-01",sep=""))
        }
dates=as.numeric(as.Date(dates))

wintering_post_2006=data.frame(sp_data_frame,dates,abundance_hivernage)

breeding_post_2006=data.frame(sp_data_frame,dates,abundance_nichage)

sp_data_frame=c()
abundance_hivernage=c()
abundance_nichage=c()
dates=c()
yy1=yy[yy<=2006]
for (i in 1:length(sp)){
        sp_data_frame=c(sp_data_frame,rep(sp[i],length(yy1)))
        abundance_hivernage=c(abundance_hivernage,array_tmp_hivnich[sp[i],"Wintering",which(yy<=2006)])
        abundance_nichage=c(abundance_nichage,array_tmp_hivnich[sp[i],"Breeding",which(yy<=2006)])
        dates=c(dates,paste(yy1,"-06-01",sep=""))
        }
dates=as.numeric(as.Date(dates))
wintering_pre_2006=data.frame(sp_data_frame,dates,abundance_hivernage)

breeding_pre_2006=data.frame(sp_data_frame,dates,abundance_nichage)

sp_data_frame=c()
abundance_hivernage=c()
abundance_nichage=c()
dates=c()
yy1=yy
for (i in 1:length(sp)){
        sp_data_frame=c(sp_data_frame,rep(sp[i],length(yy1)))
        abundance_hivernage=c(abundance_hivernage,array_tmp_hivnich[sp[i],"Wintering",])
        abundance_nichage=c(abundance_nichage,array_tmp_hivnich[sp[i],"Breeding",])
        dates=c(dates,paste(yy1,"-06-01",sep=""))
        }
dates=as.numeric(as.Date(dates))
wintering_all=data.frame(sp_data_frame,dates,abundance_hivernage)

breeding_all=data.frame(sp_data_frame,dates,abundance_nichage)

save(wintering_all,breeding_all,wintering_pre_2006,breeding_pre_2006,wintering_post_2006,breeding_post_2006,file=paste("OUT/data_ducks.RData",sep=""))

