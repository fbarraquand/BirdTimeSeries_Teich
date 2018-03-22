# Coding from the first script by C. Aluome and F. Barraquand to analyze data after creating the monthly snapshot of bird abundance. 
#Focusing on seasonal values

library("codyn")
library("lubridate")
#rm(list=ls())
graphics.off()

##################################################################################
# -- importation des données du Teich
DBt<-read.csv(file="/home/cpicoche/Documents/Data_to_treat/TRANSFERT_LIMICOLES/IN/DBWithMonthlyPhotoTeich.csv",header=TRUE,sep=",",dec=".")
DBt = subset(DBt,(DBt$Protocol==1 | DBt$Protocol==2) & DBt$Lieu_dit=="USN00-Réserve ornithologique (générique)")
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
oiseaux_Frequents_t=SpeciesL[vec_n_obs_oiseaux_t>20] #125/279
oiseaux_Frequents_t=SpeciesL[vec_n_obs_oiseaux_t>50] #75/279
oiseaux_Frequents_t=SpeciesL[vec_n_obs_oiseaux_t>75] #60/279
oiseaux_Frequents_t_F = SpeciesF[vec_n_obs_oiseaux_t>75] #60/279 #même chose mais pour les noms en français.

#############################################################################
##Définition des saisons
winter = c(1,2,3) #numero du mois qui compose chaque saison
spring = c(4,5,6)
summer = c(7,8,9)
autumn = c(10,11,12)

Hivernage = c(10,11,12,1,2)
Nichage   = c(4,5,6,7)

#Two ways of defining values per season : mean abundance or max abundance
#We can also use the three or fGour series as replicates!

sp=oiseaux_Frequents_t
#sp=SpeciesL
yy=unique(year(DBt$Date))
DBt$Nombre=as.numeric(DBt$Nombre)

array_mean_seasonal=array(0,dim=c(length(sp),4,length(yy)),dimnames=list(sp,c("Winter","Spring","Summer","Autumn"),as.character(yy)))
array_mean_hivnich=array(0,dim=c(length(sp),2,length(yy)),dimnames=list(sp,c("Wintering","Breeding"),as.character(yy)))
array_max_seasonal=array(0,dim=c(length(sp),4,length(yy)),dimnames=list(sp,c("Winter","Spring","Summer","Autumn"),as.character(yy)))
array_max_hivnich=array(0,dim=c(length(sp),2,length(yy)),dimnames=list(sp,c("Wintering","Breeding"),as.character(yy)))
for (dd in 1:length(DBt$Date)){
	d=DBt$Date[dd]
	y=as.character(year(as.Date(d)))
	mois=as.integer(month(as.Date(d)))
	tmp_DBT=subset(DBt,Date==d)
	sp1=intersect(sp,tmp_DBT$Nom_latin)
	for (s in sp1){
	if(length(DBt$Nombre[DBt$Date==d&DBt$Nom_latin==s])>0){
	#Season
	if(mois>0&mois<4){
		array_mean_seasonal[s,"Winter",y]=array_mean_seasonal[s,"Winter",y]+1/3*DBt$Nombre[DBt$Date==d&DBt$Nom_latin==s]
		array_max_seasonal[s,"Winter",y]=max(array_max_seasonal[s,"Winter",y],DBt$Nombre[DBt$Date==d&DBt$Nom_latin==s])
	}else if(mois>3&mois<7){
		array_mean_seasonal[s,"Spring",y]=array_mean_seasonal[s,"Spring",y]+1/3*DBt$Nombre[DBt$Date==d&DBt$Nom_latin==s]
		array_max_seasonal[s,"Spring",y]=max(array_max_seasonal[s,"Spring",y],DBt$Nombre[DBt$Date==d&DBt$Nom_latin==s])
	}else if (mois>6&mois<10){
		array_mean_seasonal[s,"Summer",y]=array_mean_seasonal[s,"Summer",y]+1/3*DBt$Nombre[DBt$Date==d&DBt$Nom_latin==s]
		array_max_seasonal[s,"Summer",y]=max(array_max_seasonal[s,"Summer",y],DBt$Nombre[DBt$Date==d&DBt$Nom_latin==s])
	}else{
		array_mean_seasonal[s,"Autumn",y]=array_mean_seasonal[s,"Autumn",y]+1/3*DBt$Nombre[DBt$Date==d&DBt$Nom_latin==s]
		array_max_seasonal[s,"Autumn",y]=max(array_max_seasonal[s,"Autumn",y],DBt$Nombre[DBt$Date==d&DBt$Nom_latin==s])
	}
	#Wintering
	if(mois>3&mois<8){
		array_mean_hivnich[s,"Breeding",y]=array_mean_hivnich[s,"Breeding",y]+1/4*DBt$Nombre[DBt$Date==d&DBt$Nom_latin==s]
		array_max_hivnich[s,"Breeding",y]=max(array_max_hivnich[s,"Breeding",y],DBt$Nombre[DBt$Date==d&DBt$Nom_latin==s])
	}else if(mois>9){
		array_mean_hivnich[s,"Wintering",as.character(as.integer(y)+1)]=array_mean_hivnich[s,"Wintering",as.character(as.integer(y)+1)]+1/5*DBt$Nombre[DBt$Date==d&DBt$Nom_latin==s]
		array_max_hivnich[s,"Wintering",as.character(as.integer(y)+1)]=max(array_max_hivnich[s,"Wintering",as.character(as.integer(y)+1)],DBt$Nombre[DBt$Date==d&DBt$Nom_latin==s])
	}else if(mois>0&mois<3){
		array_mean_hivnich[s,"Wintering",y]=array_mean_hivnich[s,"Wintering",y]+1/5*DBt$Nombre[DBt$Date==d&DBt$Nom_latin==s]
		array_max_hivnich[s,"Wintering",y]=max(array_max_hivnich[s,"Wintering",y],DBt$Nombre[DBt$Date==d&DBt$Nom_latin==s])
	}
	}
}
}

#}#1==0
array_tmp_seasonal=array_mean_seasonal
array_tmp_hivnich=array_mean_hivnich

sp_data_frame=c()
abundance_winter=c()
abundance_spring=c()
abundance_summer=c()
abundance_autumn=c()
abundance_hivernage=c()
abundance_nichage=c()
dates=c()
yy1=yy[yy>2006]
for (i in 1:length(sp)){
	sp_data_frame=c(sp_data_frame,rep(sp[i],length(yy1)))
	abundance_winter=c(abundance_winter,array_tmp_seasonal[sp[i],"Winter",which(yy>2006)])
	abundance_spring=c(abundance_spring,array_tmp_seasonal[sp[i],"Spring",which(yy>2006)])
	abundance_summer=c(abundance_summer,array_tmp_seasonal[sp[i],"Summer",which(yy>2006)])
	abundance_autumn=c(abundance_autumn,array_tmp_seasonal[sp[i],"Autumn",which(yy>2006)])
	abundance_hivernage=c(abundance_hivernage,array_tmp_hivnich[sp[i],"Wintering",which(yy>2006)])
	abundance_nichage=c(abundance_nichage,array_tmp_hivnich[sp[i],"Breeding",which(yy>2006)])
	dates=c(dates,paste(yy1,"-06-01",sep=""))
	}
dates=as.numeric(as.Date(dates))
winter_post_2006=data.frame(sp_data_frame,dates,abundance_winter)
tmp_sum=by(winter_post_2006$abundance_winter,winter_post_2006$dates,sum)
var_winter_com_post_2006=var(tmp_sum)
var_winter_post_2006=by(winter_post_2006$abundance_winter,winter_post_2006$sp_data_frame,sd)

spring_post_2006=data.frame(sp_data_frame,dates,abundance_spring)
tmp_sum=by(spring_post_2006$abundance_spring,spring_post_2006$dates,sum)
var_spring_com_post_2006=var(tmp_sum)
var_spring_post_2006=by(spring_post_2006$abundance_spring,spring_post_2006$sp_data_frame,sd)

summer_post_2006=data.frame(sp_data_frame,dates,abundance_summer)
tmp_sum=by(summer_post_2006$abundance_summer,summer_post_2006$dates,sum)
var_summer_com_post_2006=var(tmp_sum)
var_summer_post_2006=by(summer_post_2006$abundance_summer,summer_post_2006$sp_data_frame,sd)

autumn_post_2006=data.frame(sp_data_frame,dates,abundance_autumn)
tmp_sum=by(autumn_post_2006$abundance_autumn,autumn_post_2006$dates,sum)
var_autumn_com_post_2006=var(tmp_sum)
var_autumn_post_2006=by(autumn_post_2006$abundance_autumn,autumn_post_2006$sp_data_frame,sd)

wintering_post_2006=data.frame(sp_data_frame,dates,abundance_hivernage)
tmp_sum=by(wintering_post_2006$abundance_hivernage,wintering_post_2006$dates,sum)
var_wintering_com_post_2006=var(tmp_sum)
var_wintering_post_2006=by(wintering_post_2006$abundance_hivernage,wintering_post_2006$sp_data_frame,sd)

breeding_post_2006=data.frame(sp_data_frame,dates,abundance_nichage)
tmp_sum=by(breeding_post_2006$abundance_nichage,breeding_post_2006$dates,sum)
var_breeding_com_post_2006=var(tmp_sum)
var_breeding_post_2006=by(breeding_post_2006$abundance_nichage,breeding_post_2006$sp_data_frame,sd)

sp_data_frame=c()
abundance_winter=c()
abundance_spring=c()
abundance_summer=c()
abundance_autumn=c()
abundance_hivernage=c()
abundance_nichage=c()
dates=c()
yy1=yy[yy<=2006]
for (i in 1:length(sp)){
        sp_data_frame=c(sp_data_frame,rep(sp[i],length(yy1)))
        abundance_winter=c(abundance_winter,array_tmp_seasonal[sp[i],"Winter",which(yy<=2006)])
        abundance_spring=c(abundance_spring,array_tmp_seasonal[sp[i],"Spring",which(yy<=2006)])
        abundance_summer=c(abundance_summer,array_tmp_seasonal[sp[i],"Summer",which(yy<=2006)])
        abundance_autumn=c(abundance_autumn,array_tmp_seasonal[sp[i],"Autumn",which(yy<=2006)])
	abundance_hivernage=c(abundance_hivernage,array_tmp_hivnich[sp[i],"Wintering",which(yy<=2006)])
	abundance_nichage=c(abundance_nichage,array_tmp_hivnich[sp[i],"Breeding",which(yy<=2006)])
        dates=c(dates,paste(yy1,"-06-01",sep=""))
        }
dates=as.numeric(as.Date(dates))
winter_pre_2006=data.frame(sp_data_frame,dates,abundance_winter)
tmp_sum=by(winter_pre_2006$abundance_winter,winter_pre_2006$dates,sum)
var_winter_com_pre_2006=var(tmp_sum)
var_winter_pre_2006=by(winter_pre_2006$abundance_winter,winter_pre_2006$sp_data_frame,sd)

spring_pre_2006=data.frame(sp_data_frame,dates,abundance_spring)
tmp_sum=by(spring_pre_2006$abundance_spring,spring_pre_2006$dates,sum)
var_spring_com_pre_2006=var(tmp_sum)
var_spring_pre_2006=by(spring_pre_2006$abundance_spring,spring_pre_2006$sp_data_frame,sd)

summer_pre_2006=data.frame(sp_data_frame,dates,abundance_summer)
tmp_sum=by(summer_pre_2006$abundance_summer,summer_pre_2006$dates,sum)
var_summer_com_pre_2006=var(tmp_sum)
var_summer_pre_2006=by(summer_pre_2006$abundance_summer,summer_pre_2006$sp_data_frame,sd)

autumn_pre_2006=data.frame(sp_data_frame,dates,abundance_autumn)
tmp_sum=by(autumn_pre_2006$abundance_autumn,autumn_pre_2006$dates,sum)
var_autumn_com_pre_2006=var(tmp_sum)
var_autumn_pre_2006=by(autumn_pre_2006$abundance_autumn,autumn_pre_2006$sp_data_frame,sd)

wintering_pre_2006=data.frame(sp_data_frame,dates,abundance_hivernage)
tmp_sum=by(wintering_pre_2006$abundance_hivernage,wintering_pre_2006$dates,sum)
var_wintering_com_pre_2006=var(tmp_sum)
var_wintering_pre_2006=by(wintering_pre_2006$abundance_hivernage,wintering_pre_2006$sp_data_frame,sd)

breeding_pre_2006=data.frame(sp_data_frame,dates,abundance_nichage)
tmp_sum=by(breeding_pre_2006$abundance_nichage,breeding_pre_2006$dates,sum)
var_breeding_com_pre_2006=var(tmp_sum)
var_breeding_pre_2006=by(breeding_pre_2006$abundance_nichage,breeding_pre_2006$sp_data_frame,sd)

sp_data_frame=c()
abundance_winter=c()
abundance_spring=c()
abundance_summer=c()
abundance_autumn=c()
abundance_hivernage=c()
abundance_nichage=c()
dates=c()
yy1=yy
for (i in 1:length(sp)){
        sp_data_frame=c(sp_data_frame,rep(sp[i],length(yy1)))
        abundance_winter=c(abundance_winter,array_tmp_seasonal[sp[i],"Winter",])
        abundance_spring=c(abundance_spring,array_tmp_seasonal[sp[i],"Spring",])
        abundance_summer=c(abundance_summer,array_tmp_seasonal[sp[i],"Summer",])
        abundance_autumn=c(abundance_autumn,array_tmp_seasonal[sp[i],"Autumn",])
	abundance_hivernage=c(abundance_hivernage,array_tmp_hivnich[sp[i],"Wintering",])
	abundance_nichage=c(abundance_nichage,array_tmp_hivnich[sp[i],"Breeding",])
        dates=c(dates,paste(yy1,"-06-01",sep=""))
        }
dates=as.numeric(as.Date(dates))

winter_all=data.frame(sp_data_frame,dates,abundance_winter)
tmp_sum=by(winter_all$abundance_winter,winter_all$dates,sum)
var_winter_com_all=var(tmp_sum)
var_winter_all=by(winter_all$abundance_winter,winter_all$sp_data_frame,sd)

spring_all=data.frame(sp_data_frame,dates,abundance_spring)
tmp_sum=by(spring_all$abundance_spring,spring_all$dates,sum)
var_spring_com_all=var(tmp_sum)
var_spring_all=by(spring_all$abundance_spring,spring_all$sp_data_frame,sd)

summer_all=data.frame(sp_data_frame,dates,abundance_summer)
var_summer_all=by(summer_all$abundance_summer,summer_all$sp_data_frame,sd)
tmp_sum=by(summer_all$abundance_summer,summer_all$dates,sum)
var_summer_com_all=var(tmp_sum)

autumn_all=data.frame(sp_data_frame,dates,abundance_autumn)
var_autumn_all=by(autumn_all$abundance_autumn,autumn_all$sp_data_frame,sd)
tmp_sum=by(autumn_all$abundance_autumn,autumn_all$dates,sum)
var_autumn_com_all=var(tmp_sum)

wintering_all=data.frame(sp_data_frame,dates,abundance_hivernage)
tmp_sum=by(wintering_all$abundance_hivernage,wintering_all$dates,sum)
var_wintering_com_all=var(tmp_sum)
var_wintering_all=by(wintering_all$abundance_hivernage,wintering_all$sp_data_frame,sd)

breeding_all=data.frame(sp_data_frame,dates,abundance_nichage)
tmp_sum=by(breeding_all$abundance_nichage,breeding_all$dates,sum)
var_breeding_com_all=var(tmp_sum)
var_breeding_all=by(breeding_all$abundance_nichage,breeding_all$sp_data_frame,sd)


swinter_all=synchrony(winter_all,time.var="dates",species.var="sp_data_frame",abundance.var="abundance_winter") 
sspring_all=synchrony(spring_all,time.var="dates",species.var="sp_data_frame",abundance.var="abundance_spring") 
ssummer_all=synchrony(summer_all,time.var="dates",species.var="sp_data_frame",abundance.var="abundance_summer")
sautumn_all=synchrony(autumn_all,time.var="dates",species.var="sp_data_frame",abundance.var="abundance_autumn") 
swintering_all=synchrony(wintering_all,time.var="dates",species.var="sp_data_frame",abundance.var="abundance_hivernage")
sbreeding_all=synchrony(breeding_all,time.var="dates",species.var="sp_data_frame",abundance.var="abundance_nichage") 

swinter_pre_2006=synchrony(winter_pre_2006,time.var="dates",species.var="sp_data_frame",abundance.var="abundance_winter") 
sspring_pre_2006=synchrony(spring_pre_2006,time.var="dates",species.var="sp_data_frame",abundance.var="abundance_spring") 
ssummer_pre_2006=synchrony(summer_pre_2006,time.var="dates",species.var="sp_data_frame",abundance.var="abundance_summer")
sautumn_pre_2006=synchrony(autumn_pre_2006,time.var="dates",species.var="sp_data_frame",abundance.var="abundance_autumn") 
swintering_pre_2006=synchrony(wintering_pre_2006,time.var="dates",species.var="sp_data_frame",abundance.var="abundance_hivernage")
sbreeding_pre_2006=synchrony(breeding_pre_2006,time.var="dates",species.var="sp_data_frame",abundance.var="abundance_nichage") 

swinter_post_2006=synchrony(winter_post_2006,time.var="dates",species.var="sp_data_frame",abundance.var="abundance_winter") 
sspring_post_2006=synchrony(spring_post_2006,time.var="dates",species.var="sp_data_frame",abundance.var="abundance_spring") 
ssummer_post_2006=synchrony(summer_post_2006,time.var="dates",species.var="sp_data_frame",abundance.var="abundance_summer")
sautumn_post_2006=synchrony(autumn_post_2006,time.var="dates",species.var="sp_data_frame",abundance.var="abundance_autumn") 
swintering_post_2006=synchrony(wintering_post_2006,time.var="dates",species.var="sp_data_frame",abundance.var="abundance_hivernage")
sbreeding_post_2006=synchrony(breeding_post_2006,time.var="dates",species.var="sp_data_frame",abundance.var="abundance_nichage") 

par(mfrow=c(2,1))
plot(1:4,c(swinter_all,sspring_all,ssummer_all,sautumn_all),t="o",col="black",ylim=c(0,1))
lines(1:4,c(swinter_pre_2006,sspring_pre_2006,ssummer_pre_2006,sautumn_pre_2006),t="o",col="red")
lines(1:4,c(swinter_post_2006,sspring_post_2006,ssummer_post_2006,sautumn_post_2006),t="o",col="blue")
legend("topleft",c("All","Pre-2006","Post-2006"),col=c("black","red","blue"),lty=1,pch=16)

plot(1:2,c(swintering_all,sbreeding_all),t="o",col="black",ylim=c(0,1))
lines(1:2,c(swintering_pre_2006,sbreeding_pre_2006),t="o",col="red")
lines(1:2,c(swintering_post_2006,sbreeding_post_2006),t="o",col="blue")
legend("topleft",c("All","Pre-2006","Post-2006"),col=c("black","red","blue"),lty=1,pch=16)

