#CP function to sum groups and compute their average
sum_of_species=function(tab, sp,nom){
	require('lubridate')
	tab=subset(tab,tab$Nom_latin %in% sp, c("Nombre","Date"))
	tab_sum=aggregate(tab$Nombre,list(tab$Date),sum)
	names(tab_sum)=c("Date","Nombre")
	tab_sum$Date=as.Date(tab_sum$Date)
	yy=unique(year(tab$Date))
	nombre=c()
	dates=c()
	for(y in yy){
		min_m=1
		max_m=12
		#Do not create artefacts
		if(y==min(yy)){
		min_m=min(month(tab_sum$Date[year(tab_sum$Date)==y]))
		}else if(y==max(yy)){
		max_m=max(month(tab_sum$Date[year(tab_sum$Date)==y]))
		}
		for(m in min_m:max_m){
			dates=c(dates,as.character(paste(y,sprintf("%02d",m),"15",sep="-")))
			if(length(sum(tab_sum$Nombre[month(tab_sum$Date)==m&year(tab_sum$Date)==y]))>0){
				nombre=c(nombre,sum(tab_sum$Nombre[month(tab_sum$Date)==m&year(tab_sum$Date)==y]))
			}else{
				nombre=c(nombre,0)
			}
		}
	}
	tab_sum=cbind(rep(nom,length(nombre)),dates,nombre)
        tab_sum=as.data.frame(tab_sum)
	names(tab_sum)=c("Nom_latin","Date","Nombre")
	return(tab_sum)
}

as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

season2_average=function(tab){
#Change in code from 2018/06/10: before, cold season was of year y was november and december of y-1 and january and february of y. From now on, it's 11 and 12 of y and 1 and 2 of y+1
	require('lubridate')
	Hivernage = c(11,12,1,2)
	Nichage   = c(5,6,7,8)
	yy=unique(year(tab$Date))
	sp=unique(tab$Nom_latin)
	array_mean=array(0,dim=c(length(sp),2,length(yy)),dimnames=list(sp,c("Cold","Warm"),as.character(yy)))
	if(class(tab$Nombre)=="factor"){
        tab$Nombre=as.numeric.factor(tab$Nombre)
	}
	tab_yy=array(NA,dim=c(length(yy),2),dimnames=list(as.character(yy),c('Cold','Warm')))
	for (y in yy){
		length_warm=length(Nichage)
		length_cold=length(Hivernage)
		if(y==min(yy)){
			list_month=month(tab$Date[year(tab$Date)==y])
			list_month_next=month(tab$Date[year(tab$Date)==(y+1)])
			if(min(list_month)>min(Nichage)&max(list_month)<max(Nichage)){
				length_warm=length(intersect(Nichage,list_month))
			}else if(min(list_month)>max(Nichage)){
				length_warm=NA
			}
			length_cold=length(intersect(c(1,2),list_month_next))+length(intersect(c(11,12),list_month))
		}else if (y==max(yy)){
			list_month=month(tab$Date[year(tab$Date)==y])
			length_warm=length(intersect(Nichage,list_month))
			length_cold=length(intersect(c(11,12),list_month))
		}
		tab_yy[as.character(y),]=c(length_cold,length_warm)
	
	}
	for (dd in 1:length(tab$Date)){
        	d=tab$Date[dd]
	        y=as.character(year(as.Date(d)))
        	mois=as.integer(month(as.Date(d)))
		for (s in sp){
			if(dd==1){
			print(s)
			}
			 if(length(tab$Nombre[tab$Date==d&tab$Nom_latin==s])>0){
				if(mois>=min(Nichage)&mois<=max(Nichage)){
			               array_mean[s,"Warm",y]=array_mean[s,"Warm",y]+1/tab_yy[y,'Warm']*tab$Nombre[tab$Date==d&tab$Nom_latin==s]
        			}else if(mois>=Hivernage[1]){
	                		array_mean[s,"Cold",y]=array_mean[s,"Cold",y]+1/tab_yy[y,'Cold']*tab$Nombre[tab$Date==d&tab$Nom_latin==s]
        			}else if(mois>0&mois<=Hivernage[4]){
					if(as.integer(y)>min(yy)){ #because we can do nothing with january and february 1981
			                	array_mean[s,"Cold",as.character(as.integer(y)-1)]=array_mean[s,"Cold",as.character(as.integer(y)-1)]+1/tab_yy[as.character(as.integer(y)-1),'Cold']*tab$Nombre[tab$Date==d&tab$Nom_latin==s]
					}
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



