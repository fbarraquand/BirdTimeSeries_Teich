# Script by CPicoche 2018

rm(list=ls())
graphics.off()
source("SCRIPTS/test_synchrony_Gross.r")
library('mvcwt')
source("SCRIPTS/wmr_boot.r")
source("SCRIPTS/extract_boot.r")
source("SCRIPTS/image_mvcwt_with_iaaft.r")
source("Submission_JAE/Revisions_R2/Simulations_response/image_mvcwt_for_colormaps.r")
library('RColorBrewer')
source("SCRIPTS/iaaft.R")

set.seed(42)
log_b=F
normalizei_seq=c(T,F)
end_bio="abundances"
nrep=1000

for(normalize in normalize_seq){

db=read.csv(paste("IN/summed_",end_bio,"_v2_wtoutrarespecies.csv",sep=""),sep=";",header=T)
db$Date=as.Date(db$Date)

db=subset(db,Nom_latin %in% c("Cormorant","HeronEgret"))
db_tmp=db
names(db_tmp)=c("sp_data_frame","dates","abundance")
db_av=subset(db_tmp,dates<as.Date("2007-01-01"))
db_ap=subset(db_tmp,dates>=as.Date("2007-01-01"))
db_tmp$dates=as.numeric(db_tmp$dates)
db_av$dates=as.numeric(db_av$dates)
db_ap$dates=as.numeric(db_ap$dates)

dates=unique(db$Date)
tab=matrix(0,nrow=length(dates),ncol=2)
colnames(tab)=c("Cormorant","HeronEgret")
rownames(tab)=dates
for(id in 1:length(dates)){
        for (s in c('Cormorant',"HeronEgret")){
                id_d=which(db$Date==dates[id]&db$Nom_latin==s)
                if(length(id_d)>0){
                        if(log_b){
                                tab[id,s]=log(db$Nombre[id_d]+1)
                        }else{
                                tab[id,s]=db$Nombre[id_d]
                        }
                }
        }
}

if(normalize){
        tab[,"Cormorant"]=scale(tab[,'Cormorant'])
        tab[,"HeronEgret"]=scale(tab[,'HeronEgret'])
}

x=(dates-dates[1])/365.25
year_min=1981
#This function computes the Morlet wavelet transform for each bird species separately
mm=mvcwt(x,tab,min.scale=0.2,max.scale=10.0,nscales=100,loc=regularize(x,nsteps=length(x)/2))
ref_wmr=wmr(mm)
ref_val=ref_wmr$z[,,1]


tab_values_iaaft=array(NA,dim=c(length(mm$x),length(mm$y),nrep+1))
tab_values_iaaft[,,nrep+1]=ref_val
prog.bar = txtProgressBar(min = 0, max = nrep,style = 3)
for(i in 1:nrep){
        setTxtProgressBar(prog.bar, i)
        tab_tmp=tab
        tab_tmp[,"Cormorant"]=iaaft_surrogate(tab[,'Cormorant'])
        tab_tmp[,"HeronEgret"]=iaaft_surrogate(tab[,'HeronEgret'])
        mmtmp=mvcwt(x,tab_tmp,min.scale=0.2,max.scale=10.0,nscales=100,loc=regularize(x,nsteps=length(x)/2))
        wmr_tmp=wmr(mmtmp)
	tab_values_iaaft[,,i]=wmr_tmp$z[,,1]
}

tab_pval=array(NA,dim=c(length(mm$x),length(mm$y),1))
for(i in 1:length(mm$x)){
	for(j in 1:length(mm$y)){
		#tab_pval[i,j,1]= 2*min(sum(tab_values_iaaft[i,j,] >= ref_val[i,j]),sum(tab_values_iaaft[i,j,] < ref_val[i,j]))/(nrep+1)
		tab_pval[i,j,1]= sum(tab_values_iaaft[i,j,] <= ref_val[i,j])/(nrep+1)
		if(tab_pval[i,j,1]>1){stop()}

	}
}
pdf("Submission_JAE/Revisions_R2/iaaft_1000surrogate_triad_nocorrection.pdf")
ref_wmr$z.boot=tab_pval
layout(matrix(c(1,2),nrow=1,ncol=2,byrow=T),widths=c(10,2))
par(mar=c(3,5,2,3))
image_mvcwt_with_iaaft(ref_wmr,adj="None")
dev.off()
