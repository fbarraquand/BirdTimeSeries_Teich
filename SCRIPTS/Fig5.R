
# Script by CPicoche 2018
#The aim is to compare synchrony indices at different taxonomic levels
rm(list=ls())
graphics.off()
source("SCRIPTS/test_synchrony_Gross.r")
library('mvcwt')
source("SCRIPTS/image_mvcwt_two_panels.r") #Add to change the image function to have a nice Color Bar
library("RColorBrewer")

thresh=0.1
type_correct="BH" #was Bonferroni before
log_b=TRUE
amethod="iaaft"
anrands=100

# Grand Cormoran Phalacrocorax carbo
# Héron cendré Ardea cinerea
# Aigrette garzette Egretta garzetta

db_warm=read.csv("IN/warmseason_abundances_asdataframe_summed.csv",sep=";",header=T)
db_cold=read.csv("IN/coldseason_abundances_asdataframe_summed.csv",sep=";",header=T)

#Right format for synchrony scripts
names(db_warm)=c("dates","sp_data_frame","abundance")
names(db_cold)=c("dates","sp_data_frame","abundance")

if(log_b){
        db_warm$abundance=log(db_warm$abundance+1)
        db_cold$abundance=log(db_cold$abundance+1)
}

db_warm_all=subset(db_warm,sp_data_frame %in% c("Cormorant","HeronEgret")&dates<2016)
db_warm_pre_2006=subset(db_warm,sp_data_frame %in% c("Cormorant","HeronEgret")&dates<=2006)
db_warm_post_2006=subset(db_warm,sp_data_frame %in% c("Cormorant","HeronEgret")&dates>2006&dates<2016)#Because using the NA for 2016 makes results wrong
db_cold_all=subset(db_cold,sp_data_frame %in% c("Cormorant","HeronEgret")&dates<2016)
db_cold_pre_2006=subset(db_cold,sp_data_frame %in% c("Cormorant","HeronEgret")&dates<=2006)
db_cold_post_2006=subset(db_cold,sp_data_frame %in% c("Cormorant","HeronEgret")&dates>2006&dates<2016)

id_c=db_warm$sp_data_frame=="Cormorant"
#plot(db_warm$dates[id_c],db_warm$abundance[id_c],col="red",lty=1,lwd=1.5,t="l",ylim=c(0,1000))
id_c=db_warm$sp_data_frame=="HeronEgret"
#lines(db_warm$dates[id_c],db_warm$abundance[id_c],col="red",lty=2)
id_c=db_cold$sp_data_frame=="Cormorant"
#lines(db_cold$dates[id_c],db_cold$abundance[id_c],col="blue",lty=1,lwd=1.5)
id_c=db_cold$sp_data_frame=="HeronEgret"
#lines(db_cold$dates[id_c],db_cold$abundance[id_c],col="blue",lty=2)
#Compute synchrony values
synch_warm_all=community_sync_Gross(db_warm_all,nrands=anrands,method=amethod)
synch_warm_pre_2006=community_sync_Gross(db_warm_pre_2006,nrands=anrands,method=amethod)
synch_warm_post_2006=community_sync_Gross(db_warm_post_2006,nrands=anrands,method=amethod)

synch_cold_all=community_sync_Gross(db_cold_all,nrands=anrands,method=amethod)
synch_cold_pre_2006=community_sync_Gross(db_cold_pre_2006,nrands=anrands,method=amethod)
synch_cold_post_2006=community_sync_Gross(db_cold_post_2006,nrands=anrands,method=amethod)


#Plot everything
essai_taxo=list(synch_cold_all,synch_cold_pre_2006,synch_cold_post_2006,synch_warm_all,synch_warm_pre_2006,synch_warm_post_2006)
mat=rep(NA,length(essai_taxo))
for(v in 1:length(essai_taxo)){
        mat[v]=essai_taxo[[v]]$pval
}
mat_adj=p.adjust(mat,method=type_correct)
for(v in 1:length(essai_taxo)){
        essai_taxo[[v]]$pval=mat_adj[v]
}


color=rep(c("Black","Lightblue","Darkblue"),2)
pdf("Submission_JAE/Revisions/triad_synchrony_2panels.pdf")
layout(matrix(c(1,1,2,3),nrow=2,ncol=2,byrow=T),widths=c(10,1))

#par(mfrow=c(1,1),mar=c(3,3.5,2,.25),oma=c(1,2,1,.25),mgp=c(3,1,0),xpd=NA)
par(mar=c(4,5,2,1))




plot(0,0,t="n",ylim=c(-1.,1.0),xlim=c(0,7.5),xaxt="n",xlab="",ylab="Synchrony Index",cex.lab=1.5,cex.axis=1.5)
mtext("a)",side=2,line=-2,at=0.96,cex=1.5,outer=T,las=1)
axis(1,at=c(2,6),lab=c("Cold","Warm"),cex.axis=1.5)
for (v in 1:6){
        obs=essai_taxo[[v]]$obs
        print(obs)
#        ic=2*sd(essai_taxo[[v]]$rands[1:100])/sqrt(100)
        p_s=essai_taxo[[v]]$pval
        x=v+(v>3)
#        yHigh=as.numeric(obs)+as.numeric(ic)
#        yLow=as.numeric(obs)-as.numeric(ic)
        points(x,obs,pch=21,col=color[v],bg=color[v],cex=2)
#        arrows(x,yHigh,x,yLow,angle=90,length=0.1,code=3,col=color[v],lwd=2)
        boxplot(essai_taxo[[v]]$rands[1:100],at=x,add=T,boxwex=0.25,range=0,yaxt="n",xaxt="n")
        if (p_s<thresh){
                points(x,as.numeric(obs),pch='*',col="red",cex=2)
                }

        }
lines(c(0,7.5),c(0,0),lty=2,lwd=2)
#ll=c(essai_taxo[[2]]$obs,essai_taxo[[3]]$obs)
#lines(2:3,ll,col="black",lwd=2,lty=1)
#ll=c(essai_taxo[[5]]$obs,essai_taxo[[6]]$obs)
#lines(6:7,ll,col="black",lwd=2,lty=1)

legend("bottomleft",c("All","Pre-2006","Post-2006"),pch=NA,fill=c("black","Lightblue","Darkblue"),pt.cex=2,bty="n",cex=1.5)
############################################ WAVELETS ########################################

db=read.csv("IN/summed_abundances.csv",sep=";",header=T)
db$Date=as.Date(db$Date)

db=subset(db,Nom_latin %in% c("Cormorant","HeronEgret"))
db_tmp=db
names(db_tmp)=c("sp_data_frame","dates","abundance")
db_av=subset(db_tmp,dates<as.Date("2007-01-01"))
db_ap=subset(db_tmp,dates>=as.Date("2007-01-01"))
db_tmp$dates=as.numeric(db_tmp$dates)
db_av$dates=as.numeric(db_av$dates)
db_ap$dates=as.numeric(db_ap$dates)

synch_all=community_sync_Gross(db_tmp,nrands=100) #-0.09735, p=0.74
synch_av=community_sync_Gross(db_av,nrands=100) #-0.2028886, p=0.34
synch_ap=community_sync_Gross(db_ap,nrands=100) #0.2758426, p=0.11


dates=unique(db$Date)
tab=matrix(0,nrow=length(dates),ncol=2)
colnames(tab)=c("Cormorant","HeronEgret")
rownames(tab)=dates
for(id in 1:length(dates)){
        for (s in c('Cormorant',"HeronEgret")){
                id_d=which(db$Date==dates[id]&db$Nom_latin==s)
                if(length(id_d)>0){
                        tab[id,s]=db$Nombre[id_d]
                }
        }
}

x=(dates-dates[1])/365.25
year_min=1981
#This function computes the Morlet wavelet transform for each bird species separately
print(Sys.time())
mm=mvcwt(x,tab,min.scale=0.2,max.scale=10.0)

#This function computes the wavelet ratio of the whole community (see Keitt's paper in 2008)
mr = wmr.boot(mm, smoothing = 1,reps=anrands)
mr$x=mr$x+year_min #Change the dates to be "human-readable"

#png('OUT/wavelet_triad.png',width=800)
image_mvcwt_two_panels(mr,reset.par=F,cex.axis=4,z.fun="Mod")

#abline(v=2006,lwd=3,col="black") #This is supposed to change in 2006 with water management
print(Sys.time())
dev.off()
#
