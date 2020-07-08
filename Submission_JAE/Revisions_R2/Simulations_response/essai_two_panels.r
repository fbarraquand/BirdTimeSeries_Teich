
# Script by CPicoche 2018
#The aim is to compare synchrony indices at different taxonomic levels
#06/03/2020: Just checking the time necessary for 100 and 1000 randomization
#01/07/2020: Checking the effect of log and standardization

rm(list=ls())
graphics.off()
source("../../../SCRIPTS/test_synchrony_Gross.r")
library('mvcwt')
source("../../../SCRIPTS/image_mvcwt_two_panels.r") #Add to change the image function to have a nice Color Bar
library("RColorBrewer")
library(lubridate)

set.seed(42)

thresh=0.1
type_correct="BH" #was Bonferroni before
log_b=F
amethod="iaaft"
anrands=50
norm=c(FALSE,TRUE)

end_of_file_seq=c("4sp_pair","40sp_pair","4sp_alpha15_beta15","40sp_alpha15_beta15","4sp_alpha2_beta4","40sp_alpha2_beta4","4sp_alpha4_beta2","40sp_alpha4_beta2")
explain=c("rho=-0.8","rho=-0.8","quasi-normal","quasi-normal","compensation","compensation","synchrony","synchrony")

for(e in 1:length(end_of_file_seq)){
#for(e in 1:1){
end_of_file=end_of_file_seq[e]
#############################
tab_bm=read.table(paste("MockData_SAD_",end_of_file,".csv",sep=""),sep=";",dec=".",header=T)

rep=1

tab_bm=subset(tab_bm,Rep==rep)

tab_bm$Species=as.character(tab_bm$Species)
tab_bm$Time=as.numeric(as.character(tab_bm$Time))
tab_bm$Abundance=as.numeric(as.character(tab_bm$Abundance))

dates_tmp=unique(tab_bm$Time)

#Convert 1,2,3 to month and year

dates=rep(NA,length(dates_tmp))
yy=0
for(d in 1:length(dates)){
        if(dates_tmp[d]%%12==1){
                yy=yy+1
        }
        mm=dates_tmp[d]%%12
        if(mm==0){mm=12}
        tmp=as.Date(paste(yy,mm,"15",sep="-"))
        dates[d]=tmp
}

dates=as.Date(dates,origin="01-01-01")
sp=unique(tab_bm$Species)

tab_species=matrix(0,nrow=length(dates),ncol=length(sp))
colnames(tab_species)=sp

for(id in 1:length(dates)){
        for (s in sp){
                id_d=which(tab_bm$Time==dates_tmp[id]&tab_bm$Species==s)
                if(length(id_d)>0){
                        tab_species[id,s]=tab_bm$Abundance[id_d]
                }
        }
}

yy_seq=unique(year(dates))
tab_annual_warm=as.data.frame(matrix(NA,nrow=length(yy_seq)*length(sp),ncol=3))
tab_annual_cold=as.data.frame(matrix(NA,nrow=length(yy_seq)*length(sp),ncol=3))
names(tab_annual_warm)=names(tab_annual_cold)=c("dates","sp_data_frame","abundance")


#cold month 1:4 (we're not adding difficulties here)
#warm month 7:10

for(normalize in norm){
if(normalize){
end_nor="scaled"
}else{
end_nor="NOTscaled"
}
tab_species_tmp=tab_species
for(s in sp){
if(normalize){
        tab_species_tmp[,s]=scale(tab_species[,s])
}

}

x=0
for(s in sp){
for(y in 1:length(yy_seq)){
	x=x+1
	subset=tab_species_tmp[year(dates)==yy_seq[y],s]
	tab_annual_warm[x,"dates"]=yy_seq[y]
	tab_annual_warm[x,"sp_data_frame"]=s	
	tab_annual_warm[x,"abundance"]=mean(subset[7:10])
	tab_annual_cold[x,"dates"]=yy_seq[y]
	tab_annual_cold[x,"sp_data_frame"]=s	
	tab_annual_cold[x,"abundance"]=mean(subset[1:4])
}
}

print("Before Gross")
print(Sys.time())
synch_warm_all=community_sync_Gross(tab_annual_warm,nrands=anrands,method=amethod)

synch_cold_all=community_sync_Gross(tab_annual_cold,nrands=anrands,method=amethod)

print("After Gross")
print(Sys.time())

#if(1==0){ ##remove that after
#Plot everything
essai_taxo=list(synch_cold_all,synch_warm_all)
mat=rep(NA,length(essai_taxo))
for(v in 1:length(essai_taxo)){
        mat[v]=essai_taxo[[v]]$pval
}
mat_adj=p.adjust(mat,method=type_correct)
for(v in 1:length(essai_taxo)){
        essai_taxo[[v]]$pval=mat_adj[v]
}


color=rep(c("Black","Lightblue","Dodgerblue2"),2)

pdf(paste("wavelet_simu_",end_nor,"_",end_of_file,"_two_panels.pdf",sep=""),width=10,height=10)

layout(matrix(c(1,1,2,3),nrow=2,ncol=2,byrow=T),widths=c(10,2))

#par(mfrow=c(1,1),mar=c(3,3.5,2,.25),oma=c(1,2,1,.25),mgp=c(3,1,0),xpd=NA)
par(mar=c(4,5,2,3))

plot(0,0,t="n",ylim=c(-1.,1.0),xlim=c(0,3),xaxt="n",xlab="",ylab=expression(eta),cex.lab=1.5,cex.axis=1.5,las=1,main=explain[e])
mtext("a)",side=2,line=-2,at=0.96,cex=1.5,outer=T,las=1)
axis(1,at=c(1,2),lab=c("Cold","Warm"),cex.axis=1.5)
for (v in 1:2){
        obs=essai_taxo[[v]]$obs
        print(obs)
        p_s=essai_taxo[[v]]$pval
        points(v,obs,pch=21,col=color[v],bg=color[v],cex=2)
        boxplot(essai_taxo[[v]]$rands[1:100],at=v,add=T,boxwex=0.25,range=0,yaxt="n",xaxt="n")
        if (p_s<thresh){
                points(v,as.numeric(obs),pch='*',col="red",cex=2)
                }

        }
lines(c(0,7.5),c(0,0),lty=2,lwd=2)


x=(dates-dates[1])/365.25

#This function computes the Morlet wavelet transform for each bird species separately
print(Sys.time())
mm=mvcwt(x,tab_species_tmp,min.scale=0.2,max.scale=10.0)

#This function computes the wavelet ratio of the whole community (see Keitt's paper in 2008)
mr = wmr.boot(mm, smoothing = 1,reps=anrands)

#png('OUT/Figure3.png',width=800)
#  image_mvcwt_for_colormaps(mr,reset.par=F,cex.axis=4,z.fun="Mod",main=explain[e])
  image_mvcwt_two_panels(mr,reset.par=F,cex.axis=4,z.fun="Mod")

print(Sys.time())
dev.off()
}
}

#} #1==0, remove that after
#
