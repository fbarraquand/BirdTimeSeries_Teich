# Script by CPicoche 2018
#The aim is to compare synchrony indices at different taxonomic levels
rm(list=ls())
graphics.off()
source("SCRIPTS/test_synchrony_Gross.r")

thresh=0.05

# Grand Cormoran Phalacrocorax carbo
# Héron cendré Ardea cinerea
# Aigrette garzette Egretta garzetta

db_warm=read.csv("IN/warmseason_abundances_asdataframe.csv",sep=";",header=T)
db_cold=read.csv("IN/coldseason_abundances_asdataframe.csv",sep=";",header=T)

#Right format for synchrony scripts
names(db_warm)=c("dates","sp_data_frame","abundance")
names(db_cold)=c("dates","sp_data_frame","abundance")

db_warm_all=subset(db_warm,sp_data_frame %in% c("Cormorant","HeronEgret"))
db_warm_pre_2006=subset(db_warm,sp_data_frame %in% c("Cormorant","HeronEgret")&dates<=2006)
db_warm_post_2006=subset(db_warm,sp_data_frame %in% c("Cormorant","HeronEgret")&dates>2006)
db_cold_all=subset(db_cold,sp_data_frame %in% c("Cormorant","HeronEgret"))
db_cold_pre_2006=subset(db_cold,sp_data_frame %in% c("Cormorant","HeronEgret")&dates<=2006)
db_cold_post_2006=subset(db_cold,sp_data_frame %in% c("Cormorant","HeronEgret")&dates>2006)

#Compute synchrony values
synch_warm_all=community_sync_Gross(db_warm_all,nrands=100)
synch_warm_pre_2006=community_sync_Gross(db_warm_pre_2006,nrands=100)
synch_warm_post_2006=community_sync_Gross(db_warm_post_2006,nrands=100)

synch_cold_all=community_sync_Gross(db_cold_all,nrands=100)
synch_cold_pre_2006=community_sync_Gross(db_cold_pre_2006,nrands=100)
synch_cold_post_2006=community_sync_Gross(db_cold_post_2006,nrands=100)

#Plot everything
essai_taxo=list(synch_cold_all,synch_cold_pre_2006,synch_cold_post_2006,synch_warm_all,synch_warm_pre_2006,synch_warm_post_2006)
color=rep(c("Black","Lightblue","Darkblue"),2)
pdf("OUT/gross_triad.pdf")
par(mfrow=c(1,1),mar=c(3,2,1,.5))
plot(0,0,t="n",ylim=c(-.5,1.0),xlim=c(0,7.5),xaxt="n",xlab="",ylab="")
axis(1,at=c(2,6),lab=c("Cold","Warm"),cex.axis=2)
for (v in 1:6){
        obs=essai_taxo[[v]]$obs
        ic=2*sd(essai_taxo[[v]]$rands[1:100])/sqrt(100)
        p_s=essai_taxo[[v]]$pval
        x=v+(v>3)
        yHigh=as.numeric(obs)+as.numeric(ic)
        yLow=as.numeric(obs)-as.numeric(ic)
        points(x,obs,pch=21,col=color[v],bg=color[v],cex=2)
        arrows(x,yHigh,x,yLow,angle=90,length=0.1,code=3,col=color[v],lwd=2)
        if (p_s<thresh){
                points(x,as.numeric(obs),pch='*',col="red",cex=2)
                }

        }
abline(h=0.0,lty=2,lwd=2)
legend("topleft",c("All","Pre-2006","Post-2006"),pch=NA,fill=c("black","Lightblue","Darkblue"),pt.cex=2,bty="n",cex=2)
dev.off()