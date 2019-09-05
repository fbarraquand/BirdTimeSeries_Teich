# Script by CPicoche 2018
#The aim is to compare synchrony indices at different taxonomic levels
#TO DO : rewrite the code to use the seasonal data frame
rm(list=ls())
graphics.off()
source("SCRIPTS/test_synchrony_Gross.r")

thresh=0.05
#Anas and Calidris
db_warm=read.csv("IN/warmseason_abundances_asdataframe_summed.csv",sep=";",header=T)
db_cold=read.csv("IN/coldseason_abundances_asdataframe_summed.csv",sep=";",header=T)

#Right format for synchrony scripts
names(db_warm)=c("dates","sp_data_frame","abundance")
names(db_cold)=c("dates","sp_data_frame","abundance")

db_warm_all=subset(db_warm,sp_data_frame %in% c("Anas","Calidris")&dates<2016)
db_warm_pre_2006=subset(db_warm,sp_data_frame %in% c("Anas","Calidris")&dates<=2006)
db_warm_post_2006=subset(db_warm,sp_data_frame %in% c("Anas","Calidris")&dates>2006&dates<2016)
db_cold_all=subset(db_cold,sp_data_frame %in% c("Anas","Calidris")&dates<2016)
db_cold_pre_2006=subset(db_cold,sp_data_frame %in% c("Anas","Calidris")&dates<=2006)
db_cold_post_2006=subset(db_cold,sp_data_frame %in% c("Anas","Calidris")&dates>2006&dates<2016)

#Compute synchrony values
synch_warm_all=community_sync_Gross(db_warm_all,nrands=100)
synch_warm_pre_2006=community_sync_Gross(db_warm_pre_2006,nrands=100)
synch_warm_post_2006=community_sync_Gross(db_warm_post_2006,nrands=100)

synch_cold_all=community_sync_Gross(db_cold_all,nrands=100)
synch_cold_pre_2006=community_sync_Gross(db_cold_pre_2006,nrands=100)
synch_cold_post_2006=community_sync_Gross(db_cold_post_2006,nrands=100)

#Plot everything
upsi=0.05
essai_taxo=list(synch_cold_all,synch_cold_pre_2006,synch_cold_post_2006,synch_warm_all,synch_warm_pre_2006,synch_warm_post_2006)
color=rep(c("Black","Lightblue","Darkblue"),2)
pdf("OUT/gross_anascalidris_wader_duck.pdf",width=11,height=7,family="Times")
par(mfrow=c(1,2),mar=c(3,4.5,1,.25),oma=c(1,2,1,.25),mgp=c(3,1,0),xpd=NA)
plot(0,0,t="n",ylim=c(-.75,.66),xlim=c(0,7.5),xaxt="n",xlab="",ylab="Synchrony index",cex.lab=2,cex.axis=2)
axis(1,at=c(2,6),lab=c("Cold","Warm"),cex.axis=2)
mtext("Between",side=2,line=0.5,outer=T,cex=2)
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
                points(x,as.numeric(obs)+upsi,pch='*',col="red",cex=2)
                }

        }
lines(c(0.0,7.5),c(0,0),lty=2,lwd=2)
ll=c(essai_taxo[[2]]$obs,essai_taxo[[3]]$obs)
lines(2:3,ll,col="black",lwd=2,lty=1)
ll=c(essai_taxo[[5]]$obs,essai_taxo[[6]]$obs)
lines(6:7,ll,col="black",lwd=2,lty=1)

legend("bottomleft",c("All","Pre-2006","Post-2006"),pch=NA,fill=c("black","Lightblue","Darkblue"),pt.cex=2,bty="n",cex=2)
legend("topleft",c("Anas/Calidris"),pch=21,pt.bg=c("black"),pt.cex=2,bty="n",cex=2)

db_warm_all=subset(db_warm,sp_data_frame %in% c("Waders","Ducks")&dates<2016)
db_warm_pre_2006=subset(db_warm,sp_data_frame %in% c("Waders","Ducks")&dates<=2006)
db_warm_post_2006=subset(db_warm,sp_data_frame %in% c("Waders","Ducks")&dates>2006&dates<2016)
db_cold_all=subset(db_cold,sp_data_frame %in% c("Waders","Ducks")&dates<2016)
db_cold_pre_2006=subset(db_cold,sp_data_frame %in% c("Waders","Ducks")&dates<=2006)
db_cold_post_2006=subset(db_cold,sp_data_frame %in% c("Waders","Ducks")&dates>2006&dates<2016) #Because there is a problem for 2016 (na values)

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
plot(0,0,t="n",ylim=c(-.75,.66),xlim=c(0,7.5),xaxt="n",xlab="",ylab="",yaxt="n")
axis(1,at=c(2,6),lab=c("Cold","Warm"),cex.axis=2)
for (v in 1:6){
        obs=essai_taxo[[v]]$obs
        ic=2*sd(essai_taxo[[v]]$rands[1:100])/sqrt(100)
        p_s=essai_taxo[[v]]$pval
        x=v+(v>3)
        yHigh=as.numeric(obs)+as.numeric(ic)
        yLow=as.numeric(obs)-as.numeric(ic)
        points(x,obs,pch=22,col=color[v],bg=color[v],cex=2)
        arrows(x,yHigh,x,yLow,angle=90,length=0.1,code=3,col=color[v],lwd=2)
        if (p_s<thresh){
                points(x,as.numeric(obs)+upsi,pch='*',col="red",cex=2)
                }

        }
lines(c(0.0,7.5),c(0,0),lty=2,lwd=2)
ll=c(essai_taxo[[2]]$obs,essai_taxo[[3]]$obs)
lines(2:3,ll,col="black",lwd=2,lty=1)
ll=c(essai_taxo[[5]]$obs,essai_taxo[[6]]$obs)
lines(6:7,ll,col="black",lwd=2,lty=1)

legend("topleft",c("Waders/Ducks"),pch=c(22),pt.bg=c("black"),pt.cex=2,bty="n",cex=2)
dev.off()

