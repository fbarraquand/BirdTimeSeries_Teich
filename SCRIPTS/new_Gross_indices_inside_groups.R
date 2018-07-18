rm(list=ls())
graphics.off()

source("SCRIPTS/test_synchrony_Gross.r")

thresh=0.05

sp_to_ignore=c("Anas discors","Anas americana","Calidris melanotos","Calidris pusilla","Calidris ruficollis", "Calidris fuscicollis", "Calidris himantopus", "Burhinus oedicnemus","Phalaropus lobatus","Charadrius alexandrinus","Haematopus ostralegus","Calidris maritima","Aythya nyroca","Bucephala clangula","Melanitta nigra","Mergus serrator","Clangula hyemalis","Alopochen aegyptiaca", "Aix galericulata","Cygnus atratus","Tadorna ferruginea","Branta leucopsis","Anser fabalis","Anser albifrons","Cygnus cygnus","Mergus merganser","Anser barchyrhynchus")
#I think there are more species we should ignore, have a closer look with CHristelle

print("Anas")
db_warm=read.csv("IN/warmseason_anas_detailed.txt",sep=",",header=T)
db_cold=read.csv("IN/coldseason_anas_detailed.txt",sep=",",header=T)
db_warm=db_warm[,c(2,3,4)]
db_cold=db_cold[,c(2,3,4)]
#Right format for synchrony scripts
names(db_warm)=c("dates","sp_data_frame","abundance")
names(db_cold)=c("dates","sp_data_frame","abundance")

db_warm_all=subset(db_warm,!(sp_data_frame %in%sp_to_ignore) &dates<2016)
db_warm_pre_2006=subset(db_warm,!(sp_data_frame %in%sp_to_ignore)  &dates<=2006)
db_warm_post_2006=subset(db_warm,!(sp_data_frame %in%sp_to_ignore)  &dates>2006&dates<2016)
db_cold_all=subset(db_cold,!(sp_data_frame %in%sp_to_ignore)  &dates<2016)
db_cold_pre_2006=subset(db_cold,!(sp_data_frame %in%sp_to_ignore) &dates<=2006)
db_cold_post_2006=subset(db_cold,!(sp_data_frame %in%sp_to_ignore) &dates>2006&dates<2016)

#Compute synchrony values
synch_warm_all=community_sync_Gross(db_warm_all,nrands=100)
synch_warm_pre_2006=community_sync_Gross(db_warm_pre_2006,nrands=100)
synch_warm_post_2006=community_sync_Gross(db_warm_post_2006,nrands=100)

synch_cold_all=community_sync_Gross(db_cold_all,nrands=100)
synch_cold_pre_2006=community_sync_Gross(db_cold_pre_2006,nrands=100)
synch_cold_post_2006=community_sync_Gross(db_cold_post_2006,nrands=100)

list_anas=list(synch_cold_all,synch_cold_pre_2006,synch_cold_post_2006,synch_warm_all,synch_warm_pre_2006,synch_warm_post_2006)

print("Calidris")
#We look at only 6 species for the Calidris, we could also look at only 4
db_warm=read.csv("IN/warmseason_calidris_detailed.txt",sep=",",header=T)
db_cold=read.csv("IN/coldseason_calidris_detailed.txt",sep=",",header=T)

db_warm=db_warm[,c(2,3,4)]
db_cold=db_cold[,c(2,3,4)]
#Right format for synchrony scripts
names(db_warm)=c("dates","sp_data_frame","abundance")
names(db_cold)=c("dates","sp_data_frame","abundance")

db_warm_all=subset(db_warm,!(sp_data_frame %in%sp_to_ignore) &dates<2016)
db_warm_pre_2006=subset(db_warm,!(sp_data_frame %in%sp_to_ignore)  &dates<=2006)
db_warm_post_2006=subset(db_warm,!(sp_data_frame %in%sp_to_ignore)  &dates>2006&dates<2016)
db_cold_all=subset(db_cold,!(sp_data_frame %in%sp_to_ignore)  &dates<2016)
db_cold_pre_2006=subset(db_cold,!(sp_data_frame %in%sp_to_ignore) &dates<=2006)
db_cold_post_2006=subset(db_cold,!(sp_data_frame %in%sp_to_ignore) &dates>2006&dates<2016)

#Compute synchrony values
synch_warm_all=community_sync_Gross(db_warm_all,nrands=100)
synch_warm_pre_2006=community_sync_Gross(db_warm_pre_2006,nrands=100)
synch_warm_post_2006=community_sync_Gross(db_warm_post_2006,nrands=100)

synch_cold_all=community_sync_Gross(db_cold_all,nrands=100)
synch_cold_pre_2006=community_sync_Gross(db_cold_pre_2006,nrands=100)
synch_cold_post_2006=community_sync_Gross(db_cold_post_2006,nrands=100)

list_calidris=list(synch_cold_all,synch_cold_pre_2006,synch_cold_post_2006,synch_warm_all,synch_warm_pre_2006,synch_warm_post_2006)

print("Waders")
db_warm=read.csv("IN/warmseason_waders_detailed.txt",sep=",",header=T)
db_cold=read.csv("IN/coldseason_waders_detailed.txt",sep=",",header=T)


db_warm=db_warm[,c(2,3,4)]
db_cold=db_cold[,c(2,3,4)]
#Right format for synchrony scripts
names(db_warm)=c("dates","sp_data_frame","abundance")
names(db_cold)=c("dates","sp_data_frame","abundance")
limicoles=unique(db_warm$sp_data_frame)

db_warm_all=subset(db_warm,!(sp_data_frame %in%sp_to_ignore) &dates<2016)
db_warm_pre_2006=subset(db_warm,!(sp_data_frame %in%sp_to_ignore)  &dates<=2006)
db_warm_post_2006=subset(db_warm,!(sp_data_frame %in%sp_to_ignore)  &dates>2006&dates<2016)
db_cold_all=subset(db_cold,!(sp_data_frame %in%sp_to_ignore)  &dates<2016)
db_cold_pre_2006=subset(db_cold,!(sp_data_frame %in%sp_to_ignore) &dates<=2006)
db_cold_post_2006=subset(db_cold,!(sp_data_frame %in%sp_to_ignore) &dates>2006&dates<2016)

#Compute synchrony values
synch_warm_all=community_sync_Gross(db_warm_all,nrands=100)
synch_warm_pre_2006=community_sync_Gross(db_warm_pre_2006,nrands=100)
synch_warm_post_2006=community_sync_Gross(db_warm_post_2006,nrands=100)

synch_cold_all=community_sync_Gross(db_cold_all,nrands=100)
synch_cold_pre_2006=community_sync_Gross(db_cold_pre_2006,nrands=100)
synch_cold_post_2006=community_sync_Gross(db_cold_post_2006,nrands=100)

list_waders=list(synch_cold_all,synch_cold_pre_2006,synch_cold_post_2006,synch_warm_all,synch_warm_pre_2006,synch_warm_post_2006)

print("Freq")
db_warm=read.csv("IN/warmseason_freq_detailed.txt",sep=",",header=T)
db_cold=read.csv("IN/coldseason_freq_detailed.txt",sep=",",header=T)
db_warm=db_warm[,c(2,3,4)]
db_cold=db_cold[,c(2,3,4)]
#Right format for synchrony scripts
names(db_warm)=c("dates","sp_data_frame","abundance")
names(db_cold)=c("dates","sp_data_frame","abundance")

#db_warm_all=subset(db_warm,!(sp_data_frame %in%sp_to_ignore)&!(sp_data_frame %in%limicoles)  &dates<2016)
#db_warm_pre_2006=subset(db_warm,!(sp_data_frame %in%sp_to_ignore)&!(sp_data_frame %in%limicoles)   &dates<=2006)
#db_warm_post_2006=subset(db_warm,!(sp_data_frame %in%sp_to_ignore)&!(sp_data_frame %in%limicoles)  &dates>2006&dates<2016)
#db_cold_all=subset(db_cold,!(sp_data_frame %in%sp_to_ignore)&!(sp_data_frame %in%limicoles)  &dates<2016)
#db_cold_pre_2006=subset(db_cold,!(sp_data_frame %in%sp_to_ignore)&!(sp_data_frame %in%limicoles) &dates<=2006)
#db_cold_post_2006=subset(db_cold,!(sp_data_frame %in%sp_to_ignore)&!(sp_data_frame %in%limicoles) &dates>2006&dates<2016)

db_warm_all=subset(db_warm,!(sp_data_frame %in%sp_to_ignore)  &dates<2016)
db_warm_pre_2006=subset(db_warm,!(sp_data_frame %in%sp_to_ignore)   &dates<=2006)
db_warm_post_2006=subset(db_warm,!(sp_data_frame %in%sp_to_ignore)  &dates>2006&dates<2016)
db_cold_all=subset(db_cold,!(sp_data_frame %in%sp_to_ignore)  &dates<2016)
db_cold_pre_2006=subset(db_cold,!(sp_data_frame %in%sp_to_ignore) &dates<=2006)
db_cold_post_2006=subset(db_cold,!(sp_data_frame %in%sp_to_ignore) &dates>2006&dates<2016)


#Compute synchrony values
synch_warm_all=community_sync_Gross(db_warm_all,nrands=100)
synch_warm_pre_2006=community_sync_Gross(db_warm_pre_2006,nrands=100)
synch_warm_post_2006=community_sync_Gross(db_warm_post_2006,nrands=100)

synch_cold_all=community_sync_Gross(db_cold_all,nrands=100)
synch_cold_pre_2006=community_sync_Gross(db_cold_pre_2006,nrands=100)
synch_cold_post_2006=community_sync_Gross(db_cold_post_2006,nrands=100)

list_freq=list(synch_cold_all,synch_cold_pre_2006,synch_cold_post_2006,synch_warm_all,synch_warm_pre_2006,synch_warm_post_2006)

print("Ducks")
db_warm=read.csv("IN/warmseason_duck_detailed.txt",sep=",",header=T)
db_cold=read.csv("IN/coldseason_duck_detailed.txt",sep=",",header=T)
db_warm=db_warm[,c(2,3,4)]
db_cold=db_cold[,c(2,3,4)]
#Right format for synchrony scripts
names(db_warm)=c("dates","sp_data_frame","abundance")
names(db_cold)=c("dates","sp_data_frame","abundance")

db_warm_all=subset(db_warm,!(sp_data_frame %in%sp_to_ignore)  &dates<2016)
db_warm_pre_2006=subset(db_warm,!(sp_data_frame %in%sp_to_ignore)   &dates<=2006)
db_warm_post_2006=subset(db_warm,!(sp_data_frame %in%sp_to_ignore)  &dates>2006&dates<2016)
db_cold_all=subset(db_cold,!(sp_data_frame %in%sp_to_ignore)  &dates<2016)
db_cold_pre_2006=subset(db_cold,!(sp_data_frame %in%sp_to_ignore) &dates<=2006)
db_cold_post_2006=subset(db_cold,!(sp_data_frame %in%sp_to_ignore)&dates>2006&dates<2016)

#Compute synchrony values
synch_warm_all=community_sync_Gross(db_warm_all,nrands=100)
synch_warm_pre_2006=community_sync_Gross(db_warm_pre_2006,nrands=100)
synch_warm_post_2006=community_sync_Gross(db_warm_post_2006,nrands=100)

synch_cold_all=community_sync_Gross(db_cold_all,nrands=100)
synch_cold_pre_2006=community_sync_Gross(db_cold_pre_2006,nrands=100)
synch_cold_post_2006=community_sync_Gross(db_cold_post_2006,nrands=100)

list_duck=list(synch_cold_all,synch_cold_pre_2006,synch_cold_post_2006,synch_warm_all,synch_warm_pre_2006,synch_warm_post_2006)


#PLOT
upsi=0.05
color=rep(c("Black","Lightblue","Darkblue"),2)
if(1==0){
pdf("OUT/gross_inside_groups.pdf",width=11,height=7,family="Times")
par(mfrow=c(1,2),mar=c(3,4.5,2,.25),oma=c(1,2,1,.25),mgp=c(3,1,0),xpd=NA)
plot(0,0,t="n",ylim=c(-.75,.66),xlim=c(0,7.5),xaxt="n",xlab="",ylab="Synchrony index",cex.axis=2,cex.lab=2,main="Taxonomic groups",cex.main=2)
mtext("Within",side=2,line=0.5,outer=T,cex=2)
axis(1,at=c(2,6),lab=c("Cold","Warm"),cex.axis=2)
for (v in 1:6){
	essai_taxo=list_anas
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

        essai_taxo=list_calidris
        obs=essai_taxo[[v]]$obs
        ic=2*sd(essai_taxo[[v]]$rands[1:100])/sqrt(100)
        p_s=essai_taxo[[v]]$pval
        x=v+(v>3)
        yHigh=as.numeric(obs)+as.numeric(ic)
        yLow=as.numeric(obs)-as.numeric(ic)
        points(x,obs,pch=22,col=color[v],bg=color[v],cex=2)
        arrows(x,yHigh,x,yLow,angle=90,length=0.1,code=3,col=color[v],lwd=2,lty=2)
        if (p_s<thresh){
                points(x,as.numeric(obs)+upsi,pch='*',col="red",cex=2)
                }

        }
ll=c(list_anas[[2]]$obs,list_anas[[3]]$obs)
lines(2:3,ll,col="black",lwd=2,lty=1)
ll=c(list_anas[[5]]$obs,list_anas[[6]]$obs)
lines(6:7,ll,col="black",lwd=2,lty=1)

ll=c(list_calidris[[2]]$obs,list_calidris[[3]]$obs)
lines(2:3,ll,col="black",lwd=2,lty=2)
ll=c(list_calidris[[5]]$obs,list_calidris[[6]]$obs)
lines(6:7,ll,col="black",lwd=2,lty=2)


lines(c(0,7.5),c(0.,0.),lty=2,lwd=2)
legend("bottomleft",c("All","Pre-2006","Post-2006"),pch=NA,fill=c("black","Lightblue","Darkblue"),pt.cex=2,bty="n",cex=2)
legend("bottomright",c("Anas","Calidris"),pch=c(21,22),pt.bg=c("black"),pt.cex=2,bty="n",cex=2)

plot(0,0,t="n",ylim=c(-.75,.66),xlim=c(0,7.5),xaxt="n",xlab="",ylab="",yaxt="n",main='Functional groups',cex.main=2)
axis(1,at=c(2,6),lab=c("Cold","Warm"),cex.axis=2)
for (v in 1:6){
        essai_taxo=list_waders
        obs=essai_taxo[[v]]$obs
        ic=2*sd(essai_taxo[[v]]$rands[1:100])/sqrt(100)
        p_s=essai_taxo[[v]]$pval
        x=v+(v>3)
        yHigh=as.numeric(obs)+as.numeric(ic)
        yLow=as.numeric(obs)-as.numeric(ic)
        points(x,obs,pch=23,col=color[v],bg=color[v],cex=2)
        arrows(x,yHigh,x,yLow,angle=90,length=0.1,code=3,col=color[v],lwd=2)
        if (p_s<thresh){
                points(x,as.numeric(obs)+upsi,pch='*',col="red",cex=2)
                }

        essai_taxo=list_duck
        obs=essai_taxo[[v]]$obs
        ic=2*sd(essai_taxo[[v]]$rands[1:100])/sqrt(100)
        p_s=essai_taxo[[v]]$pval
        x=v+(v>3)
        yHigh=as.numeric(obs)+as.numeric(ic)
        yLow=as.numeric(obs)-as.numeric(ic)
        points(x,obs,pch=24,col=color[v],bg=color[v],cex=2)
        arrows(x,yHigh,x,yLow,angle=90,length=0.1,code=3,col=color[v],lwd=2,lt=2)
        if (p_s<thresh){
                points(x,as.numeric(obs)+upsi,pch='*',col="red",cex=2)
                }

        }

ll=c(list_waders[[2]]$obs,list_waders[[3]]$obs)
lines(2:3,ll,col="black",lwd=2,lty=1)
ll=c(list_waders[[5]]$obs,list_waders[[6]]$obs)
lines(6:7,ll,col="black",lwd=2,lty=1)

ll=c(list_duck[[2]]$obs,list_duck[[3]]$obs)
lines(2:3,ll,col="black",lwd=2,lty=2)
ll=c(list_duck[[5]]$obs,list_duck[[6]]$obs)
lines(6:7,ll,col="black",lwd=2,lty=2)

lines(c(0,7.5),c(0.,0.),lty=2,lwd=2)
legend("bottomright",c("Waders","Ducks"),pch=c(23,24),pt.bg="black",pt.cex=2,bty="n",cex=2)
dev.off()
}

upsi=0.02
pdf("OUT/Gross_freq_including_waders.pdf",width=11,height=7,family="Times")
par(mfrow=c(1,1),mar=c(3,4.5,2,.25),oma=c(1,2,1,.25),mgp=c(3,1,0),xpd=NA)
plot(0,0,t="n",ylim=c(-0.125,0.425),xlim=c(0,7.5),xaxt="n",xlab="",ylab="Synchrony index",cex.axis=2,cex.lab=2,main="",cex.main=2)
mtext("Within",side=2,line=0.5,outer=T,cex=2)
axis(1,at=c(2,6),lab=c("Cold","Warm"),cex.axis=2)
for (v in 1:6){
        essai_taxo=list_freq
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
lines(c(0,7.5),c(0.,0.),lty=2,lwd=2)
legend("topleft",c("All","Pre-2006","Post-2006"),pch=NA,fill=c("black","Lightblue","Darkblue"),pt.cex=2,bty="n",cex=2)
legend("bottomright",c("Freq"),pch=c(21),pt.bg="black",pt.cex=2,bty="n",cex=2)
dev.off()
