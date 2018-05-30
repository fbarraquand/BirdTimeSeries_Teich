# Script by CPicoche 2018
#The aim is to compare synchrony indices at different taxonomic levels
rm(list=ls())
graphics.off()
library("R.utils")
source("SCRIPTS/test_synchrony_Gross.r")

thresh=0.05
pdf("anascal_waderduck.pdf",width=11,height=7)

par(mfrow=c(1,2),mar=c(3,2,1,.5))

anas=loadToEnv("OUT/data_Anas.RData")
calidris=loadToEnv("OUT/data_Calidris.RData")

nn=names(anas)

plou_taxo=list()

for (nid in 1:length(nn)){
	n=nn[nid]
	anas_tmp=anas[[n]]
	anas_all=aggregate(anas_tmp[,3],list(anas_tmp$dates),sum)
	anas_all=cbind(rep("Anas",length(anas_all[[1]])),anas_all)
	names(anas_all)=c("sp_data_frame","dates","abundance")
	calidris_tmp=calidris[[n]]
	calidris_all=aggregate(calidris_tmp[,3],list(calidris_tmp$dates),sum)
	calidris_all=cbind(rep("Calidris",length(calidris_all[[1]])),calidris_all)
	names(calidris_all)=c("sp_data_frame","dates","abundance")
	if(sum(calidris_all[2]!=anas_all[2])>0){
		print(paste('Probleme in dates for',n,sep=" "))
	}else{
		new_data_frame=merge(anas_all,calidris_all,all=T)
		plou_taxo[[nid]]=community_sync_Gross(new_data_frame,nrands=100)
	}
}
names(plou_taxo)=nn

#For now, we want wintering and breeding
val=c("wintering_all","wintering_pre_2006","wintering_post_2006","breeding_all","breeding_pre_2006","breeding_post_2006")
color=rep(c("Black","Lightblue","Darkblue"),2)
plot(0,0,t="n",ylim=c(0.0,1.1),xlim=c(0,7.5),xaxt="n",xlab="",ylab="")
axis(1,at=c(2,6),lab=c("Cold","Warm"),cex.axis=2)
for (v in 1:length(val)){
	obs=plou_taxo[[val[v]]]$obs
	ic=2*sd(plou_taxo[[val[v]]]$rands[1:100])/sqrt(100)
	p_s=plou_taxo[[val[v]]]$pval
	x=v+(v>3)
        yHigh=as.numeric(obs)+as.numeric(ic)
        yLow=as.numeric(obs)-as.numeric(ic)
	points(x,obs,pch=21,col=color[v],bg=color[v],cex=2)
        arrows(x,yHigh,x,yLow,angle=90,length=0.1,code=3,col=color[v],lwd=2)
	if (p_s<thresh){
        	points(x,as.numeric(obs),pch='*',col="red",cex=2)
        	}

	}
legend("topleft",c("All","Pre-2006","Post-2006"),pch=NA,fill=c("black","Lightblue","Darkblue"),pt.cex=2,bty="n",cex=2)

duck=loadToEnv("OUT/data_freq.RData")
wader=loadToEnv("OUT/data_limicoles.RData")

plou_func=list()

for (nid in 1:length(val)){
        n=val[nid]
        duck_tmp=duck[[n]]
        duck_all=aggregate(duck_tmp[,3],list(duck_tmp$dates),sum)
        duck_all=cbind(rep("Ducks",length(duck_all[[1]])),duck_all)
        names(duck_all)=c("sp_data_frame","dates","abundance")
        wader_tmp=wader[[n]]
        wader_all=aggregate(wader_tmp[,3],list(wader_tmp$dates),sum)
        wader_all=cbind(rep("Waders",length(wader_all[[1]])),wader_all)
        names(wader_all)=c("sp_data_frame","dates","abundance")
        if(sum(wader_all[2]!=duck_all[2])>0){
                print(paste('Probleme in dates for',n,sep=" "))
        }else{  
                new_data_frame=merge(wader_all,duck_all,all=T)
                plou_func[[nid]]=community_sync_Gross(new_data_frame,nrands=100)
        }
}
names(plou_func)=val
plot(0,0,t="n",ylim=c(0.0,1.1),xlim=c(0,7.5),xaxt="n",xlab="",ylab="",yaxt="n")
axis(1,at=c(2,6),lab=c("Cold","Warm"),cex.axis=2)
for (v in 1:length(val)){
        obs=plou_func[[val[v]]]$obs
        ic=2*sd(plou_func[[val[v]]]$rands[1:100])/sqrt(100)
        p_s=plou_func[[val[v]]]$pval
        x=v+(v>3)
        yHigh=as.numeric(obs)+as.numeric(ic)
        yLow=as.numeric(obs)-as.numeric(ic)
        points(x,obs,pch=21,col=color[v],bg=color[v],cex=2)
        arrows(x,yHigh,x,yLow,angle=90,length=0.1,code=3,col=color[v],lwd=2)
        if (p_s<thresh){
                points(x,as.numeric(obs),pch='*',col="red",cex=2)
                }

        }
dev.off()
