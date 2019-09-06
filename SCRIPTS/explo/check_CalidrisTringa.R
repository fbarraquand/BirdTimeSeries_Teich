### CP 2019/09/04 Addendum for Calidris Tringa

rm(list=ls())
graphics.off()

box_index=T #If box_index is true, draw boxplots for the whole distribution for the Gross index with shifted time series. Otherwise, we just show a line for the 5%-95% percentiles
type_correct="BH" #Was bonferroni before
amethod_b="iaaft" #Could also be shift or ebisuzaki. Was shift before.
amethod_w="shift"
anrands=100 #Number of surrogates we want to test the significance. Up to now (2019/07/03), 100.
biomass=F

source("SCRIPTS/test_synchrony_Gross.r")
set.seed(42)
thresh=0.1

#Anas and Calidris
db_warm=read.csv("IN/warmseason_abundances_asdataframe_summed_v2_wtoutrarespecies.csv",sep=";",header=T)
db_cold=read.csv("IN/coldseason_abundances_asdataframe_summed_v2_wtoutrarespecies.csv",sep=";",header=T)

#Right format for synchrony scripts
names(db_warm)=c("dates","sp_data_frame","abundance")
names(db_cold)=c("dates","sp_data_frame","abundance")

db_warm_all=subset(db_warm,sp_data_frame %in% c("Tringa","Calidris")&dates<2016)
db_warm_pre_2006=subset(db_warm,sp_data_frame %in% c("Tringa","Calidris")&dates<=2006)
db_warm_post_2006=subset(db_warm,sp_data_frame %in% c("Tringa","Calidris")&dates>2006&dates<2016)
db_cold_all=subset(db_cold,sp_data_frame %in% c("Tringa","Calidris")&dates<2016)
db_cold_pre_2006=subset(db_cold,sp_data_frame %in% c("Tringa","Calidris")&dates<=2006)
db_cold_post_2006=subset(db_cold,sp_data_frame %in% c("Tringa","Calidris")&dates>2006&dates<2016)

#Compute synchrony values
synch_warm_all=community_sync_Gross(db_warm_all,nrands=anrands,method=amethod_b)
synch_warm_pre_2006=community_sync_Gross(db_warm_pre_2006,nrands=anrands,method=amethod_b)
synch_warm_post_2006=community_sync_Gross(db_warm_post_2006,nrands=anrands,method=amethod_b)

synch_cold_all=community_sync_Gross(db_cold_all,nrands=anrands,method=amethod_b)
synch_cold_pre_2006=community_sync_Gross(db_cold_pre_2006,nrands=anrands,method=amethod_b)
synch_cold_post_2006=community_sync_Gross(db_cold_post_2006,nrands=anrands,method=amethod_b)

#Plot everything
upsi=0.05
essai_taxo=list(synch_cold_all,synch_cold_pre_2006,synch_cold_post_2006,synch_warm_all,synch_warm_pre_2006,synch_warm_post_2006)

mat=rep(NA,length(essai_taxo))
for(v in 1:length(essai_taxo)){
  mat[v]=essai_taxo[[v]]$pval
}
mat_adj=p.adjust(mat,method=type_correct)
for(v in 1:length(essai_taxo)){
  essai_taxo[[v]]$pval=mat_adj[v]
}

pdf("OUT/Synchrony_Calidris_Tringa.pdf")
color=rep(c("Black","Lightblue","Darkblue"),2)
#par(mfrow=c(1,2),mar=c(3,4.5,1,.25),oma=c(1,2,1,.25),mgp=c(3,1,0),xpd=NA)
plot(0,0,t="n",ylim=c(-1,1),xlim=c(0,7.5),xaxt="n",xlab="",ylab="Synchrony index",cex.lab=2,cex.axis=2,las=1)
axis(1,at=c(2,6),lab=c("Cold","Warm"),cex.axis=2)
mtext("Between",side=2,line=0.3,outer=T,cex=2,font=2,at=0.25)
mtext("c)",side=2,line=-4,at=0.49,cex=2,outer=T,las=1)

for (v in 1:6){
  obs=essai_taxo[[v]]$obs
  #        ic=2*sd(essai_taxo[[v]]$rands[1:100])/sqrt(100)
  p_s=essai_taxo[[v]]$pval
  x=v+(v>3)
  #        yHigh=as.numeric(obs)+as.numeric(ic)
  #        yLow=as.numeric(obs)-as.numeric(ic)
  points(x,obs,pch=21,col=color[v],bg=color[v],cex=2)
  #        arrows(x,yHigh,x,yLow,angle=90,length=0.1,code=3,col=color[v],lwd=2)
  if(box_index){
    boxplot(essai_taxo[[v]]$rands[1:100],at=x,add=T,boxwex=0.25,range=0,yaxt="n",xaxt="n")
  }else{
    plou=quantile(essai_taxo[[v]]$rands[1:100],c(0.05,0.95))
    arrows(x,plou[1],x,plou[2],angle=90,length=0.1,code=3,col=color[v],lwd=2)
  }
  if (p_s<=thresh){
    points(x,as.numeric(obs)+upsi,pch='*',col="red",cex=2)
  }
  
}
lines(c(0.0,7.5),c(0,0),lty=2,lwd=2)
ll=c(essai_taxo[[2]]$obs,essai_taxo[[3]]$obs)
#lines(2:3,ll,col="black",lwd=2,lty=1)
#ll=c(essai_taxo[[5]]$obs,essai_taxo[[6]]$obs)
#lines(6:7,ll,col="black",lwd=2,lty=1)
legend("bottomleft",c("All","Pre-2006","Post-2006"),pch=NA,fill=c("black","Lightblue","Darkblue"),pt.cex=2,bty="n",cex=2)
legend("topleft",c("Tringa/Calidris"),pch=21,pt.bg=c("black"),pt.cex=2,bty="n",cex=2)
dev.off()
