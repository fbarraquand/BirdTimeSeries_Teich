####CP: original script was written in 2018 to compute the Gross index for synchrony (Gross et al. 2014) on taxonomic and functional groups of birds, and to assess their significance. 
#### 2019/06/25: Corrected an error on standard deviation. They were computed under H0 (on shifted time series) and were associated to the observed value (errobars were built around the observed value with the sd from rands)
### 2019/06/27: Added BH p-val correction. Tried IAAFT and Ebisuzaki surrogates for the "between" case, and proposed to also use them for the "within" case.
### 2019/07/03: Corrected the way rare species were ignored (not taken into account the right way before). This does not change anything, actually, as the Gross function already removed them, but this is cleaner (and removes the warnings).
### 2019/09/04: Changed the script to adapt to biomasses
### 2019/09/05: Presentation details (las=1, Anas-> Anatini, etc.)

rm(list=ls())
graphics.off()

box_index=T #If box_index is true, draw boxplots for the whole distribution for the Gross index with shifted time series. Otherwise, we just show a line for the 5%-95% percentiles
type_correct="BH" #Was bonferroni before
amethod_b="iaaft" #Could also be shift or ebisuzaki. Was shift before.
amethod_w="shift"
anrands=100 #Number of surrogates we want to test the significance. Up to now (2019/07/03), 100.
biomass=F
if(biomass){
end_bio="biomasses"
}else{
end_bio="abundances"
}

source("SCRIPTS/test_synchrony_Gross.r")
set.seed(42)
thresh=0.1

sp_to_ignore=c("Anas discors","Anas americana","Calidris melanotos","Calidris pusilla","Calidris ruficollis", "Calidris fuscicollis", "Calidris himantopus", "Burhinus oedicnemus","Phalaropus lobatus","Charadrius alexandrinus","Haematopus ostralegus","Calidris maritima","Aythya nyroca","Bucephala clangula","Melanitta nigra","Mergus serrator","Clangula hyemalis","Alopochen aegyptiaca", "Aix galericulata","Cygnus atratus","Tadorna ferruginea","Branta leucopsis","Anser fabalis","Anser albifrons","Cygnus cygnus","Mergus merganser","Anser brachyrhynchus")

print("Anas")
db_warm=read.csv(paste("IN/warmseason_anas_detailed_v2_",end_bio,"_wtoutrarespecies.txt",sep=""))
db_cold=read.csv(paste("IN/coldseason_anas_detailed_v2_",end_bio,"_wtoutrarespecies.txt",sep=""))
db_warm=db_warm[,c(2,3,4)]
db_cold=db_cold[,c(2,3,4)]
#Right format for synchrony scripts
names(db_warm)=c("dates","sp_data_frame","abundance")
names(db_cold)=c("dates","sp_data_frame","abundance")

db_warm_all=subset(db_warm,!(as.character(sp_data_frame) %in%sp_to_ignore) &dates<2016)
db_warm_pre_2006=subset(db_warm,!(as.character(sp_data_frame) %in%sp_to_ignore)  &dates<=2006)
db_warm_post_2006=subset(db_warm,!(as.character(sp_data_frame) %in%sp_to_ignore)  &dates>2006&dates<2016)
db_cold_all=subset(db_cold,!(as.character(sp_data_frame) %in%sp_to_ignore)  &dates<2016)
db_cold_pre_2006=subset(db_cold,!(as.character(sp_data_frame) %in%sp_to_ignore) &dates<=2006)
db_cold_post_2006=subset(db_cold,!(as.character(sp_data_frame) %in%sp_to_ignore) &dates>2006&dates<2016)

#Compute synchrony values
synch_warm_all=community_sync_Gross(db_warm_all,nrands=anrands,method=amethod_w)
synch_warm_pre_2006=community_sync_Gross(db_warm_pre_2006,nrands=anrands,method=amethod_w)
synch_warm_post_2006=community_sync_Gross(db_warm_post_2006,nrands=anrands,method=amethod_w)

synch_cold_all=community_sync_Gross(db_cold_all,nrands=anrands,method=amethod_w)
synch_cold_pre_2006=community_sync_Gross(db_cold_pre_2006,nrands=anrands,method=amethod_w)
synch_cold_post_2006=community_sync_Gross(db_cold_post_2006,nrands=anrands,method=amethod_w)

list_anas=list(synch_cold_all,synch_cold_pre_2006,synch_cold_post_2006,synch_warm_all,synch_warm_pre_2006,synch_warm_post_2006)

mat=rep(NA,length(list_anas))
for(v in 1:length(list_anas)){
	mat[v]=list_anas[[v]]$pval
}
mat_adj=p.adjust(mat,method=type_correct)
for(v in 1:length(list_anas)){
	list_anas[[v]]$pval=mat_adj[v]
}

print("Calidris")
#We look at only 6 species for the Calidris, we could also look at only 4
db_warm=read.csv(paste("IN/warmseason_calidris_detailed_v2_",end_bio,"_wtoutrarespecies.txt",sep=""))
db_cold=read.csv(paste("IN/coldseason_calidris_detailed_v2_",end_bio,"_wtoutrarespecies.txt",sep=""))

db_warm=db_warm[,c(2,3,4)]
db_cold=db_cold[,c(2,3,4)]
#Right format for synchrony scripts
names(db_warm)=c("dates","sp_data_frame","abundance")
names(db_cold)=c("dates","sp_data_frame","abundance")

db_warm_all=subset(db_warm,!(as.character(sp_data_frame) %in%sp_to_ignore) &dates<2016)
db_warm_pre_2006=subset(db_warm,!(as.character(sp_data_frame) %in%sp_to_ignore)  &dates<=2006)
db_warm_post_2006=subset(db_warm,!(as.character(sp_data_frame) %in%sp_to_ignore)  &dates>2006&dates<2016)
db_cold_all=subset(db_cold,!(as.character(sp_data_frame) %in%sp_to_ignore)  &dates<2016)
db_cold_pre_2006=subset(db_cold,!(as.character(sp_data_frame) %in%sp_to_ignore) &dates<=2006)
db_cold_post_2006=subset(db_cold,!(as.character(sp_data_frame) %in%sp_to_ignore) &dates>2006&dates<2016)
#Compute synchrony values

synch_warm_all=community_sync_Gross(db_warm_all,nrands=anrands,method=amethod_w)
synch_warm_pre_2006=community_sync_Gross(db_warm_pre_2006,nrands=anrands,method=amethod_w)
synch_warm_post_2006=community_sync_Gross(db_warm_post_2006,nrands=anrands,method=amethod_w)

synch_cold_all=community_sync_Gross(db_cold_all,nrands=anrands,method=amethod_w)
synch_cold_pre_2006=community_sync_Gross(db_cold_pre_2006,nrands=anrands,method=amethod_w)
synch_cold_post_2006=community_sync_Gross(db_cold_post_2006,nrands=anrands,method=amethod_w)

list_calidris=list(synch_cold_all,synch_cold_pre_2006,synch_cold_post_2006,synch_warm_all,synch_warm_pre_2006,synch_warm_post_2006)
	
mat=rep(NA,length(list_calidris))
for(v in 1:length(list_calidris)){
        mat[v]=list_calidris[[v]]$pval
}
mat_adj=p.adjust(mat,method=type_correct)
for(v in 1:length(list_calidris)){
        list_calidris[[v]]$pval=mat_adj[v]
}


print("Waders")
db_warm=read.csv(paste("IN/warmseason_waders_detailed_v2_",end_bio,"_wtoutrarespecies.txt",sep=""))
db_cold=read.csv(paste("IN/coldseason_waders_detailed_v2_",end_bio,"_wtoutrarespecies.txt",sep=""))

db_warm=db_warm[,c(2,3,4)]
db_cold=db_cold[,c(2,3,4)]
#Right format for synchrony scripts
names(db_warm)=c("dates","sp_data_frame","abundance")
names(db_cold)=c("dates","sp_data_frame","abundance")
limicoles=unique(db_warm$sp_data_frame)

db_warm_all=subset(db_warm,!(as.character(sp_data_frame) %in%sp_to_ignore) &dates<2016)
db_warm_pre_2006=subset(db_warm,!(as.character(sp_data_frame) %in%sp_to_ignore)  &dates<=2006)
db_warm_post_2006=subset(db_warm,!(as.character(sp_data_frame) %in%sp_to_ignore)  &dates>2006&dates<2016)
db_cold_all=subset(db_cold,!(as.character(sp_data_frame) %in%sp_to_ignore)  &dates<2016)
db_cold_pre_2006=subset(db_cold,!(as.character(sp_data_frame) %in%sp_to_ignore) &dates<=2006)
db_cold_post_2006=subset(db_cold,!(as.character(sp_data_frame) %in%sp_to_ignore) &dates>2006&dates<2016)

#Compute synchrony values
synch_warm_all=community_sync_Gross(db_warm_all,nrands=anrands,method=amethod_w)
synch_warm_pre_2006=community_sync_Gross(db_warm_pre_2006,nrands=anrands,method=amethod_w)
synch_warm_post_2006=community_sync_Gross(db_warm_post_2006,nrands=anrands,method=amethod_w)

synch_cold_all=community_sync_Gross(db_cold_all,nrands=anrands,method=amethod_w)
synch_cold_pre_2006=community_sync_Gross(db_cold_pre_2006,nrands=anrands,method=amethod_w)
synch_cold_post_2006=community_sync_Gross(db_cold_post_2006,nrands=anrands,method=amethod_w)

list_waders=list(synch_cold_all,synch_cold_pre_2006,synch_cold_post_2006,synch_warm_all,synch_warm_pre_2006,synch_warm_post_2006)

mat=rep(NA,length(list_waders))
for(v in 1:length(list_waders)){
        mat[v]=list_waders[[v]]$pval
}
mat_adj=p.adjust(mat,method=type_correct)
for(v in 1:length(list_waders)){
        list_waders[[v]]$pval=mat_adj[v]
}


print("Freq")
db_warm=read.csv(paste("IN/warmseason_freq_detailed_v2_",end_bio,"_wtoutrarespecies.txt",sep=""))
db_cold=read.csv(paste("IN/coldseason_freq_detailed_v2_",end_bio,"_wtoutrarespecies.txt",sep=""))
db_warm=db_warm[,c(2,3,4)]
db_cold=db_cold[,c(2,3,4)]
#Right format for synchrony scripts
names(db_warm)=c("dates","sp_data_frame","abundance")
names(db_cold)=c("dates","sp_data_frame","abundance")

db_warm_all=subset(db_warm,!(as.character(sp_data_frame) %in%sp_to_ignore) &dates<2016)
db_warm_pre_2006=subset(db_warm,!(as.character(sp_data_frame) %in%sp_to_ignore)  &dates<=2006)
db_warm_post_2006=subset(db_warm,!(as.character(sp_data_frame) %in%sp_to_ignore)  &dates>2006&dates<2016)
db_cold_all=subset(db_cold,!(as.character(sp_data_frame) %in%sp_to_ignore) &dates<2016)
db_cold_pre_2006=subset(db_cold,!(as.character(sp_data_frame) %in%sp_to_ignore) &dates<=2006)
db_cold_post_2006=subset(db_cold,!(as.character(sp_data_frame) %in%sp_to_ignore)&dates>2006&dates<2016)
#Compute synchrony values
synch_warm_all=community_sync_Gross(db_warm_all,nrands=anrands,method=amethod_w)
synch_warm_pre_2006=community_sync_Gross(db_warm_pre_2006,nrands=anrands,method=amethod_w)
synch_warm_post_2006=community_sync_Gross(db_warm_post_2006,nrands=anrands,method=amethod_w)

synch_cold_all=community_sync_Gross(db_cold_all,nrands=anrands,method=amethod_w)
synch_cold_pre_2006=community_sync_Gross(db_cold_pre_2006,nrands=anrands,method=amethod_w)
synch_cold_post_2006=community_sync_Gross(db_cold_post_2006,nrands=anrands,method=amethod_w)

list_freq=list(synch_cold_all,synch_cold_pre_2006,synch_cold_post_2006,synch_warm_all,synch_warm_pre_2006,synch_warm_post_2006)

mat=rep(NA,length(list_freq))
for(v in 1:length(list_freq)){
        mat[v]=list_freq[[v]]$pval
}
mat_adj=p.adjust(mat,method=type_correct)
for(v in 1:length(list_freq)){
        list_freq[[v]]$pval=mat_adj[v]
}

print("Ducks")
db_warm=read.csv(paste("IN/warmseason_ducks_detailed_v2_",end_bio,"_wtoutrarespecies.txt",sep=""))
db_cold=read.csv(paste("IN/coldseason_ducks_detailed_v2_",end_bio,"_wtoutrarespecies.txt",sep=""))
db_warm=db_warm[,c(2,3,4)]
db_cold=db_cold[,c(2,3,4)]
#Right format for synchrony scripts
names(db_warm)=c("dates","sp_data_frame","abundance")
names(db_cold)=c("dates","sp_data_frame","abundance")

db_warm_all=subset(db_warm,!(as.character(sp_data_frame) %in%sp_to_ignore)  &dates<2016)
db_warm_pre_2006=subset(db_warm,!(as.character(sp_data_frame) %in%sp_to_ignore)   &dates<=2006)
db_warm_post_2006=subset(db_warm,!(as.character(sp_data_frame) %in%sp_to_ignore)  &dates>2006&dates<2016)
db_cold_all=subset(db_cold,!(as.character(sp_data_frame) %in%sp_to_ignore)  &dates<2016)
db_cold_pre_2006=subset(db_cold,!(as.character(sp_data_frame) %in%sp_to_ignore) &dates<=2006)
db_cold_post_2006=subset(db_cold,!(as.character(sp_data_frame) %in%sp_to_ignore)&dates>2006&dates<2016)

#Compute synchrony values
synch_warm_all=community_sync_Gross(db_warm_all,nrands=anrands,method=amethod_w)
synch_warm_pre_2006=community_sync_Gross(db_warm_pre_2006,nrands=anrands,method=amethod_w)
synch_warm_post_2006=community_sync_Gross(db_warm_post_2006,nrands=anrands,method=amethod_w)

synch_cold_all=community_sync_Gross(db_cold_all,nrands=anrands,method=amethod_w)
synch_cold_pre_2006=community_sync_Gross(db_cold_pre_2006,nrands=anrands,method=amethod_w)
synch_cold_post_2006=community_sync_Gross(db_cold_post_2006,nrands=anrands,method=amethod_w)

list_duck=list(synch_cold_all,synch_cold_pre_2006,synch_cold_post_2006,synch_warm_all,synch_warm_pre_2006,synch_warm_post_2006)

mat=rep(NA,length(list_duck))
for(v in 1:length(list_duck)){
        mat[v]=list_duck[[v]]$pval
}
mat_adj=p.adjust(mat,method=type_correct)
for(v in 1:length(list_duck)){
        list_duck[[v]]$pval=mat_adj[v]
}


print('Freq birds output')
upsi=0.02
color=rep(c("Black","Lightblue","Darkblue"),2)
filename_pdf=paste("Gross_freq_new_",end_bio,".pdf",sep="")
pdf(paste("OUT/",filename_pdf,sep=""),width=11,height=7,family="Times")
  #if(box_index){
#pdf(paste("Submission_JAE/Revisions/Gross_freq",type_correct,amethod_b,"boxplot.pdf",sep="_"),width=11,height=7)
#}else{
#pdf(paste("Submission_JAE/Revisions/Gross_freq",type_correct,amethod_b,"line.pdf",sep="_"),width=11,height=7) 
#}
par(mfrow=c(1,1),mar=c(3,4.5,2,.25),oma=c(1,2,1,.25),mgp=c(3,1,0),xpd=NA)
plot(0,0,t="n",ylim=c(-1,1),xlim=c(0,7.5),xaxt="n",xlab="",ylab="Synchrony index",cex.axis=1.8,cex.lab=1.8,main="",cex.main=2,las=1)
#ylim=c(0.15,0.35),
#mtext("Within",side=2,line=0.5,outer=T,cex=2)
axis(1,at=c(2,6),lab=c("Cold","Warm"),cex.axis=1.8)
for (v in 1:6){
        essai_taxo=list_freq
        obs=essai_taxo[[v]]$obs
 #       ic=2*sd(essai_taxo[[v]]$rands[1:100])/sqrt(100)
        p_s=essai_taxo[[v]]$pval
        x=v+(v>3)
#        yHigh=as.numeric(obs)+as.numeric(ic)
#        yLow=as.numeric(obs)-as.numeric(ic)
        points(x,obs,pch=21,col=color[v],bg=color[v],cex=2)
#        arrows(x,yHigh,x,yLow,angle=90,length=0.1,code=3,col=color[v],lwd=2)
	if(box_index){
        boxplot(essai_taxo[[v]]$rands[1:100],at=x,add=T,boxwex=0.25,range=0,xaxt="n",yaxt="n")
	}else{
        plou=quantile(essai_taxo[[v]]$rands[1:100],c(0.05,0.95))
        arrows(x,plou[1],x,plou[2],angle=90,length=0.1,code=3,col=color[v],lwd=2)
	}
        if (p_s<=thresh){
                points(x,as.numeric(obs)+upsi,pch='*',col="red",cex=2)
                }

        }
lines(c(0,7.5),c(0.,0.),lty=2,lwd=2)
legend("topleft",c("All","Pre-2006","Post-2006"),pch=NA,fill=c("black","Lightblue","Darkblue"),pt.cex=2,bty="n",cex=2)
legend("bottomright",c("Freq"),pch=c(21),pt.bg="black",pt.cex=2,bty="n",cex=2)
dev.off()
system(paste("cp OUT/",filename_pdf," Submission_JAE/Revisions/",filename_pdf,sep=""))

print('Actual Figure 2, within groups')
#PLOT
upsi=0.05
#pdf("OUT/Fig2_new2.pdf",width=11,height=14,family="Times")
filename_pdf=paste("Fig2_new3_JAE_",end_bio,".pdf",sep="")
pdf(paste("OUT/",filename_pdf,sep=""),width=11,height=14)
#if(box_index){
#pdf(paste("Submission_JAE/Revisions/Fig2_new2_JAE",type_correct,amethod_b,"boxplot.pdf",sep="_"),width=11,height=14)
#}else{
#pdf(paste("Submission_JAE/Revisions/Fig2_new2_JAE",type_correct,amethod_b,"line.pdf",sep="_"),width=11,height=14)
#}

par(mfrow=c(2,2),mar=c(3,4.5,2,.25),oma=c(1,2,1,.25),mgp=c(3,1,0),xpd=NA)
plot(0,0,t="n",ylim=c(-1,1),xlim=c(0,7.5),xaxt="n",xlab="",ylab="Synchrony index",cex.axis=1.8,cex.lab=1.8,main="Taxonomic groups",cex.main=2,las=1)
mtext("Within",side=2,line=0.3,outer=T,cex=2,font=2,at=0.75)
mtext("a)",side=2,line=-4,at=0.99,cex=2,outer=T,las=1)
axis(1,at=c(2,6),lab=c("Cold","Warm"),cex.axis=1.8)
for (v in 1:6){
        essai_taxo=list_anas
        obs=essai_taxo[[v]]$obs
#        ic=2*sd(essai_taxo[[v]]$rands[1:100])/sqrt(100)
        p_s=essai_taxo[[v]]$pval
        x=v+(v>3)
#        yHigh=as.numeric(obs)+as.numeric(ic)
#        yLow=as.numeric(obs)-as.numeric(ic)
        points(x,obs,pch=21,col=color[v],bg=color[v],cex=2)
#        arrows(x,yHigh,x,yLow,angle=90,length=0.1,code=3,col=color[v],lwd=2)
	if(box_index){
        boxplot(essai_taxo[[v]]$rands[1:100],at=x,add=T,boxwex=0.25,range=0,xaxt="n",yaxt="n")
	}else{
        plou=quantile(essai_taxo[[v]]$rands[1:100],c(0.05,0.95))
        arrows(x,plou[1],x,plou[2],angle=90,length=0.1,code=3,col=color[v],lwd=2)
	}
        if (p_s<=thresh){
                points(x,as.numeric(obs)+upsi,pch='*',col="red",cex=2)
                }

        essai_taxo=list_calidris
        obs=essai_taxo[[v]]$obs
#        ic=2*sd(essai_taxo[[v]]$rands[1:100])/sqrt(100)
        p_s=essai_taxo[[v]]$pval
        x=v+(v>3)
#        yHigh=as.numeric(obs)+as.numeric(ic)
#        yLow=as.numeric(obs)-as.numeric(ic)
        points(x+0.2,obs,pch=22,col=color[v],bg=color[v],cex=2)
#        arrows(x,yHigh,x,yLow,angle=90,length=0.1,code=3,col=color[v],lwd=2,lty=2)
	if(box_index){
        boxplot(essai_taxo[[v]]$rands[1:100],at=x+0.2,add=T,boxwex=0.25,range=0,xaxt="n",yaxt="n")
	}else{
        plou=quantile(essai_taxo[[v]]$rands[1:100],c(0.05,0.95))
        arrows(x+0.2,plou[1],x+0.2,plou[2],angle=90,length=0.1,code=3,col=color[v],lwd=2)
	}
        if (p_s<=thresh){
                points(x+0.2,as.numeric(obs)+upsi,pch='*',col="red",cex=2)
                }

        }
#ll=c(list_anas[[2]]$obs,list_anas[[3]]$obs)
#lines(2:3,ll,col="black",lwd=2,lty=1)
#ll=c(list_anas[[5]]$obs,list_anas[[6]]$obs)
#lines(6:7,ll,col="black",lwd=2,lty=1)

#ll=c(list_calidris[[2]]$obs,list_calidris[[3]]$obs)
#lines(2:3,ll,col="black",lwd=2,lty=2)
#ll=c(list_calidris[[5]]$obs,list_calidris[[6]]$obs)
#lines(6:7,ll,col="black",lwd=2,lty=2)


lines(c(0,7.5),c(0.,0.),lty=2,lwd=2)
legend("bottomleft",c("All","Pre-2006","Post-2006"),pch=NA,fill=c("black","Lightblue","Darkblue"),pt.cex=2,bty="n",cex=2)
legend("bottomright",c("Anatini","Calidris"),pch=c(21,22),pt.bg=c("black"),pt.cex=2,bty="n",cex=2)

plot(0,0,t="n",ylim=c(-1,1),xlim=c(0,7.5),xaxt="n",xlab="",ylab="",yaxt="n",main='Functional groups',cex.main=2,las=1)
mtext("b)",side=2,line=-35,at=0.99,cex=2,outer=T,las=1)

axis(1,at=c(2,6),lab=c("Cold","Warm"),cex.axis=1.8)
for (v in 1:6){
        essai_taxo=list_waders
        obs=essai_taxo[[v]]$obs
#        ic=2*sd(essai_taxo[[v]]$rands[1:100])/sqrt(100)
        p_s=essai_taxo[[v]]$pval
        x=v+(v>3)
#        yHigh=as.numeric(obs)+as.numeric(ic)
#        yLow=as.numeric(obs)-as.numeric(ic)
        points(x,obs,pch=23,col=color[v],bg=color[v],cex=2)
#        arrows(x,yHigh,x,yLow,angle=90,length=0.1,code=3,col=color[v],lwd=2)
	if(box_index){
        boxplot(essai_taxo[[v]]$rands[1:100],at=x,add=T,boxwex=0.25,range=0,xaxt="n",yaxt="n")
        }else{
	plou=quantile(essai_taxo[[v]]$rands[1:100],c(0.05,0.95))
        arrows(x,plou[1],x,plou[2],angle=90,length=0.1,code=3,col=color[v],lwd=2)
	}
        if (p_s<=thresh){
                points(x,as.numeric(obs)+upsi,pch='*',col="red",cex=2)
                }

        essai_taxo=list_duck
        obs=essai_taxo[[v]]$obs
 #       ic=2*sd(essai_taxo[[v]]$rands[1:100])/sqrt(100)
        p_s=essai_taxo[[v]]$pval
        x=v+(v>3)
#        yHigh=as.numeric(obs)+as.numeric(ic)
#        yLow=as.numeric(obs)-as.numeric(ic)
        points(x+0.2,obs,pch=24,col=color[v],bg=color[v],cex=2)
#        arrows(x,yHigh,x,yLow,angle=90,length=0.1,code=3,col=color[v],lwd=2,lt=2)
	if(box_index){
        boxplot(essai_taxo[[v]]$rands[1:100],at=x+0.2,add=T,boxwex=0.25,range=0,yaxt="n",xaxt="n")
	}else{
        plou=quantile(essai_taxo[[v]]$rands[1:100],c(0.05,0.95))
        arrows(x+0.2,plou[1],x+0.2,plou[2],angle=90,length=0.1,code=3,col=color[v],lwd=2)
	}
        if (p_s<=thresh){
                points(x+0.2,as.numeric(obs)+upsi,pch='*',col="red",cex=2)
                }

        }

#ll=c(list_waders[[2]]$obs,list_waders[[3]]$obs)
#lines(2:3,ll,col="black",lwd=2,lty=1)
#ll=c(list_waders[[5]]$obs,list_waders[[6]]$obs)
#lines(6:7,ll,col="black",lwd=2,lty=1)

#ll=c(list_duck[[2]]$obs,list_duck[[3]]$obs)
#lines(2:3,ll,col="black",lwd=2,lty=2)
#ll=c(list_duck[[5]]$obs,list_duck[[6]]$obs)
#lines(6:7,ll,col="black",lwd=2,lty=2)

lines(c(0,7.5),c(0.,0.),lty=2,lwd=2)
legend("bottomright",c("Waders","Ducks"),pch=c(23,24),pt.bg="black",pt.cex=2,bty="n",cex=2)

print('Now computing between groups')
#Anas and Calidris
db_warm=read.csv(paste("IN/warmseason_",end_bio,"_asdataframe_summed_v2_wtoutrarespecies.csv",sep=""),sep=";",header=T)
db_cold=read.csv(paste("IN/coldseason_",end_bio,"_asdataframe_summed_v2_wtoutrarespecies.csv",sep=""),sep=";",header=T)
#db_warm=read.csv("IN/warmseason_abundances_asdataframe_summed_biomasses.csv",sep=";",header=T)
#db_cold=read.csv("IN/coldseason_abundances_asdataframe_summed_biomasses.csv",sep=";",header=T)

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


color=rep(c("Black","Lightblue","Darkblue"),2)
#par(mfrow=c(1,2),mar=c(3,4.5,1,.25),oma=c(1,2,1,.25),mgp=c(3,1,0),xpd=NA)
plot(0,0,t="n",ylim=c(-1,1),xlim=c(0,7.5),xaxt="n",xlab="",ylab="Synchrony index",cex.lab=1.8,cex.axis=1.8,las=1)
axis(1,at=c(2,6),lab=c("Cold","Warm"),cex.axis=1.8)
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
legend("topleft",c("Anatini/Calidris"),pch=21,pt.bg=c("black"),pt.cex=2,bty="n",cex=2)

db_warm_all=subset(db_warm,sp_data_frame %in% c("Waders","Ducks")&dates<2016)
db_warm_pre_2006=subset(db_warm,sp_data_frame %in% c("Waders","Ducks")&dates<=2006)
db_warm_post_2006=subset(db_warm,sp_data_frame %in% c("Waders","Ducks")&dates>2006&dates<2016)
db_cold_all=subset(db_cold,sp_data_frame %in% c("Waders","Ducks")&dates<2016)
db_cold_pre_2006=subset(db_cold,sp_data_frame %in% c("Waders","Ducks")&dates<=2006)
db_cold_post_2006=subset(db_cold,sp_data_frame %in% c("Waders","Ducks")&dates>2006&dates<2016) #Because there is a problem for 2016 (na values)

#Compute synchrony values
synch_warm_all=community_sync_Gross(db_warm_all,nrands=anrands,method=amethod_b)
synch_warm_pre_2006=community_sync_Gross(db_warm_pre_2006,nrands=anrands,method=amethod_b)
synch_warm_post_2006=community_sync_Gross(db_warm_post_2006,nrands=anrands,method=amethod_b)

synch_cold_all=community_sync_Gross(db_cold_all,nrands=anrands,method=amethod_b)
synch_cold_pre_2006=community_sync_Gross(db_cold_pre_2006,nrands=anrands,method=amethod_b)
synch_cold_post_2006=community_sync_Gross(db_cold_post_2006,nrands=anrands,method=amethod_b)

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
plot(0,0,t="n",ylim=c(-1,1),xlim=c(0,7.5),xaxt="n",xlab="",ylab="",yaxt="n",las=1)
axis(1,at=c(2,6),lab=c("Cold","Warm"),cex.axis=2)
mtext("d)",side=2,line=-35,at=0.49,cex=2,outer=T,las=1)
for (v in 1:6){
        obs=essai_taxo[[v]]$obs
#        ic=2*sd(essai_taxo[[v]]$rands[1:100])/sqrt(100)
        p_s=essai_taxo[[v]]$pval
        x=v+(v>3)
#        yHigh=as.numeric(obs)+as.numeric(ic)
#        yLow=as.numeric(obs)-as.numeric(ic)
        points(x,obs,pch=22,col=color[v],bg=color[v],cex=2)
 #       arrows(x,yHigh,x,yLow,angle=90,length=0.1,code=3,col=color[v],lwd=2)
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
#ll=c(essai_taxo[[2]]$obs,essai_taxo[[3]]$obs)
#lines(2:3,ll,col="black",lwd=2,lty=1)
#ll=c(essai_taxo[[5]]$obs,essai_taxo[[6]]$obs)
#lines(6:7,ll,col="black",lwd=2,lty=1)

legend("topleft",c("Ducks/Waders"),pch=c(22),pt.bg=c("black"),pt.cex=2,bty="n",cex=2)
dev.off()

system(paste("cp OUT/",filename_pdf," Submission_JAE/Revisions/",filename_pdf,sep=""))
