####### CP 2020/03/07 Checking time series

rm(list=ls())
graphics.off()

library("RColorBrewer")

end_of_file_seq=c("4sp_pair","40sp_pair","4sp_alpha15_beta15","40sp_alpha15_beta15","4sp_alpha2_beta4","40sp_alpha2_beta4","4sp_alpha4_beta2","40sp_alpha4_beta2")
explain=c("rho=-0.8","rho=-0.8","quasi-normal","quasi-normal","compensation","compensation","synchrony","synchrony")

rep=1

#e=1

lis_hist=list()

for(e in 1:length(end_of_file_seq)){
end_of_file=end_of_file_seq[e]
#############################
tab_bm=read.table(paste("MockData_SAD_",end_of_file,".csv",sep=""),sep=";",dec=".",header=T)
tab_bm$Species=as.character(tab_bm$Species)
tab_bm$Time=as.numeric(as.character(tab_bm$Time))
tab_bm$Abundance=as.numeric(as.character(tab_bm$Abundance))

dates_tmp=unique(tab_bm$Time)

tab_bm=subset(tab_bm,Rep==1)

species=unique(tab_bm$Species)

if(length(species)==4){
	col=c("red","blue","black","orange")
	nb=4
}else{
	col=colorRampPalette(rev(brewer.pal(length(species), "Spectral")))(length(species))
	nb=10
}

tab_species=matrix(0,nrow=length(dates_tmp),ncol=length(species))
colnames(tab_species)=species

for(id in 1:length(dates_tmp)){
        for (s in species){
                id_d=which(tab_bm$Time==dates_tmp[id]&tab_bm$Species==s)
                if(length(id_d)>0){
                        tab_species[id,s]=tab_bm$Abundance[id_d]
                }
        }
}
mean_val=apply(tab_species,2,mean)

lis_hist[[e]]=hist(mean_val,breaks=nb,plot=F,main=explain[e])

pdf(paste("time_series_",end_of_file,".pdf",sep=""),height=12,width=12)
par(mfrow=c(2,1))
plot(dates_tmp,rep(0,length(dates_tmp)),t="n",ylim=range(tab_bm$Abundance),ylab="Abundance",xlab="Date",main=explain[e])
for(s in 1:length(species)){
	lines(dates_tmp,tab_species[,species[s]],col=col[s])
}
if(length(species)==4){
legend("topleft",legend=1:4,col=col,lty=1)
}

tab_species_norm=tab_species
for(s in species){
        tab_species_norm[,s]=scale(tab_species[,s])
}
plot(dates_tmp,rep(0,length(dates_tmp)),t="n",ylim=range(c(-2,2)),ylab="Abundance",xlab="Date",main=explain[e])
for(s in 1:length(species)){
	lines(dates_tmp,tab_species_norm[,species[s]],col=col[s])
}



dev.off()
}

library("corrplot")
pdf("distrib_average_abundance_simulation.pdf",height=15)
par(mfrow=c(4,2))
se_4=seq(1,length(end_of_file_seq),by=2)
se_40=seq(2,length(end_of_file_seq),by=2)
for(e in se_4){
	plot(lis_hist[[e]],main=explain[e])
	tab=read.table(paste("MuSigma_SAD_",end_of_file_seq[e],".csv",sep=""),sep=";",header=T,dec=".")
	Sigma=as.matrix(tab[,2:ncol(tab)])
	rownames(Sigma)=colnames(Sigma)
	corrplot(Sigma,is.corr=F)
}
for(e in se_40){
        plot(lis_hist[[e]],main=explain[e])
        tab=read.table(paste("MuSigma_SAD_",end_of_file_seq[e],".csv",sep=""),sep=";",header=T,dec=".")
	Sigma=as.matrix(tab[,2:ncol(tab)])
	rownames(Sigma)=colnames(Sigma)
	corrplot(Sigma,is.corr=F)
}
dev.off()

