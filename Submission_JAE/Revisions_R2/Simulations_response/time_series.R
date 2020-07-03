####### CP 2020/03/07 Checking time series

rm(list=ls())
graphics.off()

end_of_file_seq=c("4sp_alpha4_beta2")#,"40sp_pair","4sp_alpha15_beta15","40sp_alpha15_beta15","4sp_alpha2_beta4","40sp_alpha2_beta4","4sp_alpha4_beta2","40sp_alpha4_beta2")
explain=c("rho=-0.8","rho=-0.8","quasi-normal","quasi-normal","compensation","compensation","synchrony","synchrony")

rep=1

e=1
#for(e in 1:length(end_of_file_seq)){
end_of_file=end_of_file_seq[e]
#############################
tab_bm=read.table(paste("MockData_SAD_",end_of_file,".csv",sep=""),sep=";",dec=".",header=T)
tab_bm$Species=as.character(tab_bm$Species)
tab_bm$Time=as.numeric(as.character(tab_bm$Time))
tab_bm$Abundance=as.numeric(as.character(tab_bm$Abundance))

dates_tmp=unique(tab_bm$Time)

tab_bm=subset(tab_bm,Rep==1)


species=unique(tab_bm$Species)

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

col=c("red","blue","black","orange")
pdf(paste("time_series_",end_of_file,".pdf",sep=""))
plot(dates_tmp,rep(0,length(dates_tmp)),t="n",ylim=range(tab_bm$Abundance),ylab="Abundance",xlab="Date")
for(s in 1:length(species)){
	lines(dates_tmp,tab_species[,species[s]],col=col[s])
}
dev.off()
