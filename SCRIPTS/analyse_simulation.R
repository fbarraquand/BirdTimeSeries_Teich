#CP 2018 
#Analyse_simulation

rm(list=ls())
graphics.off()
source("SCRIPTS/test_synchrony_Gross.r")

#Parameters which change
pdf("OUT/analyse_simulation.pdf")
n_species_series=c(5,15,30,50)
n_time_series=c(35,100,500)
nb_repeats=1:100
tab_synch=array(NA,dim=c(length(n_species_series),length(n_time_series),length(nb_repeats)))
par(mfrow=c(3,1))
j=0
for(n_time in n_time_series){
	j=j+1
	i=0
	plot(0,0,t="n",xlim=c(0,5),ylim=c(-1,1),xlab="Nb_species",xaxt="n",ylab="Synchrony index",main=paste("Length=",n_time,sep=""))
	axis(1,at=c(1:4),lab=c("10","30","60","100"))
	for(n_species in n_species_series){
		i=i+1
		print(paste(n_time,n_species))
		tab=read.table(paste("OUT/MockBirdsTimeSeries_species",n_species,"_time",n_time,".csv",sep=""),header=T,sep=",",dec=".")
		#nb_repeats=unique(tab$Repeat)
		tmp_synch=c()
		for(k in nb_repeats){
			data=subset(tab,Repeat==k)
			dates=as.numeric(sort(rep(1:n_time,2*n_species)))
			sp_data_frame=rep(paste("X",1:(2*n_species),sep=""),n_time)
			abundance=as.vector(t(data[,4:(4+n_species*2-1)]))
			community=data.frame(dates,sp_data_frame,abundance)
			tab_synch[i,j,k]=c(tmp_synch,community_sync_Gross(community,nrands=1)$obs)
		}
		boxplot(at=i,x=tab_synch[i,j,],add=T)
	}
}
dev.off()

#load OUT/tab_synch (for 10)
pdf("OUT/analyse_simulation_100repeats.pdf")


par(mfrow=c(3,1))
i=0
for(n_time in n_time_series){
        j=0
	i=i+1
	plot(0,0,t="n",xlim=c(0,5),ylim=c(-0.5,0.5),xlab="Nb_species",xaxt="n",ylab="Synchrony index",main=paste("Length=",n_time,sep=""))
	axis(1,at=c(1:4),lab=c("10","30","60","100"))
	for(n_species in n_species_series){
        	j=j+1
                boxplot(at=j,x=tab_synch[j,i,],add=T)
        }
	abline(h=0.0,lty=2)
}


par(mfrow=c(2,2))
i=0
for(n_species in n_species_series){
	i=i+1
	j=0
	plot(0,0,t="n",xlim=c(0,4),ylim=c(-0.5,0.5),xlab="Length time series",xaxt="n",ylab="Synchrony index",main=paste("Nb species=",n_species,sep=""))
	axis(1,at=c(1:3),lab=c("35","100","500"))
	for(n_time in n_time_series){
		j=j+1
		boxplot(at=j,x=tab_synch[i,j,],add=T)
	}
	abline(h=0.0,lty=2)
}

dev.off()
