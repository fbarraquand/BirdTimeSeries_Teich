#CP 2018 
#This script computes the Gross index without using codyn, and compare with the codyn result

rm(list=ls())
graphics.off()

source("SCRIPTS/test_synchrony_Gross.r")

compute_eta=function(tab){
	#assuming that we have the same structure than codyn, ie a data frame with time, abundance, and sp_data frame
	sp=unique(tab$sp_data_frame)
	eta=0
	for(i in 1:length(sp)){
		sum_j=aggregate(tab$abundance[tab$sp_data_frame!=sp[i]],list(tab$dates[tab$sp_data_frame!=sp[i]]),sum)
			eta=eta+cor(tab$abundance[tab$sp_data_frame==sp[i]],sum_j$x)
	}
	eta=eta/length(sp)
	return(eta)	
}

#tab_synch=array(NA,dim=c(length(n_species_series),length(n_time_series),length(b_series),length(nb_repeats)))
n_species=5
n_time=35
b=0.1

		tab=read.table(paste("OUT/MockBirdsTimeSeries_species",n_species,"_time",n_time,"_b",b,".csv",sep=""),header=T,sep=",",dec=".")
                #nb_repeats=unique(tab$Repeat)
                        k=1
			data=subset(tab,Repeat==k)
                        dates=as.numeric(sort(rep(1:n_time,2*n_species)))
                        sp_data_frame=rep(paste("X",1:(2*n_species),sep=""),n_time)
                        abundance=as.vector(t(data[,4:(4+n_species*2-1)]))
                        community=data.frame(dates,sp_data_frame,abundance)

			print(compute_eta(community))
			print(community_sync_Gross(community,nrands=0)$obs)

n_species=50
n_time=500
b=0.75

                tab=read.table(paste("OUT/MockBirdsTimeSeries_species",n_species,"_time",n_time,"_b",b,".csv",sep=""),header=T,sep=",",dec=".")
                #nb_repeats=unique(tab$Repeat)
                        k=1
                        data=subset(tab,Repeat==k)
                        dates=as.numeric(sort(rep(1:n_time,2*n_species)))
                        sp_data_frame=rep(paste("X",1:(2*n_species),sep=""),n_time)
                        abundance=as.vector(t(data[,4:(4+n_species*2-1)]))
                        community=data.frame(dates,sp_data_frame,abundance)

                        print(compute_eta(community))
                        print(community_sync_Gross(community,nrands=0)$obs)
