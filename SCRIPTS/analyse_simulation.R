#CP 2018 
#Analyse_simulation

rm(list=ls())
graphics.off()

subsample=TRUE

source("SCRIPTS/test_synchrony_Gross.r")

#Parameters which change
#pdf("OUT/analyse_simulation.pdf")
n_species_series=c(5,15,30,50)
n_time_series=c(35,100,500)
b_series=c(0.1,0.5,0.75)
nb_repeats=1:10
tab_synch=array(NA,dim=c(length(n_species_series),length(n_time_series),length(b_series),length(nb_repeats)))
#par(mfrow=c(3,1))
j=0
for(n_time in n_time_series){
	j=j+1
	i=0
#	plot(0,0,t="n",xlim=c(0,5),ylim=c(-1,1),xlab="Nb_species",xaxt="n",ylab="Synchrony index",main=paste("Length=",n_time,sep=""))
#	axis(1,at=c(1:4),lab=c("10","30","60","100"))
	for(n_species in n_species_series){
		l=0
		i=i+1
		for(b in b_series){
		l=l+1
		print(paste(n_time,n_species,b))
		if(subsample){		
			tab=read.table(paste("OUT/MockBirdsTimeSeries_species",n_species,"_time500_b",b,"_withtrend_withburnin.csv",sep=""),header=T,sep=",",dec=".")
			id=floor(seq(1,500,length.out=n_time))
			tab=subset(tab,Time_index %in% id)
			
		}else{
			tab=read.table(paste("OUT/MockBirdsTimeSeries_species",n_species,"_time",n_time,"_b",b,"_withtrend_withburnin.csv",sep=""),header=T,sep=",",dec=".")
		}
		#nb_repeats=unique(tab$Repeat)
		for(k in nb_repeats){
			data=subset(tab,Repeat==k)
			dates=as.numeric(sort(rep(1:n_time,2*n_species)))
			sp_data_frame=rep(paste("X",1:(2*n_species),sep=""),n_time)
			abundance=as.vector(t(data[,4:(4+n_species*2-1)]))
			community=data.frame(dates,sp_data_frame,abundance)
			tab_synch[i,j,l,k]=community_sync_Gross(community,nrands=0)$obs
		}
#		boxplot(at=i,x=tab_synch[i,j,],add=T)
	}
}
}
#dev.off()


#load OUT/tab_synch (for 10)
if(subsample){
pdf("OUT/analyse_simulation_varyingb_v2_subsample_withtrend_withburnin.pdf",height=12,width=11)
save(tab_synch,file="tab_synch_subsample_withtrend_withburnin.RData")
}else{
pdf("OUT/analyse_simulation_varyingb_withtrend_withburnin.pdf",height=12,width=11)
save(tab_synch,file="tab_synch_withtrend_withburnin.RData")
}
colors=c("lightblue","blue","darkblue")
par(mfrow=c(3,1),mar=c(4,6,3,1))
i=0
for(n_time in n_time_series){
	i=i+1
	plot(0,0,t="n",xlim=c(0,11),ylim=c(-0.5,0.5),xlab="Nb_species",xaxt="n",ylab="Synchrony index",main=paste("Length=",n_time,sep=""),cex=2,cex.axis=2,cex.lab=2,cex.main=2)
	axis(1,at=seq(1,12,3),lab=c("10","30","60","100"),cex.axis=2,cex.lab=2)
        j=0
	for(n_species in n_species_series){
        	j=j+1
		for(l in 1:length(b_series)){
                	boxplot(at=(j-1)*length(b_series)+l*0.5,x=tab_synch[j,i,l,],add=T,col=colors[l],yaxt="n")
		}
        }
	abline(h=0.0,lty=2)
}
legend("topleft",paste("b=",b_series,sep=""),bty="n",fill=colors,cex=2)


par(mfrow=c(2,2),mar=c(4,5,3,0.5))
i=0
for(n_species in n_species_series){
	i=i+1
	j=0
	plot(0,0,t="n",xlim=c(0,9),ylim=c(-0.5,0.5),xlab="Length time series",xaxt="n",ylab="Synchrony index",main=paste("Nb species=",2*n_species,sep=""),cex=2,cex.axis=2,cex.lab=2,cex.main=2)
	axis(1,at=seq(1,9,3),lab=c("35","100","500"),cex.axis=2,cex.lab=2)
	for(n_time in n_time_series){
		j=j+1
                for(l in 1:length(b_series)){
                        boxplot(at=(j-1)*length(b_series)+l*0.5,x=tab_synch[i,j,l,],add=T,col=colors[l],yaxt="n")
                }
	}
	abline(h=0.0,lty=2)
}
legend("topright",paste("b=",b_series,sep=""),bty="n",fill=colors,cex=2)
dev.off()
