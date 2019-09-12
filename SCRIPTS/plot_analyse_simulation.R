#CP 2018
#This script uses tab_synch*.RData to plot the variation of the Gross index with the length of the time series and the size of the community
#2019/09/12 Changed the number of species we look at to be closer to what we have in the data
graphics.off()
rm(list=ls())

tab_list=c("notrend","subsample_notrend",'subsample_withtrend')
colors=c("lightblue","blue","darkblue")
b_series=c(0.1,0.5,0.75)
#n_species_series=c(5,15,30,50)
n_species_series=c(1,5,15,30)
n_time_series=c(35,100,500)

for(id in 1:length(tab_list)){
pdf(paste("OUT/simulated_gross_values_",tab_list[id],"_20190912.pdf",sep=''),width=12,height=12)
load(paste("OUT/tab_synch_",tab_list[id],"_withburnin_20190912.RData",sep=""))
par(mfrow=c(2,2),mar=c(4,5,3,0.5))
i=0
for(n_species in n_species_series){
        i=i+1
        j=0
        plot(0,0,t="n",xlim=c(0,9),ylim=c(-.9,0.2),xlab="Length time series",xaxt="n",ylab="Synchrony index",main=paste("Nb species=",2*n_species,sep=""),cex=2,cex.axis=2,cex.lab=2,cex.main=2)
        axis(1,at=seq(1,9,3),lab=c("35","100","500"),cex.axis=2,cex.lab=2)
        for(n_time in n_time_series){
                j=j+1
                for(l in 1:length(b_series)){
                        boxplot(at=(j-1)*length(b_series)+l*0.5,x=tab_synch[i,j,l,],add=T,col=colors[l],yaxt="n")
                }
        }
        abline(h=0.0,lty=2)
}
legend("bottomright",paste("b=",b_series,sep=""),bty="n",fill=colors,cex=2)
dev.off()
}
