#CP 2018 
#This script to check some of the simulations (just plot)

rm(list=ls())
graphics.off()
source("SCRIPTS/test_synchrony_Gross.r")


#tab_synch=array(NA,dim=c(length(n_species_series),length(n_time_series),length(b_series),length(nb_repeats)))
n_species_series=c(5,15,30,50)
n_time_series=c(35,100,500)
b_series=c(0.1,0.5,0.75)
nb_repeats=1:10
		
#Doesn't make sense to have a look at 100 species at a time, so I'm just gonna sample 5 species per subgroup for all communities
cold_colors=c("cyan","darkblue","grey","darkslateblue","green")
hot_colors=c("red","orange","darkgoldenrod1","deeppink","lightsalmon")
if(1==0){ #I will take that again later
pdf("OUT/verif_visual_simulation.pdf")
for(n_species in n_species_series){
	id_species_1=paste("X",sample(1:n_species,5),sep="")
	id_species_2=paste("X",sample((n_species+1):(2*n_species),5),sep="")
	for(n_time in n_time_series){
                for(b in b_series){
			tab=read.table(paste("OUT/MockBirdsTimeSeries_species",n_species,"_time",n_time,"_b",b,".csv",sep=""),header=T,sep=",",dec=".")
			par(mfrow=c(4,3))
			for(k in nb_repeats){
				plot(0,0,t="n",ylim=c(0,2),xlim=c(0,n_time),main=paste("species",n_species,"_time",n_time,"_b",b,sep=""))
				tabbis=subset(tab,Repeat==k)
				for(l in 1:5){
					lines(tabbis[,id_species_1[l]],t="o",col=cold_colors[l])
					lines(tabbis[,id_species_2[l]],t="o",col=hot_colors[l])
				}
			}
		}
	}
}
dev.off()
}

if(1==0){
n_species=5
	id_species_1=paste("X",1:5,sep="")
	id_species_2=paste("X",6:10,sep="")
n_time_series=c(500)
b_series=c(0.1,0.5,0.75)
end=c("withtrend","notrend")
for(e in end){
pdf(paste("OUT/verif_visual_simulation_withburnin_",e,".pdf",sep=""),width=14,height=12)
	for(n_time in n_time_series){
                for(b in b_series){
                        tab=read.table(paste("OUT/MockBirdsTimeSeries_species",n_species,"_time",n_time,"_b",b,"_",e,"_withburnin.csv",sep=""),header=T,sep=",",dec=".")
                        par(mfrow=c(4,3))
                        for(k in nb_repeats){
				if(n_time==35){
					ymax=0.4
				}else{
					ymax=0.34
				}
                                plot(0,0,t="n",ylim=c(0,ymax),xlim=c(0,n_time),main=paste("time",n_time,"_b",b,sep=""))
                                tabbis=subset(tab,Repeat==k)
                                for(l in 1:5){
                                        lines(tabbis[,id_species_1[l]],t="o",col=cold_colors[l])
                                        lines(tabbis[,id_species_2[l]],t="o",col=hot_colors[l])
                                }
				par(new=TRUE)
				print(n_time)
				plot(tabbis$Abiotic_var1,t="l",col="black",ylim=c(-0.,1),yaxt="n",lwd=2,xlim=c(0,n_time),ylab="")
				print(sd(tabbis$Abiotic_var1))
				axis(4)
                        }
                }
	}
dev.off()
}
}

library('RColorBrewer')
n_species=5
cold_colors=brewer.pal(n_species,'Blues')
warm_colors=brewer.pal(n_species,'Reds')
        id_species_1=paste("X",1:5,sep="")
        id_species_2=paste("X",6:10,sep="")
n_time_series=c(500)
	b=0.75
pdf("OUT/times_series_esm.pdf",width=15)

        par(mfrow=c(2,1),mar=c(2,3,1,3),xpd=NA,oma=c(2,2.5,1,2))
end=c("withtrend","notrend")
k=9
for(e in end){
for(n_time in n_time_series){
        tab=read.table(paste("OUT/MockBirdsTimeSeries_species",n_species,"_time",n_time,"_b",b,"_",e,"_withburnin.csv",sep=""),header=T,sep=",",dec=".")
        tabbis=subset(tab,Repeat==k)
	id=floor(seq(1,500,length.out=35))
        tabbis_sample=subset(tabbis,Time_index %in% id)
                                plot(0,0,t="n",xlim=c(0,n_time),xlab="",ylim=c(0.02,0.45),xaxt="n",ylab="",cex.axis=2)
                                for(l in 1:5){
                                        lines(tabbis[,id_species_1[l]],t="o",col=cold_colors[l],pch=1)
                                        points(id,tabbis_sample[,id_species_1[l]],t="p",bg=cold_colors[l],pch=21,cex=2,col="black")
                                        lines(tabbis[,id_species_2[l]],t="o",col=warm_colors[l],pch=1)
                                        points(id,tabbis_sample[,id_species_2[l]],t="p",bg=warm_colors[l],pch=21,cex=2,col="black")
                                }
                                par(new=TRUE)
                                plot(tabbis$Abiotic_var1,t="l",col="black",xlab="",yaxt="n",lwd=2,lty=2,xlim=c(0,n_time),ylab="",xaxt="n")
                                axis(4,cex.axis=1.5)
	}
}
axis(1,cex.axis=1.5)
mtext('Time',1,line=1,outer=T,cex=2)
mtext("a)",2,at=0.95,line=-2,outer=T,cex=2,las=1)
mtext("b)",2,at=0.47,line=-2,outer=T,cex=2,las=1)
mtext("Environmental driver",4,line=1,outer=T,cex=2)
mtext("Abundance",2,line=1,outer=T,cex=2)
dev.off()

plot(0,0,t="n",xlim=c(0,n_time),ylab="",xlab="",ylim=c(-1.0,1.0))
colo=c("grey","black")
i=0
for(e in end){
	i=i+1
        tab=read.table(paste("OUT/MockBirdsTimeSeries_species",n_species,"_time",n_time,"_b",b,"_",e,"_withburnin.csv",sep=""),header=T,sep=",",dec=".")
        tabbis=subset(tab,Repeat==k)
        lines(tabbis$Abiotic_var1,t="l",col=colo[i],xlab="")
}
