#CP 2018 
#This script to check some of the simulations (just plot) and to verify the index computed by codyn

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
n_time_series=c(35,500)
b_series=c(0.1,0.5,0.75)
pdf("OUT/verif_visual_simulation_v2.pdf",width=14,height=12)
                for(b in b_series){
                        tab=read.table(paste("OUT/MockBirdsTimeSeries_species",n_species,"_time",n_time,"_b",b,".csv",sep=""),header=T,sep=",",dec=".")
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
dev.off()
}

n_species=5
        id_species_1=paste("X",1:5,sep="")
        id_species_2=paste("X",6:10,sep="")
n_time_series=c(35,500)
pdf("OUT/times_series_esm.pdf")

        par(mfrow=c(2,1),mar=c(2,2,1,2))
for(n_time in n_time_series){
	b=0.75
	tab=read.table(paste("OUT/MockBirdsTimeSeries_species",n_species,"_time",n_time,"_b",b,".csv",sep=""),header=T,sep=",",dec=".")
	k=7
                                plot(0,0,t="n",xlim=c(0,n_time),ylab="",xlab="",ylim=c(0.02,0.45))
                                tabbis=subset(tab,Repeat==k)
                                for(l in 1:5){
                                        lines(tabbis[,id_species_1[l]],t="o",col=cold_colors[l],pch=16)
                                        lines(tabbis[,id_species_2[l]],t="o",col=hot_colors[l],pch=16)
                                }
                                par(new=TRUE)
                                plot(tabbis$Abiotic_var1,t="l",col="black",xlab="",yaxt="n",lwd=2,lty=2,xlim=c(0,n_time),ylab="")
                                axis(4)
	}
dev.off()

