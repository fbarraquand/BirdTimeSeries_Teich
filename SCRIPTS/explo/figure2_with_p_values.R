# Script by CPicoche 2018
#The aim is to compare synchrony indices at different taxonomic levels
library("R.utils")
rm(list=ls())
graphics.off()


freq=loadToEnv("OUT/indices_freq_with_p_values.RData")
lim=loadToEnv("OUT/indices_limicoles_with_p_values.RData")
anas=loadToEnv("OUT/indices_Anas_with_p_values.RData")
calidris=loadToEnv("OUT/indices_Calidris_with_p_values.RData")

id_obs=c(1,5,9)
id_rands=id_obs+1
id_p=id_obs+2
thresh=0.05

pdf("tmp_figure2.pdf",width=12.5)
layout(matrix(c(1,1,3,2,2,4),2,3,byrow=T))
par(mar=c(4,5,0.5,0.5))
type_of_season="data_freq_season4_Gross"
plot(0,0,t="n",ylim=c(0.15,0.75),xlim=c(0,15),xaxt="n",xlab="",ylab="Synchrony index",cex.axis=2,cex.lab=2)
axis(1,at=c(2,6,10,14),lab=rep("",4),cex.axis=2)
#for (type_of_season in c("data_freq_season4_Gross","data_freq_season2_Gross")){
	for (season in 1:4){
		obs=calidris[[type_of_season]][[season]][id_obs]
		ic=lapply(calidris[[type_of_season]][[season]][id_rands],function(x) 2*sd(x[1:100])/sqrt(100))
		p_s=calidris[[type_of_season]][[season]][id_p]
		x=1:3+(season-1)*4
		yHigh=as.numeric(obs)+as.numeric(ic)
		yLow=as.numeric(obs)-as.numeric(ic)
		points(x,obs,t="o",lty=1,pch=21,col=c("Black","Lightblue","Darkblue"),bg=c("Black","Lightblue","Darkblue"),cex=2)
		arrows(x,yHigh,x,yLow,angle=90,length=0.1,code=3,col=c("Black","Lightblue","Darkblue"),lwd=2)
		for (i in 1:3){
			if (p_s[i]<thresh){
				points(x[i],as.numeric(obs[i]),pch='*',col="red",cex=2)
			}
		}
	}

        for (season in 1:4){
                obs=anas[[type_of_season]][[season]][id_obs]
                ic=lapply(anas[[type_of_season]][[season]][id_rands],function(x) 2*sd(x[1:100])/sqrt(100))
                p_s=anas[[type_of_season]][[season]][id_p]
                x=1:3+(season-1)*4
                yHigh=as.numeric(obs)+as.numeric(ic)
                yLow=as.numeric(obs)-as.numeric(ic)
                points(x,obs,t="o",lty=1,pch=22,col=c("Black","Lightblue","Darkblue"),bg=c("Black","Lightblue","Darkblue"),cex=2)
                arrows(x,yHigh,x,yLow,angle=90,length=0.1,code=3,col=c("Black","Lightblue","Darkblue"),lwd=2)
                for (i in 1:3){
                        if (p_s[i]<thresh){
                                points(x[i],as.numeric(obs[i]),pch='*',col="red",cex=2)
                        }
                }
        }

par(mar=c(4,5,0.5,0.5))
plot(0,0,t="n",ylim=c(0.15,0.75),xlim=c(0,15),xaxt="n",xlab="",ylab="Synchrony index",cex.axis=2,cex.lab=2)
axis(1,at=c(2,6,10,14),lab=c("Winter","Spring","Summer","Autumn"),cex.axis=2)
        for (season in 1:4){
                obs=lim[[type_of_season]][[season]][id_obs]
                ic=lapply(lim[[type_of_season]][[season]][id_rands],function(x) 2*sd(x[1:100])/sqrt(100))
                p_s=lim[[type_of_season]][[season]][id_p]
                x=1:3+(season-1)*4
                yHigh=as.numeric(obs)+as.numeric(ic)
                yLow=as.numeric(obs)-as.numeric(ic)
                points(x,obs,pch=23,t="o",lty=1,col=c("Black","Lightblue","Darkblue"),bg=c("Black","Lightblue","Darkblue"),cex=2)
                arrows(x,yHigh,x,yLow,angle=90,length=0.1,code=3,col=c("Black","Lightblue","Darkblue"),lwd=2)
                for (i in 1:3){
                        if (p_s[i]<thresh){
                                points(x[i],as.numeric(obs[i]),pch='*',col="red",cex=2)
                        }
                }
        }

        for (season in 1:4){
                obs=freq[[type_of_season]][[season]][id_obs]
                ic=lapply(freq[[type_of_season]][[season]][id_rands],function(x) 2*sd(x[1:100])/sqrt(100))
                p_s=freq[[type_of_season]][[season]][id_p]
                x=1:3+(season-1)*4
                yHigh=as.numeric(obs)+as.numeric(ic)
                yLow=as.numeric(obs)-as.numeric(ic)
                points(x,obs,pch=24,t="o",lty=1,col=c("Black","Lightblue","Darkblue"),bg=c("Black","Lightblue","Darkblue"),cex=2)
                arrows(x,yHigh,x,yLow,angle=90,length=0.1,code=3,col=c("Black","Lightblue","Darkblue"),lwd=2)
                for (i in 1:3){
                        if (p_s[i]<thresh){
                                points(x[i],as.numeric(obs[i]),pch='*',col="red",cex=2)
                        }
                }
        }


legend("topleft",c("All","Pre-2006","Post-2006"),pch=NA,fill=c("black","Lightblue","Darkblue"),pt.cex=2,bty="n",cex=2)
#dev.off()
#}		


type_of_season="data_freq_season2_Gross"
par(mar=c(4,0.5,0.5,0.5))
plot(0,0,t="n",ylim=c(0.15,0.75),xlim=c(0,7),xaxt="n",xlab="",ylab="",yaxt="n")
axis(1,at=c(2,6),lab=rep("",2),cex.axis=2)
#for (type_of_season in c("data_freq_season4_Gross","data_freq_season2_Gross")){
        for (season in 1:2){
                obs=calidris[[type_of_season]][[season]][id_obs]
                ic=lapply(calidris[[type_of_season]][[season]][id_rands],function(x) 2*sd(x[1:100])/sqrt(100))
                p_s=calidris[[type_of_season]][[season]][id_p]
                x=1:3+(season-1)*4
                yHigh=as.numeric(obs)+as.numeric(ic)
                yLow=as.numeric(obs)-as.numeric(ic)
                points(x,obs,t="o",lty=1,pch=21,col=c("Black","Lightblue","Darkblue"),bg=c("Black","Lightblue","Darkblue"),cex=2)
                arrows(x,yHigh,x,yLow,angle=90,length=0.1,code=3,col=c("Black","Lightblue","Darkblue"),lwd=2)
                for (i in 1:3){
                        if (p_s[i]<thresh){
                                points(x[i],as.numeric(obs[i]),pch='*',col="red",cex=2)
                        }
                }
        }

        for (season in 1:2){
                obs=anas[[type_of_season]][[season]][id_obs]
                ic=lapply(anas[[type_of_season]][[season]][id_rands],function(x) 2*sd(x[1:100])/sqrt(100))
                p_s=anas[[type_of_season]][[season]][id_p]
                x=1:3+(season-1)*4
                yHigh=as.numeric(obs)+as.numeric(ic)
                yLow=as.numeric(obs)-as.numeric(ic)
                points(x,obs,t="o",lty=1,pch=22,col=c("Black","Lightblue","Darkblue"),bg=c("Black","Lightblue","Darkblue"),cex=2)
                arrows(x,yHigh,x,yLow,angle=90,length=0.1,code=3,col=c("Black","Lightblue","Darkblue"),lwd=2)
                for (i in 1:3){
                        if (p_s[i]<thresh){
                                points(x[i],as.numeric(obs[i]),pch='*',col="red",cex=2)
                        }
                }
        }
legend("topleft",c("Calidris","Anas"),pch=c(21,22),pt.bg="black",pt.cex=2,bty="n",cex=2)

par(mar=c(4,0.5,0.5,0.5))
plot(0,0,t="n",ylim=c(0.15,0.75),xlim=c(0,7),xaxt="n",xlab="",ylab="",yaxt="n")
axis(1,at=c(2,6),lab=c("Cold","Warm"),cex.axis=2)
#for (type_of_season in c("data_freq_season4_Gross","data_freq_season2_Gross")){
        for (season in 1:2){
                obs=lim[[type_of_season]][[season]][id_obs]
                ic=lapply(lim[[type_of_season]][[season]][id_rands],function(x) 2*sd(x[1:100])/sqrt(100))
                p_s=lim[[type_of_season]][[season]][id_p]
                x=1:3+(season-1)*4
                yHigh=as.numeric(obs)+as.numeric(ic)
                yLow=as.numeric(obs)-as.numeric(ic)
                points(x,obs,t="o",lty=1,pch=23,col=c("Black","Lightblue","Darkblue"),bg=c("Black","Lightblue","Darkblue"),cex=2)
                arrows(x,yHigh,x,yLow,angle=90,length=0.1,code=3,col=c("Black","Lightblue","Darkblue"),lwd=2)
                for (i in 1:3){
                        if (p_s[i]<thresh){
                                points(x[i],as.numeric(obs[i]),pch='*',col="red",cex=2)
                        }
                }
        }

        for (season in 1:2){
                obs=freq[[type_of_season]][[season]][id_obs]
                ic=lapply(freq[[type_of_season]][[season]][id_rands],function(x) 2*sd(x[1:100])/sqrt(100))
                p_s=freq[[type_of_season]][[season]][id_p]
                x=1:3+(season-1)*4
                yHigh=as.numeric(obs)+as.numeric(ic)
                yLow=as.numeric(obs)-as.numeric(ic)
                points(x,obs,t="o",lty=1,pch=24,col=c("Black","Lightblue","Darkblue"),bg=c("Black","Lightblue","Darkblue"),cex=2)
                arrows(x,yHigh,x,yLow,angle=90,length=0.1,code=3,col=c("Black","Lightblue","Darkblue"),lwd=2)
                for (i in 1:3){
                        if (p_s[i]<thresh){
                                points(x[i],as.numeric(obs[i]),pch='*',col="red",cex=2)
                        }
                }
        }

legend("topleft",c("Waders","Freq"),pch=c(23,24),pt.bg="black",pt.cex=2,bty="n",cex=2)
dev.off()
