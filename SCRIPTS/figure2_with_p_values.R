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

pdf("tmp_figure2.pdf",width=20)
plot(0,0,t="n",ylim=c(0,0.8),xlim=c(0,15),xaxt="n",xlab="",ylab="")
#for (type_of_season in c("data_freq_season4_Gross","data_freq_season2_Gross")){
type_of_season="data_freq_season4_Gross"
axis(1,at=c(2,6,10,14),lab=c("Winter","Spring","Summer","Autumn"))
	for (season in 1:4){
		obs=calidris[[type_of_season]][[season]][id_obs]
		ic=lapply(calidris[[type_of_season]][[season]][id_rands],function(x) 2*sd(x[1:100])/sqrt(100))
		p_s=calidris[[type_of_season]][[season]][id_p]
		x=1:3+(season-1)*4
		yHigh=as.numeric(obs)+as.numeric(ic)
		yLow=as.numeric(obs)-as.numeric(ic)
		points(x,obs,t="o",lty=1,pch=21,col="black",bg=c("Black","Lightblue","Darkblue"),cex=2)
		arrows(x,yHigh,x,yLow,angle=90,length=0.1,code=3,col=c("Black","Lightblue","Darkblue"),lwd=2)
		for (i in 1:3){
			if (p_s[i]<thresh){
				points(x[i],as.numeric(obs[i])+0.05,pch='*',col="black",cex=2)
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
                points(x,obs,t="o",lty=1,pch=22,col="black",bg=c("Black","Lightblue","Darkblue"),cex=2)
                arrows(x,yHigh,x,yLow,angle=90,length=0.1,code=3,col=c("Black","Lightblue","Darkblue"),lwd=2)
                for (i in 1:3){
                        if (p_s[i]<thresh){
                                points(x[i],as.numeric(obs[i])+0.05,pch='*',col="black",cex=2)
                        }
                }
        }

        for (season in 1:4){
                obs=lim[[type_of_season]][[season]][id_obs]
                ic=lapply(lim[[type_of_season]][[season]][id_rands],function(x) 2*sd(x[1:100])/sqrt(100))
                p_s=lim[[type_of_season]][[season]][id_p]
                x=1:3+(season-1)*4
                yHigh=as.numeric(obs)+as.numeric(ic)
                yLow=as.numeric(obs)-as.numeric(ic)
                points(x,obs,pch=23,t="o",lty=1,col="black",bg=c("Black","Lightblue","Darkblue"),cex=2)
                arrows(x,yHigh,x,yLow,angle=90,length=0.1,code=3,col=c("Black","Lightblue","Darkblue"),lwd=2)
                for (i in 1:3){
                        if (p_s[i]<thresh){
                                points(x[i],as.numeric(obs[i])+0.05,pch='*',col="black",cex=2)
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
                points(x,obs,pch=24,t="o",lty=1,col="black",bg=c("Black","Lightblue","Darkblue"),cex=2)
                arrows(x,yHigh,x,yLow,angle=90,length=0.1,code=3,col=c("Black","Lightblue","Darkblue"),lwd=2)
                for (i in 1:3){
                        if (p_s[i]<thresh){
                                points(x[i],as.numeric(obs[i])+0.05,pch='*',col="black",cex=2)
                        }
                }
        }


legend("topright",c("Anas","Calidris","Waders","Freq"),pch=c(21,22,23,24),pt.bg="black",border="white",pt.cex=2,bty="n")
dev.off()
stop()
#}		


data_bis=anas[[which_var]]
points(c(point_bar),c(data_bis),pch=22,col="white",bg=c("Black","Lightblue","Darkblue"),cex=2)
for(i in 1:floor(length(point_bar)/3)){
lines(c(point_bar)[(1+(i-1)*3):(i*3)],c(data_bis)[(1+(i-1)*3):(i*3)])
}


stop()

pdf("OUT/comparison_taxonomic_level_2seasons.pdf",width=12)
par(mfrow=c(1,1),mar=c(2,2,1,1))

#which_var="data_freq_season2_Loreau"
#data_tmp=freq[[which_var]]
#colnames(data_tmp)=c("Cold season","Warm season")
#point_bar=barplot(data_tmp,col=c("Black","Lightblue","Darkblue"),border="white",beside=T,legend=rownames(data_tmp),args.legend=list(x='topleft',bty="n"),ylim=c(0,1.1),xlim=c(0,9),main="Loreau")
#data_bis=lim[[which_var]]
#points(c(point_bar),c(data_bis),pch=21,col="white",bg=c("Black","Lightblue","Darkblue"),cex=2)
#for(i in 1:floor(length(point_bar)/3)){
#lines(c(point_bar)[(1+(i-1)*3):(i*3)],c(data_bis)[(1+(i-1)*3):(i*3)])
#}
#data_bis=anas[[which_var]]
#points(c(point_bar),c(data_bis),pch=22,col="white",bg=c("Black","Lightblue","Darkblue"),cex=2)
#for(i in 1:floor(length(point_bar)/3)){
#lines(c(point_bar)[(1+(i-1)*3):(i*3)],c(data_bis)[(1+(i-1)*3):(i*3)])
#}
#data_bis=calidris[[which_var]]
#points(c(point_bar),c(data_bis),pch=24,col="white",bg=c("Black","Lightblue","Darkblue"),cex=2)
#for(i in 1:floor(length(point_bar)/3)){
#lines(c(point_bar)[(1+(i-1)*3):(i*3)],c(data_bis)[(1+(i-1)*3):(i*3)])
#}
#legend("topright",c("Freq","Waders","Anas","Calidris"),pch=c(NA,21,22,24),fill=c("Black",NA,NA,NA),pt.bg="black",border="white",pt.cex=2,bty="n")

which_var="data_freq_season2_Gross"
data_tmp=freq[[which_var]]
point_bar=barplot(data_tmp,col=c("Black","Lightblue","Darkblue"),border="white",beside=T,legend=rownames(data_tmp),args.legend=list(x='topleft',bty="n"),ylim=c(0.,0.7),xlim=c(0,9),main="")
data_bis=lim[[which_var]]
points(c(point_bar),c(data_bis),pch=21,col="white",bg=c("Black","Lightblue","Darkblue"),cex=2)
for(i in 1:floor(length(point_bar)/3)){
lines(c(point_bar)[(1+(i-1)*3):(i*3)],c(data_bis)[(1+(i-1)*3):(i*3)])
}
data_bis=anas[[which_var]]
points(c(point_bar),c(data_bis),pch=22,col="white",bg=c("Black","Lightblue","Darkblue"),cex=2)
for(i in 1:floor(length(point_bar)/3)){
lines(c(point_bar)[(1+(i-1)*3):(i*3)],c(data_bis)[(1+(i-1)*3):(i*3)])
}
data_bis=calidris[[which_var]]
points(c(point_bar),c(data_bis),pch=24,col="white",bg=c("Black","Lightblue","Darkblue"),cex=2)
for(i in 1:floor(length(point_bar)/3)){
lines(c(point_bar)[(1+(i-1)*3):(i*3)],c(data_bis)[(1+(i-1)*3):(i*3)])
}
legend("topright",c("Freq","Waders","Anas","Calidris"),pch=c(NA,21,22,24),pt.bg="black",fill=c("black",NA,NA,NA),border="white",pt.cex=2,bty="n")
dev.off()

pdf("OUT/comparison_taxonomic_level_4seasons.pdf",width=14)
par(mfrow=c(1,1),mar=c(2,2,1,1))
which_var="data_freq_season4_Gross"
data_tmp=freq[[which_var]]
point_bar=barplot(data_tmp,col=c("Black","Lightblue","Darkblue"),border="white",beside=T,legend=rownames(data_tmp),args.legend=list(x='topleft',bty="n"),ylim=c(0.,0.7),xlim=c(0,17),main="")
data_bis=lim[[which_var]]
points(c(point_bar),c(data_bis),pch=21,col="white",bg=c("Black","Lightblue","Darkblue"),cex=2)
for(i in 1:floor(length(point_bar)/3)){
lines(c(point_bar)[(1+(i-1)*3):(i*3)],c(data_bis)[(1+(i-1)*3):(i*3)])
}
data_bis=anas[[which_var]]
points(c(point_bar),c(data_bis),pch=22,col="white",bg=c("Black","Lightblue","Darkblue"),cex=2)
for(i in 1:floor(length(point_bar)/3)){
lines(c(point_bar)[(1+(i-1)*3):(i*3)],c(data_bis)[(1+(i-1)*3):(i*3)])
}
data_bis=calidris[[which_var]]
points(c(point_bar),c(data_bis),pch=24,col="white",bg=c("Black","Lightblue","Darkblue"),cex=2)
for(i in 1:floor(length(point_bar)/3)){
lines(c(point_bar)[(1+(i-1)*3):(i*3)],c(data_bis)[(1+(i-1)*3):(i*3)])
}
legend("topright",c("Freq","Waders","Anas","Calidris"),pch=c(NA,21,22,24),pt.bg="black",fill=c("black",NA,NA,NA),border="white",pt.cex=2,bty="n")
dev.off()

