# Script by CPicoche 2018
#The aim is to compare synchrony indices at different taxonomic levels
library("R.utils")
rm(list=ls())
graphics.off()


freq=loadToEnv("OUT/indices_freq.RData")
lim=loadToEnv("OUT/indices_limicoles.RData")
anas=loadToEnv("OUT/indices_Anas.RData")
calidris=loadToEnv("OUT/indices_Calidris.RData")

pdf("OUT/comparison_taxonomic_level_2seasons.pdf")
par(mfrow=c(2,1),mar=c(2,2,1,1))
which_var="data_freq_season2_Loreau"
data_tmp=freq[[which_var]]
colnames(data_tmp)=c("Cold season","Warm season")
point_bar=barplot(data_tmp,col=c("Black","Lightblue","Darkblue"),border="white",beside=T,legend=rownames(data_tmp),args.legend=list(x='topleft',bty="n"),ylim=c(0,1.1),xlim=c(0,9),main="Loreau")
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
legend("topright",c("Freq","Waders","Anas","Calidris"),pch=c(NA,21,22,24),fill=c("Black",NA,NA,NA),pt.bg="black",border="white",pt.cex=2,bty="n")

which_var="data_freq_season2_Gross"
data_tmp=freq[[which_var]]
point_bar=barplot(data_tmp,col=c("Black","Lightblue","Darkblue"),border="white",beside=T,legend=rownames(data_tmp),args.legend=list(x='topleft',bty="n"),ylim=c(0.,1.1),xlim=c(0,9),main="Gross")
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
