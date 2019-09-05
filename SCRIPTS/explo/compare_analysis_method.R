# Script by CPicoche 2018
#The aim is to compare synchrony indices with different definitions of the abundance over a season: mean value or max value
library("R.utils")
rm(list=ls())
graphics.off()

mean_value=loadToEnv("OUT/indices_freq.RData")
max_value=loadToEnv("OUT/indices_freq_using_max_values.RData")

pdf("OUT/comparison_meanmax_level_4seasons.pdf",width=16)
par(mfrow=c(2,1),mar=c(2,2,1,1))
which_var="data_freq_season4_Loreau"
data_tmp=mean_value[[which_var]]
#colnames(data_tmp)=c("Cold season","Warm season")
point_bar=barplot(data_tmp,col=c("Black","Lightblue","Darkblue"),border="white",beside=T,legend=rownames(data_tmp),args.legend=list(x='topleft',bty="n"),ylim=c(0,1.1),xlim=c(0,17),main="Loreau")
data_bis=max_value[[which_var]]
points(c(point_bar),c(data_bis),pch=21,col="white",bg=c("Black","Lightblue","Darkblue"),cex=2)
for(i in 1:floor(length(point_bar)/3)){
lines(c(point_bar)[(1+(i-1)*3):(i*3)],c(data_bis)[(1+(i-1)*3):(i*3)])
}
legend("topright",c("Mean value","Max value"),pch=c(NA,21),fill=c("Black",NA),pt.bg="black",border="white",pt.cex=2,bty="n")

which_var="data_freq_season4_Gross"
data_tmp=mean_value[[which_var]]
point_bar=barplot(data_tmp,col=c("Black","Lightblue","Darkblue"),border="white",beside=T,legend=rownames(data_tmp),args.legend=list(x='topleft',bty="n"),ylim=c(0,1.1),xlim=c(0,17),main="Gross")
data_bis=max_value[[which_var]]
points(c(point_bar),c(data_bis),pch=21,col="white",bg=c("Black","Lightblue","Darkblue"),cex=2)
for(i in 1:floor(length(point_bar)/3)){
lines(c(point_bar)[(1+(i-1)*3):(i*3)],c(data_bis)[(1+(i-1)*3):(i*3)])
}
legend("topright",c("Mean value","Max value"),pch=c(NA,21),fill=c("Black",NA),pt.bg="black",border="white",pt.cex=2,bty="n")
dev.off()
