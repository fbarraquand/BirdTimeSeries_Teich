### Deprecated

rm(list=ls())
graphics.off()
library("R.utils")


anas=loadToEnv("OUT/data_Anas.RData")
calidris=loadToEnv("OUT/data_Calidris.RData")
duck=loadToEnv("OUT/data_freq.RData")
wader=loadToEnv("OUT/data_limicoles.RData")

val=c("wintering_all","wintering_pre_2006","wintering_post_2006","breeding_all","breeding_pre_2006","breeding_post_2006")

#Il y a 44 ans dans tous les cas
yy=1973:2016


anas_all_wintering=aggregate(anas$wintering_all[,3],list(anas$wintering_all$dates),sum)
anas_all_breeding=aggregate(anas$breeding_all[,3],list(anas$breeding_all$dates),sum)

calidris_all_wintering=aggregate(calidris$wintering_all[,3],list(calidris$wintering_all$dates),sum)
calidris_all_breeding=aggregate(calidris$breeding_all[,3],list(calidris$breeding_all$dates),sum)

duck_all_wintering=aggregate(duck$wintering_all[,3],list(duck$wintering_all$dates),sum)
duck_all_breeding=aggregate(duck$breeding_all[,3],list(duck$breeding_all$dates),sum)

wader_all_wintering=aggregate(wader$wintering_all[,3],list(wader$wintering_all$dates),sum)
wader_all_breeding=aggregate(wader$breeding_all[,3],list(wader$breeding_all$dates),sum)

mini=log10(1+min(c(anas_all_wintering$x,anas_all_breeding$x,calidris_all_wintering$x,calidris_all_breeding$x)))
maxi=log10(1+max(c(anas_all_wintering$x,anas_all_breeding$x,calidris_all_wintering$x,calidris_all_breeding$x)))

pdf("OUT/timeseries_season.pdf",width=12.5,height=7)
plot(yy,log10(anas_all_wintering$x+1),pch=16,col="blue",ylim=c(mini,maxi),t="o",ylab="Log10(abundance)",xlab="Year",cex=2)
points(yy,log10(calidris_all_wintering$x+1),pch=16,col="darkblue",t="o",cex=2)
points(yy+0.25,log10(anas_all_breeding$x+1),pch=16,col="orange",t="o",cex=2)
points(yy+0.25,log10(calidris_all_breeding$x+1),pch=16,col="red",t="o",cex=2)
abline(v=2007,lty=2)
legend("topleft",c("Anas cold","Calidris cold",'Anas warm','Calidris warm'),col=c('blue','darkblue','orange','red'),pch=16)

mini=log10(1+min(c(wader_all_wintering$x,wader_all_breeding$x,duck_all_wintering$x,duck_all_breeding$x)))
maxi=log10(1+max(c(wader_all_wintering$x,wader_all_breeding$x,duck_all_wintering$x,duck_all_breeding$x)))
plot(yy,log10(wader_all_wintering$x+1),pch=16,col="blue",ylim=c(mini,maxi),t="o",ylab="Log10(abundance)",xlab="Year",cex=2)
points(yy,log10(duck_all_wintering$x+1),pch=16,col="darkblue",t="o",cex=2)
points(yy+0.25,log10(wader_all_breeding$x+1),pch=16,col="orange",t="o",cex=2)
points(yy+0.25,log10(duck_all_breeding$x+1),pch=16,col="red",t="o",cex=2)
abline(v=2007,lty=2)
legend("topleft",c("Wader cold","Duck cold",'Wader warm','Duck warm'),col=c('blue','darkblue','orange','red'),pch=16)
dev.off()
