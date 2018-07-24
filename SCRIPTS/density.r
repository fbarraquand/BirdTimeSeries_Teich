#CP 2018
#This script checks the evolution of the density of the birds community with time

rm(list=ls())
graphics.off()

pdf("OUT/density.pdf")
par(mfrow=c(3,3))
print("Anas")
db_warm=read.csv("IN/warmseason_anas_detailed.txt",sep=",",header=T)
db_cold=read.csv("IN/coldseason_anas_detailed.txt",sep=",",header=T)


tab_sum_warm=aggregate(db_warm$abundance_warm,list(db_warm$dates),sum)
tab_sum_cold=aggregate(db_cold$abundance_cold,list(db_cold$dates),sum)

ylimite=log10(1+c(min(c(tab_sum_warm$x,tab_sum_cold$x),na.rm=T),max(c(tab_sum_warm$x,tab_sum_cold$x),na.rm=T)))

plot(tab_sum_warm$Group.1,log10(tab_sum_warm$x+1),ylab="Log10(Anas)",col="red",t="l",lty=1,ylim=ylimite,xlab="")
abline(h=mean(log10(tab_sum_warm$x+1),na.rm=T),col="red",lty=3)
lines(tab_sum_cold$Group.1,log10(tab_sum_cold$x+1),col="blue")
abline(h=mean(log10(tab_sum_cold$x+1),na.rm=T),col="blue",lty=3)

print("Calidris")
db_warm=read.csv("IN/warmseason_calidris_detailed.txt",sep=",",header=T)
db_cold=read.csv("IN/coldseason_calidris_detailed.txt",sep=",",header=T)


tab_sum_warm=aggregate(db_warm$abundance_warm,list(db_warm$dates),sum)
tab_sum_cold=aggregate(db_cold$abundance_cold,list(db_cold$dates),sum)

ylimite=log10(1+c(min(c(tab_sum_warm$x,tab_sum_cold$x),na.rm=T),max(c(tab_sum_warm$x,tab_sum_cold$x),na.rm=T)))

plot(tab_sum_warm$Group.1,log10(tab_sum_warm$x+1),ylab="Log10(Calidris)",col="red",t="l",lty=1,ylim=ylimite,xlab="")
abline(h=mean(log10(tab_sum_warm$x+1),na.rm=T),col="red",lty=3)
lines(tab_sum_cold$Group.1,log10(tab_sum_cold$x+1),col="blue")
abline(h=mean(log10(tab_sum_cold$x+1),na.rm=T),col="blue",lty=3)

print("Anas+Calidris")
db_warm=read.csv("IN/warmseason_abundances_asdataframe_summed.csv",sep=";",header=T)
db_cold=read.csv("IN/coldseason_abundances_asdataframe_summed.csv",sep=";",header=T)

db_warm=subset(db_warm,Species %in% c("Anas","Calidris"))
db_cold=subset(db_cold,Species %in% c("Anas","Calidris"))

tab_sum_warm=aggregate(db_warm$Abundance,list(db_warm$Date),sum)
tab_sum_cold=aggregate(db_cold$Abundance,list(db_cold$Date),sum)

ylimite=log10(1+c(min(c(tab_sum_warm$x,tab_sum_cold$x),na.rm=T),max(c(tab_sum_warm$x,tab_sum_cold$x),na.rm=T)))

plot(tab_sum_warm$Group.1,log10(tab_sum_warm$x+1),ylab="Log10(Anas+Calidris)",col="red",t="l",lty=1,ylim=ylimite,xlab="")
abline(h=mean(log10(tab_sum_warm$x+1),na.rm=T),col="red",lty=3)
lines(tab_sum_cold$Group.1,log10(tab_sum_cold$x+1),col="blue")
abline(h=mean(log10(tab_sum_cold$x+1),na.rm=T),col="blue",lty=3)


print("Waders")
db_warm=read.csv("IN/warmseason_waders_detailed.txt",sep=",",header=T)
db_cold=read.csv("IN/coldseason_waders_detailed.txt",sep=",",header=T)


tab_sum_warm=aggregate(db_warm$abundance_warm,list(db_warm$dates),sum)
tab_sum_cold=aggregate(db_cold$abundance_cold,list(db_cold$dates),sum)

ylimite=log10(1+c(min(c(tab_sum_warm$x,tab_sum_cold$x),na.rm=T),max(c(tab_sum_warm$x,tab_sum_cold$x),na.rm=T)))

plot(tab_sum_warm$Group.1,log10(tab_sum_warm$x+1),ylab="Log10(Waders)",col="red",t="l",lty=1,ylim=ylimite,xlab="")
abline(h=mean(log10(tab_sum_warm$x+1),na.rm=T),col="red",lty=3)
lines(tab_sum_cold$Group.1,log10(tab_sum_cold$x+1),col="blue")
abline(h=mean(log10(tab_sum_cold$x+1),na.rm=T),col="blue",lty=3)

print("Ducks")
db_warm=read.csv("IN/warmseason_duck_detailed.txt",sep=",",header=T)
db_cold=read.csv("IN/coldseason_duck_detailed.txt",sep=",",header=T)


tab_sum_warm=aggregate(db_warm$abundance_warm,list(db_warm$dates),sum)
tab_sum_cold=aggregate(db_cold$abundance_cold,list(db_cold$dates),sum)

ylimite=log10(1+c(min(c(tab_sum_warm$x,tab_sum_cold$x),na.rm=T),max(c(tab_sum_warm$x,tab_sum_cold$x),na.rm=T)))

plot(tab_sum_warm$Group.1,log10(tab_sum_warm$x+1),ylab="Log10(Ducks)",col="red",t="l",lty=1,ylim=ylimite,xlab="")
abline(h=mean(log10(tab_sum_warm$x+1),na.rm=T),col="red",lty=3)
lines(tab_sum_cold$Group.1,log10(tab_sum_cold$x+1),col="blue")
abline(h=mean(log10(tab_sum_cold$x+1),na.rm=T),col="blue",lty=3)

print("Waders+Ducks")
db_warm=read.csv("IN/warmseason_abundances_asdataframe_summed.csv",sep=";",header=T)
db_cold=read.csv("IN/coldseason_abundances_asdataframe_summed.csv",sep=";",header=T)

db_warm=subset(db_warm,Species %in% c("Waders","Ducks"))
db_cold=subset(db_cold,Species %in% c("Waders","Ducks"))

tab_sum_warm=aggregate(db_warm$Abundance,list(db_warm$Date),sum)
tab_sum_cold=aggregate(db_cold$Abundance,list(db_cold$Date),sum)

ylimite=log10(1+c(min(c(tab_sum_warm$x,tab_sum_cold$x),na.rm=T),max(c(tab_sum_warm$x,tab_sum_cold$x),na.rm=T)))

plot(tab_sum_warm$Group.1,log10(tab_sum_warm$x+1),ylab="Log10(Waders+Ducks)",col="red",t="l",lty=1,ylim=ylimite,xlab="")
abline(h=mean(log10(tab_sum_warm$x+1),na.rm=T),col="red",lty=3)
lines(tab_sum_cold$Group.1,log10(tab_sum_cold$x+1),col="blue")
abline(h=mean(log10(tab_sum_cold$x+1),na.rm=T),col="blue",lty=3)


print("Freq")
db_warm=read.csv("IN/warmseason_freq_detailed.txt",sep=",",header=T)
db_cold=read.csv("IN/coldseason_freq_detailed.txt",sep=",",header=T)


tab_sum_warm=aggregate(db_warm$abundance_warm,list(db_warm$dates),sum)
tab_sum_cold=aggregate(db_cold$abundance_cold,list(db_cold$dates),sum)

ylimite=log10(1+c(min(c(tab_sum_warm$x,tab_sum_cold$x),na.rm=T),max(c(tab_sum_warm$x,tab_sum_cold$x),na.rm=T)))

plot(tab_sum_warm$Group.1,log10(tab_sum_warm$x+1),ylab="Log10(Freq)",col="red",t="l",lty=1,ylim=ylimite,xlab="")
abline(h=mean(log10(tab_sum_warm$x+1),na.rm=T),col="red",lty=3)
lines(tab_sum_cold$Group.1,log10(tab_sum_cold$x+1),col="blue")
abline(h=mean(log10(tab_sum_cold$x+1),na.rm=T),col="blue",lty=3)
dev.off()
