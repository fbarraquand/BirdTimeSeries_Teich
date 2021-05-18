rm(list=ls())
graphics.off()

library("WaveletComp")


biomass=F
if(biomass){
end_bio="biomasses"
}else{
end_bio="abundances"
}

log_b=F
if(log_b){
        end_log="log"
}else{
        end_log="NOlog"
}
normalize=F


#Load data
db=read.csv(paste("../../IN/summed_",end_bio,"_v2_wtoutrarespecies.csv",sep=""),sep=";",header=T)
db$Date=as.Date(db$Date)

db=subset(db,Nom_latin %in% c("Cormorant","HeronEgret"))

#Build a table with the right format to be analysed
dates=unique(db$Date)
tab=matrix(0,nrow=length(dates),ncol=2)
colnames(tab)=c("Cormorant","HeronEgret")
rownames(tab)=dates
for(id in 1:length(dates)){
        for (s in c('Cormorant',"HeronEgret")){
                id_d=which(db$Date==dates[id]&db$Nom_latin==s)
                if(length(id_d)>0){
                        if(log_b){
                                tab[id,s]=log(db$Nombre[id_d]+1)
                        }else{
                                tab[id,s]=db$Nombre[id_d]
                        }
                }
        }
}

if(normalize){
        tab[,"Cormorant"]=scale(tab[,'Cormorant'])
        tab[,"HeronEgret"]=scale(tab[,'HeronEgret'])
}

x_dates=(dates-dates[1])/365.25

year_min=1981
#This function computes the Morlet wavelet transform for each bird species separately

my.data <- data.frame(date=dates,Cormorant = tab[,"Cormorant"], HeronEgret = tab[,"HeronEgret"])

my.wc <- analyze.coherency(my.data, my.pair = c("Cormorant","HeronEgret"),loess.span = 0,dt = 1/12, dj = 1/100,make.pval = TRUE, n.sim = 1000,date.format="%d-%m-%Y",upperPeriod=6.0)
pdf("phase_difference_cormorant_heron_6elike.pdf",width=10,height=5)
phase_difference_period12_sig=wc.sel.phases(my.wc, sel.period = 1,only.sig = F,which.sig = "wc",siglvl = 0.1,phaselim = c(-pi,+pi*1.25),legend.coords = "topright", legend.horiz = TRUE,main = "", sub = "", timelab = "",show.date=T)
dev.off()
phase_difference_period12_notsig=wc.sel.phases(my.wc, sel.period = 1,only.sig = T,which.sig = "wc",siglvl = 0.1,phaselim = c(-pi,+pi),legend.coords = "topright", legend.horiz = FALSE,main = "", sub = "", timelab = "")
pdf("phase_difference_cormorant_heron_egret.pdf",width=10,height=5)
par(mfrow=c(1,2),oma=c(0.25,0.25,0.25,1.75))
plot(dates,phase_difference_period12_notsig$Angle,ylab="Annual phase difference",t="l",yaxt="n",ylim=c(-1*pi,pi),xlab="")
axis(2,at=seq(-1*pi,pi,by=pi/2),labels=c(expression(-pi),expression(-pi/2),"0",expression(pi/2),expression(pi)))
points(dates,phase_difference_period12_sig$Angle,col="red",pch=16)
#wc.sel.phases(my.wc, sel.period = 12,only.sig = F,which.sig = "wc",siglvl = 0.05,phaselim = c(-pi,+pi),legend.coords = "topright", legend.horiz = FALSE,main = "", sub = "", timelab = "")
wc.phasediff.image(my.wc, which.contour = "wc", use.sAngle = TRUE,n.levels = 250, siglvl = 0.1,legend.params = list(lab = "phase difference levels",lab.line = 3),timelab = "",show.date=T)
dev.off()
