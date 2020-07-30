# Script by CPicoche

rm(list=ls())
graphics.off()
source("../test_synchrony_Gross.r")
library('mvcwt')
source("wmr_boot.r")
source("extract_boot.r")
source("../image_mvcwt_for_pvalues.r")
source("../image_mvcwt_for_colormaps.r")
library('RColorBrewer')
source("../iaaft.R")

set.seed(42)
log_b=F
normalize=F
end_bio="abundances"

db=read.csv(paste("IN/summed_",end_bio,"_v2_wtoutrarespecies.csv",sep=""),sep=";",header=T)
db$Date=as.Date(db$Date)

db=subset(db,Nom_latin %in% c("Cormorant","HeronEgret"))
db_tmp=db
names(db_tmp)=c("sp_data_frame","dates","abundance")
db_av=subset(db_tmp,dates<as.Date("2007-01-01"))
db_ap=subset(db_tmp,dates>=as.Date("2007-01-01"))
db_tmp$dates=as.numeric(db_tmp$dates)
db_av$dates=as.numeric(db_av$dates)
db_ap$dates=as.numeric(db_ap$dates)

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

x=(dates-dates[1])/365.25
year_min=1981
#This function computes the Morlet wavelet transform for each bird species separately
mm=mvcwt(x,tab,min.scale=0.2,max.scale=10.0,nscales=212,loc=regularize(x,nsteps=length(x)/2))

print(paste("Start 100",Sys.time()))
distrib100_func=wmr_boot(mm, smoothing = 1,reps=100,output="function")
distrib100_quant=wmr_boot(mm, smoothing = 1,reps=100,output="quantile")
print(paste("Start 1000",Sys.time()))
distrib1000_func=wmr_boot(mm, smoothing = 1,reps=1000,output="function")
distrib1000_quant=wmr_boot(mm, smoothing = 1,reps=1000,output="quantile")
print(paste("Stop 1000",Sys.time()))

save(distrib100_func,distrib100_quant,distrib1000_func,distrib1000_quant,file="distribecdf.RData")


pdf("compare_distrib.pdf",width=10,height=10)
par(mfrow=c(3,3))
id=1:length(distrib100_func$z.boot)
dim(id)=c(length(distrib100_func$x),length(distrib100_func$y))
n_list=sample(1:length(distrib100_func$z.boot),9)
for(i in 1:9){
	plot(distrib100_func$z.boot[[n_list[i]]],col="blue",xlim=c(0,1))
	lines(distrib1000_func$z.boot[[n_list[i]]],col="red",xlim=c(0,1))
	id_obs=which(id==n_list[[i]],arr.ind=T)
	abline(v=distrib100_func$z.boot[[n_list[[i]]]](distrib100_func$z[id_obs[1],id_obs[2],1]),lty=2,col="blue")
	abline(v=distrib1000_func$z.boot[[n_list[[i]]]](distrib1000_func$z[id_obs[1],id_obs[2],1]),lty=2,col="red")
}
dev.off()

pdf("pvalue_for100rands.pdf")
layout(matrix(c(1,2,3,4,5,6),nrow=3,ncol=2,byrow=T),widths=c(10,2))
par(mar=c(3,5,5,3))
image_mvcwt_for_pvalues(distrib100_quant,reset.par=F,cex.axis=4,z.fun="Mod",adj="None",amain="100, Not corrected")
par(mar=c(3,5,2,3))
image_mvcwt_for_pvalues(distrib100_quant,reset.par=F,cex.axis=4,z.fun="Mod",adj="BH",amain="100, BH-corrected")
par(mar=c(3,5,2,3))
image_mvcwt_for_pvalues(distrib100_quant,reset.par=F,cex.axis=4,z.fun="Mod",adj="BY",amain="100, BY-corrected")
dev.off()

pdf("pvalue_for1000rands.pdf")
layout(matrix(c(1,2,3,4,5,6),nrow=3,ncol=2,byrow=T),widths=c(10,2))
par(mar=c(3,5,2,3))
image_mvcwt_for_pvalues(distrib1000_quant,reset.par=F,cex.axis=4,z.fun="Mod",adj="None",amain="1000, Not corrected")
par(mar=c(3,5,2,3))
image_mvcwt_for_pvalues(distrib1000_quant,reset.par=F,cex.axis=4,z.fun="Mod",adj="BH",amain="1000, BH corrected")
par(mar=c(3,5,2,3))
image_mvcwt_for_pvalues(distrib1000_quant,reset.par=F,cex.axis=4,z.fun="Mod",adj="BY",amain="1000, BY corrected")
dev.off()

pdf("CompareBHBY_index_100rands.pdf")
layout(matrix(c(1,2,3,4,5,6),nrow=3,ncol=2,byrow=T),widths=c(10,2))
par(mar=c(3,5,2,3))
image_mvcwt_for_colormaps(distrib100_quant,reset.par=F,cex.axis=4,z.fun="Mod",adj="None",amain="100, Not corrected")
par(mar=c(3,5,2,3))
image_mvcwt_for_colormaps(distrib100_quant,reset.par=F,cex.axis=4,z.fun="Mod",adj="BH",amain="100, BH corrected")
par(mar=c(3,5,2,3))
image_mvcwt_for_colormaps(distrib100_quant,reset.par=F,cex.axis=4,z.fun="Mod",adj="BY",amain="100, BY corrected")
dev.off()

pdf("CompareBHBY_index_1000rands.pdf")
layout(matrix(c(1,2,3,4,5,6),nrow=3,ncol=2,byrow=T),widths=c(10,2))
par(mar=c(3,5,2,3))
image_mvcwt_for_colormaps(distrib1000_quant,reset.par=F,cex.axis=4,z.fun="Mod",adj="None",amain="1000, Not corrected")
par(mar=c(3,5,2,3))
image_mvcwt_for_colormaps(distrib1000_quant,reset.par=F,cex.axis=4,z.fun="Mod",adj="BH",amain="1000, BH corrected")
par(mar=c(3,5,2,3))
image_mvcwt_for_colormaps(distrib1000_quant,reset.par=F,cex.axis=4,z.fun="Mod",adj="BY",amain="1000, BY corrected")
dev.off()

##There are places where we find 0's in the pvalues.
#Find them
pval=1-abs(1-2*distrib100_quant$z.boot)
zb = p.adjust(as.vector(pval), method = "BH") #Adjust all p-values
#Zb<0.1 if and only if pval==0, leading to zb=0
#sum(which(pval==0)==which(zb==0))
#sum(zb==0)
#sum(zb<=0.1)
id_0=which(zb==0)
#To what values of the index do they correspond
id_index0=distrib100_quant$z[id_0]
#Are they THAT specific
#sample_val=sample(1:length(id_0),100)
sample_val=1:length(id_0)
pdf("values_for_pval0_100rands.pdf",height=15,width=15)
par(mfrow=c(5,5),mar=c(0.5,1,4,1))
for(VAL in sample_val){
id_0_array=which(id==id_0[VAL],arr.ind=T)
plot(distrib100_func$z.boot[[id_0[VAL]]],xlim=c(0,1),main=paste("obs",format(distrib100_func$z[id_0_array[1],id_0_array[2],1],digits=3),"// 5% ecdf",format(quantile(distrib100_func$z.boot[[id_0[VAL]]],0.05),digits=3),"// 95% ecdf",format(quantile(distrib100_func$z.boot[[id_0[VAL]]],0.95),digits=3)))
abline(v=distrib100_func$z.boot[[id_0[[VAL]]]](distrib100_func$z[id_0_array[1],id_0_array[2],1]),lty=2)
}
dev.off()

pdf("diff_for_pval0_100rands.pdf")
par(mfrow=c(1,1))
diff=rep(NA,length(id_0))
for(VAL in 18:length(id_0)){
	id_0_array=which(id==id_0[VAL],arr.ind=T)
	valz=distrib100_func$z[id_0_array[1],id_0_array[2],1]
	valmin=quantile(distrib100_func$z.boot[[id_0[VAL]]],0.05)
	valmax=quantile(distrib100_func$z.boot[[id_0[VAL]]],0.95)
	if(valmin>valz){
		diff[VAL]=valmin-valz
	}else if(valmax<valz){
		diff[VAL]=valmax-valz
	}
}
hist(diff,freq=T)
dev.off()
stop()
id_not0=which(zb>0)
count=0
diff_not0=rep(NA,length(id_not0))
pval_val=c()
counter=0
for(VAL in 1:length(id_not0)){
        id_not0_array=which(id==id_not0[VAL],arr.ind=T)
        valz=distrib100_func$z[id_not0_array[1],id_not0_array[2],1]
        valmin=quantile(distrib100_func$z.boot[[id_not0[VAL]]],0.05)
        valmax=quantile(distrib100_func$z.boot[[id_not0[VAL]]],0.95)
        if((valmin>valz)|(valmax<valz)){
		counter=counter+1
		pval_val=c(zb_val,zb[id_not0[VAL]])
	        if(valmin>valz){
        	        diff_not0[VAL]=valmin-valz
        	}else if(valmax<valz){
                	diff_not0[VAL]=valmax-valz
        	}
		if(counter<=4){
		plot(distrib100_func$z.boot[[id_not0[VAL]]],xlim=c(0,1),main=paste("obs",format(distrib100_func$z[id_not0_array[1],id_not0_array[2],1],digits=3),"// 5% ecdf",format(quantile(distrib100_func$z.boot[[id_not0[VAL]]],0.05),digits=3),"// 95% ecdf",format(quantile(distrib100_func$z.boot[[id_not0[VAL]]],0.95),digits=3)))
		abline(v=distrib100_func$z.boot[[id_not0[[VAL]]]](distrib100_func$z[id_not0_array[1],id_not0_array[2],1]),lty=2)

		}
#        	print("There is a problem with your reasoning")
	}
}
#hist(diff_not0,freq=T)
#dev.off()

par(mfrow=c(1,2))
hist(zb[zb>0],freq=F,main="All pval>0")
hist(zb_val,freq=F,main="pval for z outside IC")

pdf(paste("image_for_surrogate.pdf",sep=""),height=12)
layout(matrix(c(1,2,3,4,5,6),nrow=3,ncol=2,byrow=T),widths=c(10,2))
par(mar=c(3,5,2,3))
image_mvcwt_for_colormaps(wmr(mm),amain="reference")
for(i in 1:5){
	xx=extract_boot(mm)
	xx$z=xx$z.boot
	xx$z.boot=NULL
	par(mar=c(3,5,2,3))
	image_mvcwt_for_colormaps(xx)
}
dev.off()

#What happens with IAAFT
pdf(paste("image_for_surrogate_IAAFT.pdf",sep=""),height=12)
layout(matrix(c(1,2,3,4,5,6),nrow=3,ncol=2,byrow=T),widths=c(10,2))
par(mar=c(3,5,2,3))
image_mvcwt_for_colormaps(wmr(mm),amain="reference")
for(i in 1:5){
	tab_tmp=tab
        tab_tmp[,"Cormorant"]=iaaft_surrogate(tab[,'Cormorant'])
        tab_tmp[,"HeronEgret"]=iaaft_surrogate(tab[,'HeronEgret'])
	mmtmp=mvcwt(x,tab_tmp,min.scale=0.2,max.scale=10.0,nscales=212,loc=regularize(x,nsteps=length(x)/2))
        par(mar=c(3,5,2,3))
        image_mvcwt_for_colormaps(wmr(mmtmp))
}
dev.off()


