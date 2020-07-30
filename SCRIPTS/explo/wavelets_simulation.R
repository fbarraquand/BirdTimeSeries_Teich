rm(list=ls())
graphics.off()
source("SCRIPTS/image_mvcwt.r") #Add to change the image function to have a nice Color Bar

n_species_series=c(5,15,30,50)
n_time_series=c(500)
b_series=c(0.1,0.5,0.75)
nb_repeats=1
#par(mfrow=c(3,1))
j=0
pdf("OUT/worstcase_scenario_wavelets.pdf")
for(n_time in n_time_series){
        j=j+1
        i=0
        for(n_species in 5:5){
                l=0
                i=i+1
                for(b in c(0.1)){
                l=l+1
                print(paste(n_time,n_species,b))
                tab=read.table(paste("OUT/MockBirdsTimeSeries_species",n_species,"_time",n_time,"_b",b,"_withtrend_withburnin.csv",sep=""),header=T,sep=",",dec=".")
                for(k in nb_repeats){
                        data=subset(tab,Repeat==k)
			dates=data$Time_index/12
			mm=mvcwt(dates,data[,4:13],min.scale=0.2,max.scale=10.0)

			mr = wmr.boot(mm, smoothing = 1,reps=100)
			image_mvcwt(mr,reset.par=F,cex.axis=4,z.fun="Mod")

			abundance_tmp=data[,4:13]	
			for(ss in 1:(n_species*2)){
				abundance_tmp[sample(1:500,50),ss]=rep(0,50)
			}
			mmbis=mvcwt(dates,abundance_tmp,min.scale=0.2,max.scale=10.0)
			mrbis = wmr.boot(mmbis, smoothing = 1,reps=100)
			image_mvcwt(mrbis,reset.par=F,cex.axis=4,z.fun="Mod")
			

                }
        }
}
}
dev.off()


