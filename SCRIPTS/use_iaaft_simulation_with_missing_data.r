#CP 30/07/2020 This script apply wavelet anaysis with IAAFT surrogates to simulation in order to check the effect of missing data

rm(list=ls())
graphics.off()
source("SCRIPTS/image_mvcwt_for_colormaps.r") #Better graphics for image.mvcwt
source("SCRIPTS/iaaft.R") #Better graphics for image.mvcwt
library("RColorBrewer")
library("mvcwt")

#n_species_series=c(5,15,30,50) #Number of species in the community is actually multiplied by 2 in the simulation of the data.
n_species_series=c(5) #As the community which appears in the Appendix is made of 10 species, we focus on this case
n_time_series=c(500)
#b_series=c(0.1,0.5,0.75)
b_series=c(0.1) #As the signal which appears in the Appendix corresponds to a small impact of the driver, we focus on this case
nb_repeats=1
anrands=1000

doyouload=T #True if analyses are taken from previous files, False if we want to relaunch an analysis

j=0
pdf(paste("OUT/worstcase_scenario_wavelets_with",anrands,"_IAAFT.pdf",sep=""),height=15,width=12)
par(mfrow=c(2,1))
for(n_time in n_time_series){
        j=j+1
        i=0
        for(n_species in n_species_series){
                l=0
                i=i+1
                for(b in b_series){
                l=l+1
                print(paste(n_time,n_species,b))
                tab=read.table(paste("../Teich_resultsLFS/simulated_timeseries_Ricker_forGrossIndex/MockBirdsTimeSeries_species",n_species,"_time",n_time,"_b",b,"_withtrend_withburnin.csv",sep=""),header=T,sep=",",dec=".")
                for(k in nb_repeats){
                        data=subset(tab,Repeat==k)
			dates=data$Time_index/12 #There is one point per month, we don't need to regularize

			#Case where there is no missing data
			print("Begin no missing data")
			#data[,4:13] -> Hard-coding is bad practice, but it does use only the species (data[,1] is the row number, data[,2] is the time index,data[,3] is the abiotic variable and data[,14] is index of the repeat
			abundances=data[,4:13]
			print("WARNING: you only use the first 10 species here, check that it's actually what you want to do")
			if(n_species*2!=ncol(abundances)){
				stop("Check again the number of species you are selecting in the table")
			}

			#First, compute the wavelet decomposition for each simulated time series
			mm=mvcwt(dates,abundances,min.scale=0.2,max.scale=10.0,nscales=100) #we reduce the number of scales as compared to the previous version of the article

			if(!doyouload){
			#Then, compute the corresponding wavelet modulus ratio and store it as the reference value
			ref_wmr = wmr(mm)
			ref_val=ref_wmr$z[,,1]
			
			#Now, begin computing the IAAFT surrogates
			tab_values_iaaft=array(NA,dim=c(length(mm$x),length(mm$y),anrands+1))
			tab_values_iaaft[,,anrands+1]=ref_val
			prog.bar = txtProgressBar(min = 0, max = anrands,style = 3)
			for(i in 1:anrands){
			        setTxtProgressBar(prog.bar, i)
			        tab_tmp=abundances
			        for(s in 1:(2*n_species)){
				        tab_tmp[,s]=iaaft_surrogate(abundances[,s])
				}
				#Compute the wavelet decomposition for each surrogate dataset
				mmtmp=mvcwt(dates,tab_tmp,min.scale=0.2,max.scale=10.0,nscales=100)
			        wmr_tmp=wmr(mmtmp)
		        	tab_values_iaaft[,,i]=wmr_tmp$z[,,1]
        		}

			#We compute Pr(X<=x)
			tab_pval=array(NA,dim=c(length(mm$x),length(mm$y),1))
			for(i in 1:length(mm$x)){
			        for(j in 1:length(mm$y)){
			                tab_pval[i,j,1]= sum(tab_values_iaaft[i,j,] <= ref_val[i,j])/(anrands+1)
			                if(tab_pval[i,j,1]>1){stop()}

		        	}
			}
			ref_wmr$z.boot=tab_pval

			#We store all Pr(X<=x) as well as scales, locations and the observed values
			if(length(ref_wmr$x)>length(ref_wmr$y)){
			        yy=c(ref_wmr$y,rep(NA,length(ref_wmr$x)-length(ref_wmr$y)))
			        xx=ref_wmr$x
			}else{
			        xx=c(ref_wmr$x,rep(NA,length(ref_wmr$y)-length(ref_wmr$x)))
			        yy=ref_wmr$y
			}

			tab_xy=cbind(xx,yy)
			colnames(tab_xy)=c("x","y")
			write.table(tab_xy,paste("../Teich_resultsLFS/simulated_timeseries_Ricker_forGrossIndex/tab_xy_mr_simulated_data_with",anrands,"_IAAFT.csv",sep=""),sep=";",dec=".",col.names=T,row.names=F)

			tab_z=ref_wmr$z
			write.table(as.matrix(tab_z[,,1]),paste("../Teich_resultsLFS/simulated_timeseries_Ricker_forGrossIndex/tab_z_mr_simulated_data_with",anrands,"_IAAFT.csv",sep=""),sep=";",dec=".",col.names=F,row.names=F)

			tab_z.boot=ref_wmr$z.boot
			write.table(as.matrix(tab_z.boot[,,1]),paste("../Teich_resultsLFS/simulated_timeseries_Ricker_forGrossIndex/tab_zboot_mr_simulated_data_with",anrands,"_IAAFT.csv",sep=""),sep=";",dec=".",col.names=F,row.names=F)


			}else{ #We load the data instead of relaunching the analyses

			#We need to create a wmr object to store the values, so that it is recognized by the function image_mvcwt. We then store the values from previous analyses in it.
			ref_wmr = wmr(mm)

			tmp_xy=read.csv(paste("../Teich_resultsLFS/simulated_timeseries_Ricker_forGrossIndex/tab_xy_mr_simulated_data_with",anrands,"_IAAFT.csv",sep=""),header=T,sep=";",dec=".")
			ref_wmr$x=tmp_xy[!is.na(tmp_xy[,"x"]),"x"]
			ref_wmr$y=tmp_xy[!is.na(tmp_xy[,"y"]),"y"]

			tmp_z=as.matrix(read.csv(paste("../Teich_resultsLFS/simulated_timeseries_Ricker_forGrossIndex/tab_z_mr_simulated_data_with",anrands,"_IAAFT.csv",sep=""),header=F,sep=";",dec="."))
			tmp_array_z=array(0,dim=c(dim(tmp_z),1))
			tmp_array_z[,,1]=tmp_z
			ref_wmr$z=tmp_array_z

			tmp_z.boot=as.matrix(read.csv(paste("../Teich_resultsLFS/simulated_timeseries_Ricker_forGrossIndex/tab_zboot_mr_simulated_data_with",anrands,"_IAAFT.csv",sep=""),header=F,sep=";",dec="."))
			tmp_array_z.boot=array(0,dim=c(dim(tmp_z.boot),1))
			tmp_array_z.boot[,,1]=tmp_z.boot
			ref_wmr$z.boot=tmp_array_z.boot
			}

			#Begin plot
			layout(matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=T),widths=c(10,2))
			par(mar=c(6,6,3,3))
			image_mvcwt_for_colormaps(ref_wmr,reset.par=F,cex.axis=4,z.fun="Mod",adj="None",amain="")
			mtext("a)",side=2,line=-2,at=0.98,cex=1.5,outer=T,las=1)
			mtext("Scale (years)", side=2, line=-2,at=0.775, outer = TRUE,cex=1.5)

			#Case where there are missing data
			print("Begin with missing data")
			abundance_missingdata=abundances
			for(ss in 1:(n_species*2)){
				abundance_missingdata[sample(1:500,50),ss]=rep(0,50)
			}

			mm_missingdata=mvcwt(dates,abundance_missingdata,min.scale=0.2,max.scale=10.0,nscales=100)
                        
			if(!doyouload){
                        #Then, compute the corresponding wavelet modulus ratio and store it as the reference value
                        ref_wmr_missingdata = wmr(mm_missingdata)
                        ref_val_missingdata=ref_wmr_missingdata$z[,,1]

                        #Now, begin computing the IAAFT surrogates
                        tab_values_iaaft_missingdata=array(NA,dim=c(length(mm_missingdata$x),length(mm_missingdata$y),anrands+1))
                        tab_values_iaaft_missingdata[,,anrands+1]=ref_val_missingdata
                        prog.bar = txtProgressBar(min = 0, max = anrands,style = 3)
                        for(i in 1:anrands){
                                setTxtProgressBar(prog.bar, i)
                                tab_tmp=abundance_missingdata
                                for(s in 1:(2*n_species)){
                                        tab_tmp[,s]=iaaft_surrogate(abundance_missingdata[,s])
                                }
                        #Compute the wavelet decomposition for each surrogate dataset
                        mmtmp=mvcwt(dates,tab_tmp,min.scale=0.2,max.scale=10.0,nscales=100)
                        wmr_tmp=wmr(mmtmp)
                        tab_values_iaaft_missingdata[,,i]=wmr_tmp$z[,,1]
                        }

                        #We compute Pr(X<=x)
                        tab_pval_missingdata=array(NA,dim=c(length(mm_missingdata$x),length(mm_missingdata$y),1))
                        for(i in 1:length(mm_missingdata$x)){
                                for(j in 1:length(mm_missingdata$y)){
                                        tab_pval_missingdata[i,j,1]= sum(tab_values_iaaft_missingdata[i,j,] <= ref_val_missingdata[i,j])/(anrands+1)
                                        if(tab_pval_missingdata[i,j,1]>1){stop()}

                                }
                        }
                        ref_wmr_missingdata$z.boot=tab_pval_missingdata

                        #We store all Pr(X<=x) as well as scales, locations and the observed values
                        if(length(ref_wmr_missingdata$x)>length(ref_wmr_missingdata$y)){
                                yy=c(ref_wmr_missingdata$y,rep(NA,length(ref_wmr_missingdata$x)-length(ref_wmr_missingdata$y)))
                                xx=ref_wmr_missingdata$x
                        }else{
                                xx=c(ref_wmr_missingdata$x,rep(NA,length(ref_wmr_missingdata$y)-length(ref_wmr_missingdata$x)))
                                yy=ref_wmr_missingdata$y
                        }
                        tab_xy=cbind(xx,yy)
                        colnames(tab_xy)=c("x","y")
                        write.table(tab_xy,paste("../Teich_resultsLFS/simulated_timeseries_Ricker_forGrossIndex/tab_xy_mr_simulated_data_missing_with",anrands,"_IAAFT.csv",sep=""),sep=";",dec=".",col.names=T,row.names=F)

                        tab_z=ref_wmr_missingdata$z
                        write.table(as.matrix(tab_z[,,1]),paste("../Teich_resultsLFS/simulated_timeseries_Ricker_forGrossIndex/tab_z_mr_simulated_data_missing_with",anrands,"_IAAFT.csv",sep=""),sep=";",dec=".",col.names=F,row.names=F)

                        tab_z.boot=ref_wmr_missingdata$z.boot
                        write.table(as.matrix(tab_z.boot[,,1]),paste("../Teich_resultsLFS/simulated_timeseries_Ricker_forGrossIndex/tab_zboot_mr_simulated_data_missing_with",anrands,"_IAAFT.csv",sep=""),sep=";",dec=".",col.names=F,row.names=F)


                        }else{ #We load the data instead of relaunching the analyses

                        #We need to create a wmr object to store the values, so that it is recognized by the function image_mvcwt. We then store the values from previous analyses in it.
                        ref_wmr_missingdata = wmr(mm_missingdata)

                        tmp_xy=read.csv(paste("../Teich_resultsLFS/simulated_timeseries_Ricker_forGrossIndex/tab_xy_mr_simulated_data_missing_with",anrands,"_IAAFT.csv",sep=""),header=T,sep=";",dec=".")
                        ref_wmr_missingdata$x=tmp_xy[!is.na(tmp_xy[,"x"]),"x"]
                        ref_wmr_missingdata$y=tmp_xy[!is.na(tmp_xy[,"y"]),"y"]

                        tmp_z=as.matrix(read.csv(paste("../Teich_resultsLFS/simulated_timeseries_Ricker_forGrossIndex/tab_z_mr_simulated_data_missing_with",anrands,"_IAAFT.csv",sep=""),header=F,sep=";",dec="."))
                        tmp_array_z=array(0,dim=c(dim(tmp_z),1))
                        tmp_array_z[,,1]=tmp_z
                        ref_wmr_missingdata$z=tmp_array_z

                        tmp_z.boot=as.matrix(read.csv(paste("../Teich_resultsLFS/simulated_timeseries_Ricker_forGrossIndex/tab_zboot_mr_simulated_data_missing_with",anrands,"_IAAFT.csv",sep=""),header=F,sep=";",dec="."))
                        tmp_array_z.boot=array(0,dim=c(dim(tmp_z.boot),1))
                        tmp_array_z.boot[,,1]=tmp_z.boot
                        ref_wmr_missingdata$z.boot=tmp_array_z.boot
                        }

			par(mar=c(6,6,3,3))
			image_mvcwt_for_colormaps(ref_wmr_missingdata,reset.par=F,cex.axis=4,z.fun="Mod",adj="None",amain="")

			mtext("b)",side=2,line=-2,at=0.49,cex=1.5,outer=T,las=1)
			mtext("Scale (years)", side=2, line=-2,at=0.28, outer = TRUE,cex=1.5)
			mtext("Years",side=1, line=-2, at=0.425,outer = TRUE,cex=1.5)
			

                }
        }
}
}
dev.off()


