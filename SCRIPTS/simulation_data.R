#CP 2018 
#Simulates community with different size, and different time series length to check the evolution of the Gross synchrony index
#2019/09/12 Added the 2-species case

rm(list=ls())
graphics.off()

#Parameters which change
burn_in=500
#n_species_series=c(5,15,30,50)
n_species_series=c(1)
n_time_series=burn_in+c(35,100,500)
n_repeats=10
b_series=c(0.1,0.5,0.75)

for(n_species in n_species_series){

	for(n_time in n_time_series){
		
		for(b in b_series){
#Parameters which do not change
alpha=0.5
alpha_matrix=diag(2*n_species)+alpha
K=(1+alpha*(n_species*2-1))/(2*n_species)
#b=0.1 #Don't think that's ok
sigma_noise=0.1

for (k.repeats in 1:n_repeats)
{
set.seed(k.repeats)
x<-matrix(0, nrow=n_time,n_species*2,byrow=TRUE) # Matrix of log-abundances or log-densities
epsilon<-matrix(0, nrow=n_time,n_species*2,byrow=TRUE) # Noise on growth rates
x1<-matrix(1/(2*n_species),2*n_species,byrow=TRUE) # Vector of initial abundances - could be changed
x[1,]<-x1

r=rnorm(2*n_species,mean=1,sd=0.25)

#Environmental signal
ymin=0
ymax=1
y0noise=arima.sim(model=list(ar=c(0.1, 0.2, 0.1,0.5,-0.1)), n=n_time,sd=sqrt(0.1) )
y0noise=y0noise/max(abs(y0noise))
y1noise=y0noise
#y1noise<-ymin+(ymax-ymin)*(0.5*y0noise+0.5*(seq(1:n_time)-1)/(n_time-1))

for (t in 1:(n_time-1))
        {
                epsilon[t,]<-rnorm(2*n_species,mean = 0, sd = sigma_noise)# Noise
                x[t+1,1:n_species]<-x[t,1:n_species]*exp(r[1:n_species]*(1-alpha_matrix[1:n_species,1:n_species]%*%x[t,1:n_species]/K)+b*y1noise[t]+epsilon[t,1:n_species])
                x[t+1,(n_species+1):(2*n_species)]<-x[t,(n_species+1):(2*n_species)]*exp(r[(n_species+1):(2*n_species)]*(1-alpha_matrix[(n_species+1):(2*n_species),(n_species+1):(2*n_species)]%*%x[t,(n_species+1):(2*n_species)]/K)-b*y1noise[t]+epsilon[t,(n_species+1):(2*n_species)])
	}
DataPlankton=data.frame(1:(n_time-burn_in),as.numeric(y1noise[(burn_in+1):n_time]),x[(burn_in+1):n_time,])
DataPlankton$Repeat=k.repeats
names(DataPlankton)[1:2]=c("Time_index","Abiotic_var1")

if(k.repeats==1){
	DataPlankton_all=DataPlankton
}else{
	DataPlankton_all=rbind(DataPlankton_all,DataPlankton)
}

} #end k
write.csv(DataPlankton_all,file=paste("OUT/MockBirdsTimeSeries_species",n_species,"_time",(n_time-burn_in),"_b",b,"_notrend_withburnin.csv",sep=""))
} #end b_series
} #end time series
} #end n_species
