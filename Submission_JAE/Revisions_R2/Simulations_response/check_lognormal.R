#### 02/07/2020 This script checks that average distribution of birds have an approximately log-normal distribution


rm(list=ls())
graphics.off()

set.seed(190293)

DBt<-read.csv(file="../../../IN/DBWithMonthlyPhotoTeich_completed.csv",header=TRUE,sep=",",dec=".")

DBt = subset(DBt,(((DBt$Protocol==1 | DBt$Protocol==2) & DBt$Lieu_dit=="USN00-Réserve ornithologique (générique)") | ((DBt$Protocol==1 | DBt$Protocol==2)  & DBt$Lieu_dit=="USN01-Artigues-Réserve ornithologique")))

DBt$Date=as.Date(as.character(DBt$Date))
minAnnee = as.numeric(format(min(DBt$Date), format = "%Y"))
maxAnnee = as.numeric(format(max(DBt$Date), format = "%Y"))

DBt<-DBt[!colnames(DBt) %in% c("X","Ref")] #edit du 14/04/2017 : there are 3 duplicates
DBt = unique.matrix(DBt) # delete with this line,  but I must check the monthly picture
# duplicates are present in the original file, it's not a code pb
# to delete with unique.matrix I must delete the first two columns.

##################################################################################
SpeciesL = as.character(unique(DBt$Nom_latin)) #-> 280
SpeciesF = as.character(unique(DBt$Nom_espece)) #same thing but with french names.

#minimum of observation to keep a species = Frequent Species
vec_n_obs_oiseaux_t=rep(0,length(SpeciesL))
for (i in 1:length(SpeciesL)){
  vec_n_obs_oiseaux_t[i]=sum(as.character(DBt$Nom_latin)==SpeciesL[i],na.rm = TRUE)
}
oiseaux_Frequents_t=SpeciesL[vec_n_obs_oiseaux_t>75] #60/279
oiseaux_Frequents_t_F = SpeciesF[vec_n_obs_oiseaux_t>75] 

DB_frequent=subset(DBt,DBt$Nom_latin %in% oiseaux_Frequents_t)

x=aggregate(DB_frequent$Nombre,list(DB_frequent$Nom_latin),mean)

pdf("distrib_average_abundance.pdf",height=4)
par(mfrow=c(1,2))
hist(x$x,xlab="",ylab="",main="Raw abundance",breaks=9)
hist(log(x$x),xlab="",ylab="",main="Log(abundance)",ylim=c(0,22))
m=mean(log(x$x))
log_mu=rnorm(60,m)
hist(log_mu,xlab="",ylab="",main="Simulated log(ab)",add=T,col=rgb(0,0,1,alpha=0.2),breaks=9)
dev.off()

