#### ATTENTION 
#FAUT METTRE A JOUR "DIRECTORY_ORIGIN" POUR UTILISER CE SCRIPT
#analyse des données après la création de la photo mensuelle.
rm(list=ls())
graphics.off()
library("lattice")
library ("RColorBrewer") #pour la génération automatique de couleur
library("corrplot")
library("lubridate")
library("stringr")
library("codyn")
library ("pracma")
DIRECTORY_ORIGIN = "/home/caluome/git/BirdTimeSeries_Teich/"
#DIRECTORY_ORIGIN = "/home/caluome/Documents/LabExCOTE/DATA/Graphes/Axe1/TRANSFERT_LIMICOLES"
setwd(paste(DIRECTORY_ORIGIN,"IN/",sep=""))
options(nwarnings = 300) #nb de messages de warnings conservés

#update.packages(ask = F) #update all packages

##################################################################################
# -- importation des données du Teich
# DBt<-read.csv(file="/home/caluome/Documents/DATA/DATA/le_Teich/DBWithMonthlyPhotoTeich.csv",header=TRUE,sep=",",dec=".")
# DBt<-read.csv(file="DBWithMonthlyPhotoTeich.csv",header=TRUE,sep=",",dec=".")
DBt<-read.csv(file="DBWithMonthlyPhotoTeich_completed.csv",header=TRUE,sep=",",dec=".")
#DBt<-read.csv(file="DBWithMonthlyPhotoTeich_maxsum.csv",header=TRUE,sep=",",dec=".")

#DBt = subset(DBt,(DBt$Protocol==1 | DBt$Protocol==2) & DBt$Lieu_dit=="USN00-Réserve ornithologique (générique)")
#DBt = subset(DBt,((DBt$Protocol==1 | DBt$Protocol==2) & (DBt$Lieu_dit=="USN00-Réserve ornithologique (générique)") | (DBt$Annee==1987 & DBt$Lieu_dit=="USN01-Artigues-Réserve ornithologique")))
#DBt = subset(DBt,((DBt$Protocol==1 | DBt$Protocol==2) & (DBt$Lieu_dit=="USN00-Réserve ornithologique (générique)") |  DBt$Lieu_dit=="USN01-Artigues-Réserve ornithologique"))
#DBt = subset(DBt,((DBt$Protocol==1 | DBt$Protocol==2) & (DBt$Lieu_dit=="USN00-Réserve ornithologique (générique)") | (DBt$Annee==1987 & DBt$Lieu_dit=="USN01-Artigues-Réserve ornithologique")| (DBt$Annee==2003 & DBt$Lieu_dit=="USN01-Artigues-Réserve ornithologique")))
#DBt = subset(DBt,((DBt$Protocol==1 | DBt$Protocol==2) & (DBt$Lieu_dit=="USN00-Réserve ornithologique (générique)") |  DBt$Lieu_dit=="USN01-Artigues-Réserve ornithologique"))

DBt = subset(DBt,(((DBt$Protocol==1 | DBt$Protocol==2) & DBt$Lieu_dit=="USN00-Réserve ornithologique (générique)") | ((DBt$Protocol==1 | DBt$Protocol==2)  & DBt$Lieu_dit=="USN01-Artigues-Réserve ornithologique")))

DBt$Date=as.Date(as.character(DBt$Date))
minAnnee = as.numeric(format(min(DBt$Date), format = "%Y"))
maxAnnee = as.numeric(format(max(DBt$Date), format = "%Y"))

DBt<-DBt[!colnames(DBt) %in% c("X","Ref")] #edit du 14/04/2017 : there are 3 duplicates
DBt = unique.matrix(DBt) # delete with this line,  but I must check the monthly picture
# duplicates are present in the original file, it's not a code pb
# to delete with unique.matrix I must delete the first two columns.

##################################################################################
#How many species ? 
SpeciesL = as.character(unique(DBt$Nom_latin)) #-> 280
SpeciesF = as.character(unique(DBt$Nom_espece)) #same thing but with french names.

#minimum of observation to keep a species = Frequent Species
vec_n_obs_oiseaux_t=rep(0,length(SpeciesL)) 
for (i in 1:length(SpeciesL)){
  vec_n_obs_oiseaux_t[i]=sum(as.character(DBt$Nom_latin)==SpeciesL[i],na.rm = TRUE)
}
oiseaux_Frequents_t=SpeciesL[vec_n_obs_oiseaux_t>20] #125/279
oiseaux_Frequents_t=SpeciesL[vec_n_obs_oiseaux_t>50] #75/279
oiseaux_Frequents_t=SpeciesL[vec_n_obs_oiseaux_t>75] #60/279
oiseaux_Frequents_t_F = SpeciesF[vec_n_obs_oiseaux_t>75] #60/279 #same thing but with french names.



# wader birds
limicoles = c("Recurvirostra avosetta","Limosa limosa","Limosa lapponica","Calidris temminckii","Calidris canutus",
              "Calidris alba","Calidris alpina","Calidris minuta","Calidris maritima" ,"Gallinago gallinago",
              "Tringa flavipes","Tringa nebularia","Tringa erythropus","Tringa ochropus","Tringa totanus",
              "Tringa glareola","Actitis hypoleucos","Philomachus pugnax","Numenius arquata","Numenius phaeopus",
              "Himantopus himantopus","Charadrius hiaticula","Charadrius alexandrinus","Haematopus ostralegus",
              "Burhinus oedicnemus","Charadrius dubius","Phalaropus lobatus","Pluvialis squatarola",
              "Pluvialis apricaria","Arenaria interpres","Vanellus vanellus")


birds_to_remove=c("Anas discors","Anas americana","Calidris melanotos","Calidris pusilla","Calidris ruficollis","Calidris fuscicollis","Calidris himantopus","Burhinus oedicnemus","Phalaropus lobatus","Charadrius alexandrinus","Haematopus ostralegus")

DBt=subset(DBt,!DBt$Nom_latin %in% birds_to_remove)
summed_abundances<-read.csv(file="summed_abundances.csv",header=TRUE,sep=";",dec=".") 
calidris_summed_abundances<-read.csv(file="coldseason_calidris_detailed.txt",header=TRUE,sep=",",dec=".") 
setwd(paste(DIRECTORY_ORIGIN,"OUT/",sep=""))

# ---------------------------------------------------------------------
#
#                              SCRIPT T-36
#
# ---------------------------------------------------------------------
#### T36 : SYNCHRONY CALCUL by MONTH
SynchronyMonth=function(matrice,file_out,titre,loreau=TRUE, gross=TRUE,ymin=-1,ymax=1){ 
  pdf(file_out,width=14,height=8)
  par(mar=c(4,4.5,5,2.5))
  par(oma = c(4.2, 0.5, 0.5, 0.5))
  vec_synchrony_Loreau = rep(0,12) #final results
  vec_synchrony_Gross  = rep(0,12) #final results
  for (m in 1:12){
    subdata   = subset(matrice,as.numeric(format(matrice$Date, format = "%m"))==m)
    date      = subdata$Date
    date      = c(as.numeric(date))
    abundance = c(as.numeric(subdata$Nombre))
    species   = c(as.character(subdata$Nom_latin))
    dataF_DBt = data.frame(date,abundance,species) #Codyn need dataframe to work
    colnames(dataF_DBt)=c("date","abundance","species") 
    vec_synchrony_Loreau[m]= synchrony(dataF_DBt, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Loreau")
    vec_synchrony_Gross[m]=synchrony(dataF_DBt, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Gross") 
    
  }
  if(loreau==TRUE & gross==TRUE){
  plot(1:12,vec_synchrony_Loreau,type="o",ylim=c(ymin,ymax),pch=19, 
       main=titre,ylab="Value of synchrony",xlab="month",cex.main=2,cex.lab=2,
       cex=2,cex.axis=2)
  lines(1:12,vec_synchrony_Gross,col="orange",type = "o",pch=19,cex=2)
  }
  if(loreau == TRUE & gross==FALSE){
    plot(1:12,vec_synchrony_Loreau,type="o",ylim=c(ymin,ymax),pch=19, 
         main=titre,ylab="Value of synchrony",xlab="month",cex.main=2,cex.lab=2,
         cex=2,cex.axis=2)
  }
  if(loreau==FALSE & gross ==TRUE){
    plot(1:12,vec_synchrony_Gross,type="o",ylim=c(ymin,ymax),pch=19, 
         main=titre,ylab="Value of synchrony",xlab="month",cex.main=2,cex.lab=2,
         cex=2,cex.axis=2,col="orange")
  }
  if(loreau ==FALSE & gross==FALSE){print("ERREUR, l'une des métriques doit être égale à true")
    exit()}
  abline(h=0,col='black',lty=3)
  print ("Loreau toute la série")
  print (vec_synchrony_Loreau)
  print ("Gross toute la série")
  print (vec_synchrony_Gross)
  #-----------------------
  # PERIODE >= 2006
	vec_synchrony_Loreau = rep(0,12) #Final results
  vec_synchrony_Gross = rep(0,12) #Final results
  for (m in 1:12){
    #subdata   = subset(DBt,as.numeric(format(DBt$Date, format = "%m"))==m & as.numeric(format(DBt$Date, format = "%Y"))>=2006 & DBt$Nom_latin %in% limicoles)
    subdata   = subset(DBt,as.numeric(format(DBt$Date, format = "%m"))==m & as.numeric(format(DBt$Date, format = "%Y"))>=2007)
    date      = subdata$Date
    date      = c(as.numeric(date))
    abundance = c(as.numeric(subdata$Nombre))
    species   = c(as.character(subdata$Nom_latin))
    dataF_DBt = data.frame(date,abundance,species)#Codyn need dataframe to work
    colnames(dataF_DBt)=c("date","abundance","species") 
    vec_synchrony_Loreau[m]= synchrony(dataF_DBt, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Loreau")
    vec_synchrony_Gross[m]=synchrony(dataF_DBt, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Gross") 
    
  }
  if(loreau==TRUE){lines(1:12,vec_synchrony_Loreau,col="blue",type = "o",pch=19,cex=2)}
  if (gross==TRUE){lines(1:12,vec_synchrony_Gross,col="red",type = "o",pch=19,cex=2)}
  print ("Loreau après 2006")
  print (vec_synchrony_Loreau)
  print ("Gross  après 2006")
  print (vec_synchrony_Gross)
  #-----------------------
  # PERIODE < 2006
  vec_synchrony_Loreau = rep(0,12) #Final results
  vec_synchrony_Gross = rep(0,12) #Final results
  for (m in 1:12){
    subdata   = subset(DBt,as.numeric(format(DBt$Date, format = "%m"))==m & as.numeric(format(DBt$Date, format = "%Y"))<2007 )
    date      = subdata$Date
    date      = c(as.numeric(date))
    abundance = c(as.numeric(subdata$Nombre))
    species   = c(as.character(subdata$Nom_latin))
    dataF_DBt = data.frame(date,abundance,species) #Codyn need dataframe to work
    colnames(dataF_DBt)=c("date","abundance","species") 
    vec_synchrony_Loreau[m]= synchrony(dataF_DBt, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Loreau")
    vec_synchrony_Gross[m]=synchrony(dataF_DBt, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Gross") 
    
  }
  if(loreau==TRUE){lines(1:12,vec_synchrony_Loreau,col="lightblue",type = "o",pch=19,cex=2)}
  if(gross==TRUE){lines(1:12,vec_synchrony_Gross,col="lightgreen",type = "o",pch=19,cex=2)}
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  if(loreau == TRUE & gross==TRUE){
    legend("bottom",c("Synchony by Loreau","Synchrony by Gross","Synchony by Loreau from 2006","Synchrony by Gross from 2006","Synchony by Loreau before 2006","Synchrony by Gross before 2006"), col=c("black","orange","blue","red","lightblue","lightgreen"),pch=c(16,16), xpd = TRUE, 
           horiz=TRUE,
           inset = c(0,0),bty = "n",cex=0.7,pt.cex=1.5,lty=1,pt.lwd=1.5,lwd=1.5)
  }
  if(loreau==TRUE & gross==FALSE){
    legend("bottom",c("Synchony by Loreau","Synchony by Loreau from 2006","Synchony by Loreau before 2006"), col=c("black","blue","lightblue"),pch=c(16,16), xpd = TRUE,  horiz=TRUE,inset = c(0,0),bty = "n",cex=0.7,pt.cex=1.5,lty=1,pt.lwd=1.5,lwd=1.5)
  }
  
  if(loreau==FALSE & gross==TRUE){
    legend("bottom",c("Synchrony by Gross","Synchrony by Gross from 2006","Synchrony by Gross before 2006"), col=c("orange","red","lightgreen"),pch=c(16,16), xpd = TRUE,  horiz=TRUE,inset = c(0,0),bty = "n",cex=0.7,pt.cex=1.5,lty=1,pt.lwd=1.5,lwd=1.5)
    
  }
  dev.off()
  print ("Loreau Avant 2006")
  print (vec_synchrony_Loreau)
  print ("Gross  Avant 2006")
  print (vec_synchrony_Gross)
}

## - Calls of the function n°T36
#SynchronyMonth(DBt,"t36-Axe1-Teich-SynchronyComparison_byMonth_AllSpecies.pdf","Evolution of Synchrony from Loreau and Gross by month (all species)")
tempDBt = subset(DBt,as.character(DBt$Nom_latin)%in% oiseaux_Frequents_t )
SynchronyMonth(tempDBt,"t38-Axe1-Teich-SynchronyComparison_byMonth_frequentSpecies.pdf","Evolution of Synchrony from Loreau and Gross by month (frequent species)")
tempDBt = subset(DBt,as.character(DBt$Nom_latin)%in% limicoles )
SynchronyMonth(tempDBt,"t39-Axe1-Teich-SynchronyComparison_byMonth_onlywadingBirds_test.pdf","Evolution of Synchrony from Loreau and Gross by month (Wading Birds)")


# ---------------------------------------------------------------------
#
#                              SCRIPT T-37
#
# ---------------------------------------------------------------------
#### T37 : SYNCHRONY BY MONTH, ONLY ANAS AND CALIDRIS
Genre_synchronie =c('Anas','Calidris')

pdf("t37-Axe1-Teich-SynchronyComparison_byMonth_selectedGenres.pdf",width=14,height=8)
for (g in 1:length(Genre_synchronie)){
  par(mar=c(4,4.5,3,2.5))
  par(oma = c(4.2, 0.5, 0.5, 0.5))
  vec_synchrony_Loreau = rep(0,12)  #final results
  vec_synchrony_Gross  = rep(0,12)  #final results
  for (m in 1:12){
    #subdata   = subset(DBt,as.numeric(format(DBt$Date, format = "%m"))==m & grepl("^Anas",DBt$Nom_latin,ignore.case = TRUE))
    subdata   = subset(DBt,as.numeric(format(DBt$Date, format = "%m"))==m & grepl(paste("^",Genre_synchronie[g],sep=""),DBt$Nom_latin,ignore.case = TRUE))
    date      = subdata$Date
    date      = c(as.numeric(date))
    abundance = c(as.numeric(subdata$Nombre))
    species   = c(as.character(subdata$Nom_latin))
    dataF_DBt = data.frame(date,abundance,species) #Codyn need a dataframe to work
    colnames(dataF_DBt)=c("date","abundance","species") 
    vec_synchrony_Loreau[m]= synchrony(dataF_DBt, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Loreau")
    vec_synchrony_Gross[m]=synchrony(dataF_DBt, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Gross") 
    
  }
  plot(1:12,vec_synchrony_Loreau,type="o",ylim=c(-1,1),pch=19, 
       main=paste("Evolution of Synchrony from Loreau and Gross by month (",Genre_synchronie[g],")"),ylab="Value of synchrony",xlab="month",cex.main=2,cex.lab=2,
       cex=2,cex.axis=2)
  lines(1:12,vec_synchrony_Gross,col="orange",type = "o",pch=19,cex=2)
  #-----------------------
  # PERIODE >= 2006
  vec_synchrony_Loreau = rep(0,12) #Final results
  vec_synchrony_Gross  = rep(0,12) #Final results
  for (m in 1:12){
    subdata   = subset(DBt,as.numeric(format(DBt$Date, format = "%m"))==m & as.numeric(format(DBt$Date, format = "%Y"))>=2006 & grepl(paste("^",Genre_synchronie[g],sep=""),DBt$Nom_latin,ignore.case = TRUE))
    date      = subdata$Date
    date      = c(as.numeric(date))
    abundance = c(as.numeric(subdata$Nombre))
    species   = c(as.character(subdata$Nom_latin))
    dataF_DBt = data.frame(date,abundance,species) 
    colnames(dataF_DBt)=c("date","abundance","species") 
    vec_synchrony_Loreau[m]= synchrony(dataF_DBt, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Loreau")
    vec_synchrony_Gross[m]=synchrony(dataF_DBt, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Gross") 
    
  }
  lines(1:12,vec_synchrony_Loreau,col="blue",type = "o",pch=19,cex=2)
  lines(1:12,vec_synchrony_Gross,col="red",type = "o",pch=19,cex=2)
  #-----------------------
  # PERIODE < 2006
  vec_synchrony_Loreau = rep(0,12)  #Final results
  vec_synchrony_Gross  = rep(0,12)  #Final results
  for (m in 1:12){
    subdata   = subset(DBt,as.numeric(format(DBt$Date, format = "%m"))==m & as.numeric(format(DBt$Date, format = "%Y"))<2006 & grepl(paste("^",Genre_synchronie[g],sep=""),DBt$Nom_latin,ignore.case = TRUE))
    date      = subdata$Date
    date      = c(as.numeric(date))
    abundance = c(as.numeric(subdata$Nombre))
    species   = c(as.character(subdata$Nom_latin))
    dataF_DBt = data.frame(date,abundance,species) 
    colnames(dataF_DBt)=c("date","abundance","species") 
    vec_synchrony_Loreau[m]= synchrony(dataF_DBt, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Loreau")
    vec_synchrony_Gross[m]=synchrony(dataF_DBt, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Gross") 
  }
  lines(1:12,vec_synchrony_Loreau,col="lightblue",type = "o",pch=19,cex=2)
  lines(1:12,vec_synchrony_Gross,col="lightgreen",type = "o",pch=19,cex=2)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom",c("Synchony by Loreau","Synchrony by Gross","Synchony by Loreau from 2006","Synchrony by Gross from 2006","Synchony by Loreau before 2006","Synchrony by Gross before 2006"), col=c("black","orange","blue","red","lightblue","lightgreen"),pch=c(16,16), xpd = TRUE, 
         horiz=TRUE,
         inset = c(0,0),bty = "n",cex=0.7,pt.cex=1,lty=1,pt.lwd=1.5,lwd=1.5)
}

dev.off()


# ---------------------------------------------------------------------
#
#                              SCRIPT T-40
#
# ---------------------------------------------------------------------
#### T40 : SYNCHRONY BY SEASON 
SynchronySeason = function(matrice,file_out,titre,saison,type){
  #saison=list(c(11,12,1,2,3),c(5,6,7,8,9))
  vec_synchrony_Loreau = rep(0,length(saison)*3) #final results
  vec_synchrony_Gross = rep(0,length(saison)*3) #final results
  for (s in 1:length(saison)){
    subdata   = subset(matrice,as.numeric(format(matrice$Date, format = "%m")) %in% saison[[s]])
    date      = subdata$Date
    date      = c(as.numeric(date))
    abundance = c(as.numeric(subdata$Nombre))
    species   = c(as.character(subdata$Nom_latin))
    dataF_DBt = data.frame(date,abundance,species) #Codyn need dataframe
    colnames(dataF_DBt)=c("date","abundance","species") 
    vec_synchrony_Loreau[s]= synchrony(dataF_DBt, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Loreau")
    vec_synchrony_Gross[s]=synchrony(dataF_DBt, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Gross") 
    
  }
  #-----------------------
  # PERIODE >= 2006
  for (m in 1:length(saison)){
    subdata   = subset(matrice,as.numeric(format(matrice$Date, format = "%m")) %in% saison[[m]] & as.numeric(format(matrice$Date, format = "%Y"))>=2006)
    date      = subdata$Date
    date      = c(as.numeric(date))
    abundance = c(as.numeric(subdata$Nombre))
    species   = c(as.character(subdata$Nom_latin))
    dataF_DBt = data.frame(date,abundance,species) #Codyn need dataframe
    colnames(dataF_DBt)=c("date","abundance","species") 
    vec_synchrony_Loreau[m+length(saison)]= synchrony(dataF_DBt, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Loreau")
    vec_synchrony_Gross[m+length(saison)]=synchrony(dataF_DBt, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Gross") 
    
  }
  #-----------------------
  # PERIODE < 2006
  for (m in 1:length(saison)){
    subdata   = subset(matrice,as.numeric(format(matrice$Date, format = "%m")) %in% saison[[m]] & as.numeric(format(matrice$Date, format = "%Y"))<2006)
    date      = subdata$Date
    date      = c(as.numeric(date))
    abundance = c(as.numeric(subdata$Nombre))
    species   = c(as.character(subdata$Nom_latin))
    dataF_DBt = data.frame(date,abundance,species) #Codyn need dataframe
    colnames(dataF_DBt)=c("date","abundance","species") 
    vec_synchrony_Loreau[m+(length(saison)*2)]= synchrony(dataF_DBt, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Loreau")
    vec_synchrony_Gross[m+(length(saison)*2)]=synchrony(dataF_DBt, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Gross") 
  }
  
  data_synchro =c(c(vec_synchrony_Loreau),c(vec_synchrony_Gross))
  max_value = max(data_synchro,na.rm=TRUE)+0.002
  data_synchro = matrix(data_synchro,nc=length(type), byrow=T)
  #type =c("winter-all","summer-all","winter>=2006","summer>=2006","winter<2006","summer<2006")
  colnames(data_synchro) = type
  print (vec_synchrony_Loreau)
  print (vec_synchrony_Gross)
  pdf(file_out,width=20,height=8)
  par(mar=c(4,4.5,3,2.5))
  par(oma = c(4.2, 0.5, 0.5, 0.5))
  barplot(data_synchro, xlab="X",ylab="Y", beside=T, col=c("#FFFFAA","#AAFFAA"), ylim=c(0,max_value), main=titre,cex.lab=1,cex.axis = 1,cex.main=2) # Tracer les barres
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom",c("Synchony of Loreau","Synchrony of Gross"), col=c("#FFFFAA","#AAFFAA"),pch=c(16,16), xpd = TRUE, 
         horiz=TRUE,
         inset = c(0,0),bty = "n",cex=0.8,pt.cex=1.5,lty=1,pt.lwd=1.5,lwd=1.5)
  dev.off()
}
## - Calls of the function n° T40
saison=list(c(11,12,1,2,3),c(5,6,7,8,9))
type =c("winter-all","summer-all","winter>=2006","summer>=2006","winter<2006","summer<2006")
#SynchronySeason(DBt,"t40-Axe1-Teich-Comparison_synchrony_by_season_all_species.pdf","Comparison Synchrony Season - All species",saison,type)
winter = c(1,2,3) #number of month for each season
spring = c(4,5,6)
summer = c(7,8,9)
autumn = c(10,11,12)
saison2=list(winter,spring,summer,autumn)
type =c("winter-all","spring-all","summer-all","autumn-all","winter>=2006","spring>=2006","summer>=2006","autumn>=2006","winter<2006","spring<2006","summer<2006","autumn<2006")
SynchronySeason(DBt,"t40-Axe1-Teich-Comparison_synchrony_4seasons_all_species.pdf","Comparison of the synchrony between the 4 seasons for all bird species",saison2,type)
# ---------------------------------------------------------------------
#
#                              SCRIPT T-40 - V2
#
# ---------------------------------------------------------------------
###### SYNCHRONy BY SEASON, other way of calculating, using the overlap of winter over two years.
#### T40 V2  
SynchronySeason2 = function(matrice,file_out,titre,Loreau = TRUE, Gross =TRUE,max_value = 1,min_value=-1){
  # saison=list(c(11,12,1,2,3),c(5,6,7,8,9))
  # saison =  list(hiver_n=c(11, 12), hiver_n1=c(1,2,3),ete=c(5,6,7,8,9))
  saison = list(hiver_n=c(11,12),hiver_n1=c(1,2),ete=c(5,6,7,8))
  all_species = unique(as.character(matrice$Nom_latin))
  vec_synchrony_Loreau=rep(0,6) #6 for 2season by 3 periods (all,before 2006 and from 2006)
  vec_synchrony_Gross=rep(0,6)
  abondance_ete=c()
  date_ete=c()
  species_ete=c()
  abondance_hiver=c()
  date_hiver=c()
  species_hiver=c()
  for(y in min(matrice$Annee):(max(matrice$Annee)-1)){
      subdata_ete = subset(matrice,as.numeric(format(matrice$Date, format = "%m")) %in% saison[['ete']] & as.numeric(format(matrice$Date, format = "%Y"))==y)
      subdata_hiver = subset(matrice,(as.numeric(format(matrice$Date, format = "%m")) %in% saison[['hiver_n']] & as.numeric(format(matrice$Date, format = "%Y"))==y) | (as.numeric(format(matrice$Date, format = "%m")) %in% saison[['hiver_n1']] & as.numeric(format(matrice$Date, format = "%Y"))==y+1))
      #print (paste("y",y,dim(subdata_ete)[1]))
      #print (paste("y",y,dim(subdata_hiver)[1]))
      if(dim(subdata_ete)[1]>0){
        species_en_cours = c(unique(as.character(subdata_ete$Nom_latin)))
        for (e in 1:length(species_en_cours)){
          #print (paste("sepcies",species_en_cours[e]))
          #print (paste("mean",mean(subdata_ete$Nombre[subdata_ete$Nom_latin==species_en_cours[e]])))
          #abondance_ete=c(abondance_ete,c(mean(subdata_ete$Nombre[subdata_ete$Nom_latin==species_en_cours[e]])))
          moyenne = sum(subdata_ete$Nombre[subdata_ete$Nom_latin==species_en_cours[e]])/4
          abondance_ete=c(abondance_ete,c(moyenne))
          date_ete=c(date_ete,c(as.Date(as.character(subdata_ete$Date[1]),format="%Y-%m-%d"))) 
          species_ete =c(species_ete,c(species_en_cours[e]))
        }
      }
      if(dim(subdata_ete)[1]==0){
        for (s in 1:length(all_species)){
          date = as.Date(paste(y,"-07-15",sep=""))
          abondance_ete=c(abondance_ete,c(0))
          date_ete=c(date_ete,c(date)) 
          species_ete =c(species_ete,c(all_species[s]))
        }
        
      }
      if(dim(subdata_hiver)[1]>0){
        species_en_cours = c(unique(as.character(subdata_hiver$Nom_latin)))
        for (e in 1:length(species_en_cours)){
          #print (paste("sepcies",species_en_cours[e]))
          #print (paste("mean",mean(subdata_hiver$Nombre[subdata_hiver$Nom_latin==species_en_cours[e]])))
          #abondance_hiver=c(abondance_hiver,c(mean(subdata_hiver$Nombre[subdata_hiver$Nom_latin==species_en_cours[e]])))
          moyenne = sum(subdata_hiver$Nombre[subdata_hiver$Nom_latin==species_en_cours[e]])/4
          abondance_hiver=c(abondance_hiver,c(moyenne))
          date_hiver=c(date_hiver,c(as.Date(as.character(subdata_hiver$Date[1]),format="%Y-%m-%d"))) 
          species_hiver =c(species_hiver,c(species_en_cours[e]))
        }
      }
      if(dim(subdata_hiver)[1]==0){  # I create inputs to 0 artificially when there is 0 data for one year, otherwise it impacts the calculation of synchrony
        for (s in 1:length(all_species)){
          date = as.Date(paste(y,"-12-15",sep=""))
          abondance_hiver=c(abondance_hiver,c(0))
          date_hiver=c(date_hiver,c(date)) 
          species_hiver =c(species_hiver,c(all_species[s]))
        }
        
      }
  }
  date      = c(as.numeric(date_ete))
  abundance = c(as.numeric(abondance_ete))
  species   = c(as.character(species_ete))
  dataF_ete = data.frame(date,abundance,species) #Codyn need dataframe
  colnames(dataF_ete)=c("date","abundance","species") 
  if(length(dataF_ete$abundance[dataF_ete$abundance>0])>0){
    vec_synchrony_Loreau[1]= synchrony(dataF_ete, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Loreau")
    vec_synchrony_Gross[1]=synchrony(dataF_ete, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Gross")
  }
  else{
    vec_synchrony_Loreau[1]=NA
    vec_synchrony_Gross[1]=NA
  }
  date      = c(as.numeric(date_hiver))
  abundance = c(as.numeric(abondance_hiver))
  species   = c(as.character(species_hiver))
  dataF_hiver = data.frame(date,abundance,species) 
  colnames(dataF_hiver)=c("date","abundance","species") 
  if(length(dataF_hiver$abundance[dataF_hiver$abundance>0])>0){
    vec_synchrony_Loreau[2]= synchrony(dataF_hiver, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Loreau")
    vec_synchrony_Gross[2]=synchrony(dataF_hiver, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Gross")
  }
  else{
    vec_synchrony_Loreau[2]=NA
    vec_synchrony_Gross[2]=NA
  }
  #-----------------------
  # PERIODE >= 2007
  abondance_ete=c()
  date_ete=c()
  species_ete=c()
  abondance_hiver=c()
  date_hiver=c()
  species_hiver=c()
  for(y in 2007:(max(matrice$Annee)-1)){
    subdata_ete = subset(matrice,as.numeric(format(matrice$Date, format = "%m")) %in% saison[['ete']] & as.numeric(format(matrice$Date, format = "%Y"))==y)
    subdata_hiver = subset(matrice,(as.numeric(format(matrice$Date, format = "%m")) %in% saison[['hiver_n']] & as.numeric(format(matrice$Date, format = "%Y"))==y) | (as.numeric(format(matrice$Date, format = "%m")) %in% saison[['hiver_n1']] & as.numeric(format(matrice$Date, format = "%Y"))==y+1))
    
    if(dim(subdata_ete)[1]>0){
      species_en_cours = c(unique(as.character(subdata_ete$Nom_latin)))
      for (e in 1:length(species_en_cours)){
        #print (paste("sepcies",species_en_cours[e]))
        #print (paste("mean",mean(subdata_ete$Nombre[subdata_ete$Nom_latin==species_en_cours[e]])))
        #abondance_ete=c(abondance_ete,c(mean(subdata_ete$Nombre[subdata_ete$Nom_latin==species_en_cours[e]])))
        moyenne = sum(subdata_ete$Nombre[subdata_ete$Nom_latin==species_en_cours[e]])/4
        abondance_ete=c(abondance_ete,c(moyenne))
        date_ete=c(date_ete,c(as.Date(as.character(subdata_ete$Date[1]),format="%Y-%m-%d"))) 
        species_ete =c(species_ete,c(species_en_cours[e]))
      }
    }
    if(dim(subdata_ete)[1]==0){
      for (s in 1:length(all_species)){
        date = as.Date(paste(y,"-07-15",sep=""))
        abondance_ete=c(abondance_ete,c(0))
        date_ete=c(date_ete,c(date)) 
        species_ete =c(species_ete,c(all_species[s]))
      }
      
    }
    if(dim(subdata_hiver)[1]>0){
      species_en_cours = c(unique(as.character(subdata_hiver$Nom_latin)))
      for (e in 1:length(species_en_cours)){
        #print (paste("sepcies",species_en_cours[e]))
        #print (paste("mean",mean(subdata_hiver$Nombre[subdata_hiver$Nom_latin==species_en_cours[e]])))
        #abondance_hiver=c(abondance_hiver,c(mean(subdata_hiver$Nombre[subdata_hiver$Nom_latin==species_en_cours[e]])))
        moyenne = sum(subdata_hiver$Nombre[subdata_hiver$Nom_latin==species_en_cours[e]])/4
        abondance_hiver=c(abondance_hiver,c(moyenne))
        date_hiver=c(date_hiver,c(as.Date(as.character(subdata_hiver$Date[1]),format="%Y-%m-%d"))) 
        species_hiver =c(species_hiver,c(species_en_cours[e]))
      }

    }
    if(dim(subdata_hiver)[1]==0){ # I create inputs to 0 artificially when there is 0 data for one year, otherwise it impacts the calculation of synchrony
      for (s in 1:length(all_species)){
        date = as.Date(paste(y,"-12-15",sep=""))
        abondance_hiver=c(abondance_hiver,c(0))
        date_hiver=c(date_hiver,c(date)) 
        species_hiver =c(species_hiver,c(all_species[s]))
      }
    }
  }
  date      = c(as.numeric(date_ete))
  abundance = c(as.numeric(abondance_ete))
  species   = c(as.character(species_ete))
  dataF_ete = data.frame(date,abundance,species) #Codyn need dataframe
  colnames(dataF_ete)=c("date","abundance","species") 
  if(length(dataF_ete$abundance[dataF_ete$abundance>0])>0){
    vec_synchrony_Loreau[3]= synchrony(dataF_ete, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Loreau")
    vec_synchrony_Gross[3]=synchrony(dataF_ete, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Gross")
  }
  else{
    vec_synchrony_Loreau[3]=NA
    vec_synchrony_Gross[3]=NA
  }
  date      = c(as.numeric(date_hiver))
  abundance = c(as.numeric(abondance_hiver))
  species   = c(as.character(species_hiver))
  dataF_hiver = data.frame(date,abundance,species)
  colnames(dataF_hiver)=c("date","abundance","species") 
  if(length(dataF_hiver$abundance[dataF_hiver$abundance>0])>0){
    vec_synchrony_Loreau[4]= synchrony(dataF_hiver, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Loreau")
    vec_synchrony_Gross[4]=synchrony(dataF_hiver, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Gross")
  }
  else{
    vec_synchrony_Loreau[4]=NA
    vec_synchrony_Gross[4]=NA
  }
  #-----------------------
  # PERIODE < 2006
  abondance_ete=c()
  date_ete=c()
  species_ete=c()
  abondance_hiver=c()
  date_hiver=c()
  species_hiver=c()
  for(y in min(matrice$Annee):2006){
    subdata_ete = subset(matrice,as.numeric(format(matrice$Date, format = "%m")) %in% saison[['ete']] & as.numeric(format(matrice$Date, format = "%Y"))==y)
    subdata_hiver = subset(matrice,(as.numeric(format(matrice$Date, format = "%m")) %in% saison[['hiver_n']] & as.numeric(format(matrice$Date, format = "%Y"))==y) | (as.numeric(format(matrice$Date, format = "%m")) %in% saison[['hiver_n1']] & as.numeric(format(matrice$Date, format = "%Y"))==y+1))
    #print (paste("y",y,dim(subdata_hiver)[1]))
    if(dim(subdata_ete)[1]>0){
      species_en_cours = c(unique(as.character(subdata_ete$Nom_latin)))
      for (e in 1:length(species_en_cours)){
        #print (paste("sepcies",species_en_cours[e]))
        #print (paste("mean",mean(subdata_ete$Nombre[subdata_ete$Nom_latin==species_en_cours[e]])))
        #abondance_ete=c(abondance_ete,c(mean(subdata_ete$Nombre[subdata_ete$Nom_latin==species_en_cours[e]])))
        moyenne = sum(subdata_ete$Nombre[subdata_ete$Nom_latin==species_en_cours[e]])/4
        abondance_ete=c(abondance_ete,c(moyenne))
        date_ete=c(date_ete,c(as.Date(as.character(subdata_ete$Date[1]),format="%Y-%m-%d"))) 
        species_ete =c(species_ete,c(species_en_cours[e]))
      }
    }
    if(dim(subdata_ete)[1]==0){
      for (s in 1:length(all_species)){
        date = as.Date(paste(y,"-07-15",sep=""))
        abondance_ete=c(abondance_ete,c(0))
        date_ete=c(date_ete,c(date)) 
        species_ete =c(species_ete,c(all_species[s]))
      }
      
    }
    if(dim(subdata_hiver)[1]>0){
      species_en_cours = c(unique(as.character(subdata_hiver$Nom_latin)))
      for (e in 1:length(species_en_cours)){
        #print (paste("sepcies",species_en_cours[e]))
        #print (paste("mean",mean(subdata_hiver$Nombre[subdata_hiver$Nom_latin==species_en_cours[e]])))
        #abondance_hiver=c(abondance_hiver,c(mean(subdata_hiver$Nombre[subdata_hiver$Nom_latin==species_en_cours[e]])))
        moyenne = sum(subdata_hiver$Nombre[subdata_hiver$Nom_latin==species_en_cours[e]])/4
        abondance_hiver=c(abondance_hiver,c(moyenne))
        date_hiver=c(date_hiver,c(as.Date(as.character(subdata_hiver$Date[1]),format="%Y-%m-%d"))) 
        species_hiver =c(species_hiver,c(species_en_cours[e]))
      }
    }
    if(dim(subdata_hiver)[1]==0){  # I create inputs to 0 artificially when there is 0 data for one year, otherwise it impacts the calculation of synchrony
      for (s in 1:length(all_species)){
        date = as.Date(paste(y,"-12-15",sep=""))
        abondance_hiver=c(abondance_hiver,c(0))
        date_hiver=c(date_hiver,c(date)) 
        species_hiver =c(species_hiver,c(all_species[s]))
      }
    }
  }
  date      = c(as.numeric(date_ete))
  abundance = c(as.numeric(abondance_ete))
  species   = c(as.character(species_ete))
  dataF_ete = data.frame(date,abundance,species) #Codyn need dataframe
  colnames(dataF_ete)=c("date","abundance","species") 
  if(length(dataF_ete$abundance[dataF_ete$abundance>0])>0){
    vec_synchrony_Loreau[5]= synchrony(dataF_ete, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Loreau")
    vec_synchrony_Gross[5]=synchrony(dataF_ete, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Gross")
  }
  else{
    vec_synchrony_Loreau[5]=NA
    vec_synchrony_Gross[5]=NA
  }
  date      = c(as.numeric(date_hiver))
  abundance = c(as.numeric(abondance_hiver))
  species   = c(as.character(species_hiver))
  dataF_hiver = data.frame(date,abundance,species)
  colnames(dataF_hiver)=c("date","abundance","species") 
  
  if(length(dataF_hiver$abundance[dataF_hiver$abundance>0])>0){
    vec_synchrony_Loreau[6]= synchrony(dataF_hiver, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Loreau")
    vec_synchrony_Gross[6]=synchrony(dataF_hiver, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Gross")
  }
  else{
    vec_synchrony_Loreau[6]=NA
    vec_synchrony_Gross[6]=NA
  }
  data_synchro=c()
  
  if (Loreau == TRUE & Gross == TRUE){data_synchro =c(c(vec_synchrony_Loreau),c(vec_synchrony_Gross))
  }
  if(Loreau ==TRUE & Gross ==FALSE){data_synchro =c(vec_synchrony_Loreau)}
  if(Loreau ==FALSE & Gross ==TRUE){data_synchro =c(vec_synchrony_Gross)
  print ("Gross :")
  print ("Sum.-All    Wint.-All   Sum.>=2007    Wint.>=2007   Sum.<2007    wint.<2007")
  print (vec_synchrony_Gross)}
  if(Loreau == FALSE & Gross ==FALSE){exit()}
  #max_value = max(data_synchro,na.rm=TRUE)+0.1
  #min_value = min(data_synchro,na.rm=TRUE)-0.1
  data_synchro = matrix(data_synchro,nc=6, byrow=T)
  type =c("Summer-all","winter-all","summer>=2007","winter>=2007","summer<2007","winter<2007")
  colnames(data_synchro) = type
  pdf(file_out,width=14,height=8)
  par(mar=c(4,4.5,3,2.5))
  par(oma = c(4.2, 0.5, 0.5, 0.5))
  if (Loreau == TRUE & Gross == TRUE){
  barplot(data_synchro, xlab="X",ylab="Y", beside=T, col=c("#FFFFAA","#AAFFAA"), ylim=c(min_value,max_value), main=titre,cex.lab=1,cex.axis = 1,cex.main=1.2)}
  if (Loreau == TRUE & Gross == FALSE){
    barplot(data_synchro, xlab="X",ylab="Y", beside=T, col=c("#FFFFAA"), ylim=c(min_value,max_value), main=titre,cex.lab=1,cex.axis = 1,cex.main=1.2)}
  if (Loreau == FALSE & Gross == TRUE){
    barplot(data_synchro, xlab="X",ylab="Y", beside=T, col=c("#AAFFAA"), ylim=c(min_value,max_value), main=titre,cex.lab=1,cex.axis = 1,cex.main=1.2)
    }
  
  
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  if (Loreau == TRUE & Gross == TRUE){legend("bottom",c("Synchony by Loreau","Synchrony by Gross"), col=c("#FFFFAA","#AAFFAA"),pch=c(16,16), xpd = TRUE, 
         horiz=TRUE,
         inset = c(0,0),bty = "n",cex=0.8,pt.cex=1.5,lty=1,pt.lwd=1.5,lwd=1.5)}
  if (Loreau == TRUE & Gross == FALSE){legend("bottom",c("Synchony by Loreau"), col=c("#FFFFAA"),pch=c(16,16), xpd = TRUE, 
                                             horiz=TRUE,
                                             inset = c(0,0),bty = "n",cex=0.8,pt.cex=1.5,lty=1,pt.lwd=1.5,lwd=1.5)}
  if (Loreau == FALSE & Gross == TRUE){legend("bottom",c("Synchrony by Gross"), col=c("#AAFFAA"),pch=c(16,16), xpd = TRUE, 
                                             horiz=TRUE,
                                             inset = c(0,0),bty = "n",cex=0.8,pt.cex=1.5,lty=1,pt.lwd=1.5,lwd=1.5)}
  
  dev.off()
}
## - call function
SynchronySeason2(DBt,"t40-Axe1-Teich-Comparison_synchrony_season_all_species_v2.pdf","Comparison Synchrony Season - All species")
tempDBt = subset(DBt,as.character(DBt$Nom_latin)%in% oiseaux_Frequents_t )
SynchronySeason2(tempDBt,"t41-Axe1-Teich-Synchrony_season_only_frequent_species.pdf","Comparison Synchrony Season - Frequent Species")

##################################################################
# TEST manual verification of values, to compare with coralie results
saison =  list(hiver_n=c(11, 12), hiver_n1=c(1,2,3),ete=c(5,6,7,8,9))
for(y in min(DBt$Annee):max(DBt$Annee)-1){
  subdata_ete = subset(DBt,as.numeric(format(DBt$Date, format = "%m")) %in% saison[['ete']] & as.numeric(format(DBt$Date, format = "%Y"))==y)
  subdata_hiver = subset(DBt,(as.numeric(format(DBt$Date, format = "%m")) %in% saison[['hiver_n']] & as.numeric(format(DBt$Date, format = "%Y"))==y) | (as.numeric(format(DBt$Date, format = "%m")) %in% saison[['hiver_n1']] & as.numeric(format(DBt$Date, format = "%Y"))==y+1))

  if(dim(subdata_ete)[1]>0){
    #var_ete=by(subdata_ete$Nombre,as.character(subdata_ete$Nom_latin),sd)
    var_ete=by(subdata_ete$Nombre,as.character(subdata_ete$Nom_latin),var)

    print (paste("année",y))
    print (paste("Dimension subdata : ",dim(subdata_ete)[1]))


    tmp_ete_sum=by(subdata_ete$Nombre,subdata_ete$Date,sum)
    var_ete_sum_all=var(tmp_ete_sum)

    print (var_ete_sum_all)
  }
}

##################################################################
# TEST same thing, but on calculation for all years
saison =  list(hiver_n=c(11, 12), hiver_n1=c(1,2,3),ete=c(5,6,7,8,9))
subdata_ete = subset(DBt,as.numeric(format(DBt$Date, format = "%m")) %in% saison[['ete']])
subdata_hiver = subset(DBt,as.numeric(format(DBt$Date, format = "%m")) %in% saison[['hiver_n']]| as.numeric(format(DBt$Date, format = "%m")) %in% saison[['hiver_n1']])
subdata_ete=subset(subdata_ete,subdata_ete$Nombre !=0)      
# variance of communities
# summer
tmp_ete_sum=by(subdata_ete$Nombre,subdata_ete$Date,sum)
var_ete_sum_all=var(tmp_ete_sum)
print(paste("Variance of communities - summer :",var_ete_sum_all))
# winter
tmp_hiver_sum=by(subdata_hiver$Nombre,subdata_hiver$Date,sum)
var_hiver_sum_all=var(tmp_hiver_sum)
print(paste("Variance of communities - winter :",var_hiver_sum_all))

# variance species by species
# summer
var_ete=by(subdata_ete$Nombre,as.character(subdata_ete$Nom_latin),var)
var_ete_com_all = sum(sqrt(var_ete),na.rm=TRUE)
print(paste("Variance species by species - summer :",var_ete_com_all))
#winter
var_hiver=by(subdata_hiver$Nombre,as.character(subdata_hiver$Nom_latin),var)
var_hiver_com_all = sum(sqrt(var_hiver),na.rm=TRUE) # @COCO : So, you have to take the square root before you sum on all species..
print(paste("Variance species by species - winter :",var_hiver_com_all))

# ---------------------------------------------------------------------
#
#                              SCRIPT T-42
#
# ---------------------------------------------------------------------

SynchronySeason_moving_month=function(matrice,file_out,titre){ 
	# idea: Study of Synchronie Comparison over rolling months, to see if a group of months at values really different from the others groups
	# -> Too boring to calculate the seasons, I make rather lists.
	ete =list(c(4,5,6),c(5,6,7),c(6,7,8),c(7,8,9),c(4,5,6,7),c(5,6,7,8),c(6,7,8,9),c(4,5,6,7,8),c(5,6,7,8,9),c(4,5,6,7,8,9))
	hivern = list(c(10,11,12),c(11,12),c(12),c(),c(10,11,12),c(11,12),c(12),c(10,11,12),c(11,12),c(10,11,12))
	hivern1 = list(c(),c(1),c(1,2),c(1,2,3),c(1),c(1,2),c(1,2,3),c(1,2),c(1,2,3),c(1,2,3))

	vec_synchrony_Loreau_ete=c()
	vec_synchrony_Gross_ete=c()
	vec_synchrony_Loreau_hiver=c()
	vec_synchrony_Gross_hiver=c()
	for (l in 1:length(ete)){
	  ### summer
	  abondance_ete=c()
	  date_ete=c()
	  species_ete=c()
	  for(y in min(matrice$Annee):max(matrice$Annee)-1){
		 subdata_ete = subset(matrice,as.numeric(format(matrice$Date, format = "%m")) %in% ete[[l]] & as.numeric(format(matrice$Date, format = "%Y"))==y)
		 if(dim(subdata_ete)[1]>0){
		   species_en_cours = c(unique(as.character(subdata_ete$Nom_latin)))
		   for (e in 1:length(species_en_cours)){
		     if(length(subdata_ete$Nombre[subdata_ete$Nom_latin==species_en_cours[e]])>0){
		       abondance_ete=c(abondance_ete,c(mean(subdata_ete$Nombre[subdata_ete$Nom_latin==species_en_cours[e]])))
		       date_ete=c(date_ete,c(as.Date(as.character(subdata_ete$Date[1]),format="%Y-%m-%d"))) 
		       species_ete =c(species_ete,c(species_en_cours[e]))
		     }
		   }
		 }
	  }
	  # summer
	  date      = c(as.numeric(date_ete))
	  abundance = c(as.numeric(abondance_ete))
	  species   = c(as.character(species_ete))
	  dataF_ete = data.frame(date,abundance,species)
	  colnames(dataF_ete)=c("date","abundance","species") 
	  vec_synchrony_Loreau_ete[l]= synchrony(dataF_ete, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Loreau")
	  vec_synchrony_Gross_ete[l]=synchrony(dataF_ete, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Gross")
	  #### winter
	  abondance_hiver=c()
	  date_hiver=c()
	  species_hiver=c()
	  for(y in min(matrice$Annee):max(matrice$Annee)-1){
		 subdata_hiver = subset(matrice,(as.numeric(format(matrice$Date, format = "%m")) %in% hivern[[l]] & as.numeric(format(matrice$Date, format = "%Y"))==y) | (as.numeric(format(matrice$Date, format = "%m")) %in% hivern1[[l]] & as.numeric(format(matrice$Date, format = "%Y"))==y+1))
		 if(dim(subdata_hiver)[1]>0){
		   species_en_cours = c(unique(as.character(subdata_hiver$Nom_latin)))
		   for (e in 1:length(species_en_cours)){
		     if(length(subdata_hiver$Nombre[subdata_hiver$Nom_latin==species_en_cours[e]])>0){
		       abondance_hiver=c(abondance_hiver,c(mean(subdata_hiver$Nombre[subdata_hiver$Nom_latin==species_en_cours[e]])))
		       date_hiver=c(date_hiver,c(as.Date(as.character(subdata_hiver$Date[1]),format="%Y-%m-%d"))) 
		       species_hiver =c(species_hiver,c(species_en_cours[e]))
		     }
		   }
		 }
	  }
	 
	  # winter
	  date      = c(as.numeric(date_hiver))
	  abundance = c(as.numeric(abondance_hiver))
	  species   = c(as.character(species_hiver))
	  dataF_hiver = data.frame(date,abundance,species)
	  colnames(dataF_hiver)=c("date","abundance","species") 
	  vec_synchrony_Loreau_hiver[l]= synchrony(dataF_hiver, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Loreau")
	  vec_synchrony_Gross_hiver[l]=synchrony(dataF_hiver, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Gross")
	}
	print (vec_synchrony_Loreau_ete)
	print (vec_synchrony_Gross_ete)
	print (vec_synchrony_Loreau_hiver)
	print (vec_synchrony_Gross_hiver)

	# Data compilation
	data_synchro =c(c(vec_synchrony_Loreau_ete),c(vec_synchrony_Gross_ete),c(vec_synchrony_Loreau_hiver),c(vec_synchrony_Gross_hiver))

	#The names of the columns! probably we can do this automatically, but I have not found the trick yet.
	type =list("E:4,5,6\nH:10,11,12","E:5,6,7\nH:11,12,1","E:6,7,8\nH:12,1,2","E:7,8,9\nH:1,2,3","E:4,5,6,7\nH:10,11,12,1",
		        "E:5,6,7,8\nH:11,12,1,2","E:6,7,8,9\nH:12,1,2,3","E:4,5,6,7,8\nH:10,11,12,1,2","E:5,6,7,8,9\nH:11,12,1,2,3",
		        "E:4,5,6,7,8,9\nH:10,11,12,1,2,3")

	# type =list("E:4,5,6","H:10,11,12", "E:5,6,7","H:11,12,1","E:6,7,8","H:12,1,2","E:7,8,9","H:1,2,3","E:4,5,6,7","H:10,11,12,1",
	#            "E:5,6,7,8","H:11,12,1,2","E:6,7,8,9","H:12,1,2,3","E:4,5,6,7,8","H:10,11,12,1,2","E:5,6,7,8,9","H:11,12,1,2,3",
	#            "E:4,5,6,7,8,9","H:10,11,12,1,2,3")
	# type =list("4,5,6","10,11,12", "5,6,7","11,12,1", "6,7,8","12,1,2","7,8,9","1,2,3","4,5,6,7","10,11,12,1",
	#            "5,6,7,8","11,12,1,2","6,7,8,9","12,1,2,3","4,5,6,7,8","10,11,12,1,2","5,6,7,8,9","11,12,1,2,3",
	#            "4,5,6,7,8,9","10,11,12,1,2,3")

	max_value = max(data_synchro,na.rm=TRUE)+0.002
	#data_synchro = matrix(data_synchro,nc=length(ete)*2, byrow=T)
	data_synchro = matrix(data_synchro,nc=length(ete), byrow=T)
	colnames(data_synchro) = type
	plot.new()
	pdf(file_out,width=25,height=10)
	par(mar=c(4,4.5,3,2.5))
	par(oma = c(4.2, 0.5, 0.5, 0.5))
	#bp = barplot(data_synchro, beside=T, col=c("#FFFFAA","#AAFFAA","red","blue"),ylim=c(0,10))
	barplot(data_synchro, xlab="X",ylab="Y", beside=T, col=c("#FFFFAA","#AAFFAA","red","blue"), ylim=c(0,max_value), main=titre,cex.lab=1,cex.axis = 1,cex.main=2) # Draw the bars
	box() #Frame the diagram
	par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
	plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
	legend("bottom",c("Loreau-Summer","Gross-Summer","Loreau-Winter","Gross-Winter"), col=c("#FFFFAA","#AAFFAA","red","blue"), xpd = TRUE,  horiz=TRUE,
		    inset=.02,bty = "n",cex=1.5, pch=15)
	dev.off()
	return (data_synchro)
}
#-------------------------------------------------
## - call function number T-42
type =c("E:4,5,6\nH:10,11,12","E:5,6,7\nH:11,12,1","E:6,7,8\nH:12,1,2","E:7,8,9\nH:1,2,3","E:4,5,6,7\nH:10,11,12,1",
           "E:5,6,7,8\nH:11,12,1,2","E:6,7,8,9\nH:12,1,2,3","E:4,5,6,7,8\nH:10,11,12,1,2","E:5,6,7,8,9\nH:11,12,1,2,3",
           "E:4,5,6,7,8,9\nH:10,11,12,1,2,3")

Tab_synchronySaison=matrix(ncol = 10)
colnames(Tab_synchronySaison)=type
# All birds, all years
data_synchro = SynchronySeason_moving_month(DBt,"T42-Axe1-Teich-Synchrony_Season_Comparison_All_Species_moving_month.pdf","Synchrony by season (month moving) all species considered (1973-2016)")
rownames(data_synchro)= c('SummerLoreau_AllSpecieAllDate','SummerGross_AllSpecieAllDate','WinterLoreau_AllSpecieAllDate','WinterGross_AllSpecieAllDate')
Tab_synchronySaison = rbind(Tab_synchronySaison,data_synchro)

# All birds from 2006
tempDBt = subset(DBt,as.numeric(format(DBt$Date, format = "%Y"))>=2006)
data_synchro = SynchronySeason_moving_month(tempDBt,"T42-Axe1-Teich-Synchrony_Season_Comparison_All_Species_moving_month_From2006.pdf","Synchrony by season (month moving) all species considered (2006-2016)")
rownames(data_synchro)= c('SummerLoreau_AllSpecieFrom2006','SummerGross_AllSpecieFrom2006','WinterLoreau_AllSpecieFrom2006','WinterGross_AllSpecieFrom2006')
Tab_synchronySaison = rbind(Tab_synchronySaison,data_synchro)

# All Birds before 2006
tempDBt = subset(DBt,as.numeric(format(DBt$Date, format = "%Y"))<2006)
data_synchro = SynchronySeason_moving_month(tempDBt,"T42-Axe1-Teich-Synchrony_Season_Comparison_All_Species_moving_month_Before2006.pdf","Synchrony by season (month moving) all species considered (1973-2005)")
rownames(data_synchro)= c('SummerLoreau_AllSpecieBefore2006','SummerGross_AllSpecieBefore2006','WinterLoreau_AllSpecieBefore2006','WinterGross_AllSpecieBefore2006')
Tab_synchronySaison = rbind(Tab_synchronySaison,data_synchro)

# All Frequent Birds
tempDBt = subset(DBt,as.character(DBt$Nom_latin)%in% oiseaux_Frequents_t )
data_synchro = SynchronySeason_moving_month(tempDBt,"T42-Axe1-Teich-Synchrony_Season_Comparison_frequentSpecies_moving_month.pdf","Synchrony by season (month moving), frequent species(1973-2016)")
rownames(data_synchro)= c('SummerLoreau_FreqBirdsAllDate','SummerGross_FreqBirdsAllDate','WinterLoreau_FreqBirdsAllDate','WinterGross_FreqBirdsAllDate')
Tab_synchronySaison = rbind(Tab_synchronySaison,data_synchro)

# Frequent Birds from 2006
tempDBt = subset(DBt,as.character(DBt$Nom_latin)%in% oiseaux_Frequents_t & as.numeric(format(DBt$Date, format = "%Y"))>=2006)
data_synchro = SynchronySeason_moving_month(tempDBt,"T42-Axe1-Teich-Synchrony_Season_Comparison_frequentSpecies_moving_month_From2006.pdf","Synchrony by season (month moving), frequent species(2006-2016)")
rownames(data_synchro)= c('SummerLoreau_FreqBirdsFrom2006','SummerGross_FreqBirdsFrom2006','WinterLoreau_FreqBirdsFrom2006','WinterGross_FreqBirdsFrom2006')
Tab_synchronySaison = rbind(Tab_synchronySaison,data_synchro)

# Frequent Birds Before 2006
tempDBt = subset(DBt,as.character(DBt$Nom_latin)%in% oiseaux_Frequents_t & as.numeric(format(DBt$Date, format = "%Y"))<2006)
data_synchro = SynchronySeason_moving_month(tempDBt,"T42-Axe1-Teich-Synchrony_Season_Comparison_frequentSpecies_moving_month_Before2006.pdf","Synchrony by season (month moving), frequent species(1973-2006)")
rownames(data_synchro)= c('SummerLoreau_FreqBirdsBefore2006','SummerGross_FreqBirdsBefore2006','WinterLoreau_FreqBirdsBefore2006','WinterGross_FreqBirdsBefore2006')
Tab_synchronySaison = rbind(Tab_synchronySaison,data_synchro)

# All the waders
tempDBt = subset(DBt,as.character(DBt$Nom_latin)%in% limicoles )
data_synchro = SynchronySeason_moving_month(tempDBt,"T42-Axe1-Teich-Synchrony_Season_Comparison_wadingBirds_moving_month.pdf","Synchrony by season (month moving), wading Birds(1973-2016)")
rownames(data_synchro)= c('SummerLoreauLimicole','SummerGrossLimicole','WinterLoreauLimicole','WinterGrossLimicole')
Tab_synchronySaison = rbind(Tab_synchronySaison,data_synchro)

# Waders from 2006
tempDBt = subset(DBt,as.character(DBt$Nom_latin)%in% limicoles & as.numeric(format(DBt$Date, format = "%Y"))>=2006)
data_synchro = SynchronySeason_moving_month(tempDBt,"T42-Axe1-Teich-Synchrony_Season_Comparison_wadingBirds_moving_month_From2006.pdf","Synchrony by season (month moving), wading Birds(2006-2016)")
rownames(data_synchro)= c('SummerLoreauLimicole_From2006','SummerGrossLimicole_From2006','WinterLoreauLimicole_From2006','WinterGrossLimicole_From2006')
Tab_synchronySaison = rbind(Tab_synchronySaison,data_synchro)

# Waders before 2006
tempDBt = subset(DBt,as.character(DBt$Nom_latin)%in% limicoles & as.numeric(format(DBt$Date, format = "%Y"))<2006 )
data_synchro = SynchronySeason_moving_month(tempDBt,"T42-Axe1-Teich-Synchrony_Season_Comparison_wadingBirds_moving_month_Before2006.pdf","Synchrony by season (month moving), wading Birds(1973-2006)")
rownames(data_synchro)= c('SummerLoreauLimicole_Before2006','SummerGrossLimicole_Before2006','WinterLoreauLimicole_Before2006','WinterGrossLimicole_Before2006')
Tab_synchronySaison = rbind(Tab_synchronySaison,data_synchro)

##write all the values in a table
write.csv(Tab_synchronySaison, file = "T42-Axe1-Teich-Synchrony_Season_Comparison_Tab_synchronySaison.csv")

# ---------------------------------------------------------------------
#
#                              Values in a table
#
# ---------------------------------------------------------------------
Tab_synchronySaisonHiver = subset(Tab_synchronySaison,grepl('Winter',rownames(Tab_synchronySaison))    )
type_winter =c("H:10,11,12","H:11,12,1","H:12,1,2","H:1,2,3","H:10,11,12,1",
        "H:11,12,1,2","H:12,1,2,3","H:10,11,12,1,2","H:11,12,1,2,3",
        "H:10,11,12,1,2,3")
colnames(Tab_synchronySaisonHiver)=type_winter
barplot(Tab_synchronySaisonHiver, xlab="X",ylab="Y", beside=T, main='titre',cex.lab=1,cex.axis = 1,cex.main=2) 
Tab_synchronySaisonEte = subset(Tab_synchronySaison,grepl('Summer',rownames(Tab_synchronySaison))    )
type_summer =c("E:4,5,6","E:5,6,7","E:6,7,8","E:7,8,9","E:4,5,6,7",
        "E:5,6,7,8","E:6,7,8,9","E:4,5,6,7,8","E:5,6,7,8,9",
        "E:4,5,6,7,8,9")
colnames(Tab_synchronySaisonEte)=type_summer
barplot(Tab_synchronySaisonEte, xlab="X",ylab="Y", beside=T, main='titre',cex.lab=1,cex.axis = 1,cex.main=2) 
couleur = c('green','red','black','blue','#F0C300','pink','lightgreen','violet','lightblue','darkred','#6600FF','orange','#BABABA','#40826D','#3A9D23','#C71585','#F88E55','#EDFF0C')



# ---------------------------------------------------------------------
#
#                              SCRIPT T-43
#
# ---------------------------------------------------------------------
Tab_synchronySaisonLoreau = subset(Tab_synchronySaison,grepl('Loreau',rownames(Tab_synchronySaison))    )
Tab_synchronySaisonGross = subset(Tab_synchronySaison,grepl('Gross',rownames(Tab_synchronySaison))    )
# Only Frequent Birds
Tab_synchronySaisonLoreau = subset(Tab_synchronySaison,grepl('Loreau',rownames(Tab_synchronySaison)) & grepl('FreqBirds',rownames(Tab_synchronySaison)) )
Tab_synchronySaisonGross = subset(Tab_synchronySaison,grepl('Gross',rownames(Tab_synchronySaison)) & grepl('FreqBirds',rownames(Tab_synchronySaison))    )
pdf("T43-Axe1-Teich-Synchrony_Season_Comparison_Tab_synchronySaison_frequent_birds.pdf",width=30,height=15)
par(mar=c(4,4.5,3,2.5))
par(oma = c(4.2, 0.5, 0.5, 0.5))
barplot(Tab_synchronySaisonLoreau, xlab="X",ylab="Value of Synchrony",col=(couleur[1:length(rownames(Tab_synchronySaisonLoreau))]),beside=T, main='Evolution of synchrony values for Loreau, for different months, and different datasets, only for frequent birds',cex.lab=1,cex.axis = 1,cex.main=2,space = c(0,3)) 
box() # Frame the diagram
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom",rownames(Tab_synchronySaisonLoreau), col=couleur[1:length(rownames(Tab_synchronySaisonLoreau))], xpd = TRUE, ncol=6,inset=.02,bty = "n",cex=1.2, pch=15)

# GROSS
par(mar=c(4,4.5,3,2.5))
par(oma = c(4.2, 0.5, 0.5, 0.5))
barplot(Tab_synchronySaisonGross, xlab="X",ylab="Value of Synchrony",col=(couleur[1:length(rownames(Tab_synchronySaisonGross))]), main='Evolution of synchrony values for Gross, for different months, and different datasets, only for frequent birds',cex.lab=1,cex.axis = 1,cex.main=2,space = c(0,3),beside=T) 
box() # Frame the diagram
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom",rownames(Tab_synchronySaisonGross), col=(couleur[1:length(rownames(Tab_synchronySaisonGross))]), xpd = TRUE, ncol=6,inset=.02,bty = "n",cex=1.2, pch=15)
dev.off()

# ---------------------------------------------------------------------
#
#                              SCRIPT T-44
#
# ---------------------------------------------------------------------
library(pracma)
# abondance, calidris, anas, waders.
abondance_anas = subset(DBt,grepl("^Anas",DBt$Nom_latin,ignore.case = TRUE))
abondance_calidris = subset(DBt,grepl("^Calidris",DBt$Nom_latin,ignore.case = TRUE))
abondance_waders = subset(DBt, DBt$Nom_latin %in% limicoles)
pdf("T44-Axe1-Teich-Abondance_of_waders_anas_and_calidris.pdf",width=15,height=10)
# rather make three graphs!!!
par(mfrow=c(3,1)) 
minD = min(abondance_waders$Date,abondance_anas$Date,abondance_calidris$Date)
maxD = max(abondance_waders$Date,abondance_anas$Date,abondance_calidris$Date)
ticks = date(paste(unique(as.numeric(format(unique(DBt$Date), format = "%Y"))),"-01-01",sep=""))
ticks = ticks[seq(1,length(ticks),by = 2)]
nameTicks=unique(as.numeric(format(unique(DBt$Date), format = "%Y")))
nameTicks = nameTicks[seq(1,length(nameTicks),by = 2)]
plot(abondance_waders$Date,abondance_waders$Nombre,type="o",col="black",xlab="Date",ylab="Abondance of all Waders",xlim=c(minD,maxD),main='Abundance of waders in the dataset',cex.lab=1.2,cex.main=1.5,xaxt="n")
axis(1, at=ticks, labels=nameTicks)
plot(abondance_calidris$Date,abondance_calidris$Nombre,type="o",col="blue",xlab="Date",ylab="Abondance of all Calidris",xlim=c(minD,maxD),main='Abundance of Calidris in the dataset',cex.lab=1.2,cex.main=1.5,xaxt="n")
axis(1, at=ticks, labels=nameTicks)
plot(abondance_anas$Date,abondance_anas$Nombre,type="o",col="green",xlab="Date",ylab="Abondance of all Anas",xlim=c(minD,maxD),main='Abundance of Anas in the dataset',cex.lab=1.2,cex.main=1.5,xaxt="n")
axis(1, at=ticks, labels=nameTicks)

## with log
# rather make three graphs!!!
par(mfrow=c(3,1)) 

plot(abondance_waders$Date,log10(abondance_waders$Nombre),type="o",col="black",xlab="Date",ylab="Abondance of all Waders",xlim=c(minD,maxD),main='Log of abundance of waders in the dataset',cex.lab=1.2,cex.main=1.5,xaxt="n")
axis(1, at=ticks, labels=nameTicks)
plot(abondance_calidris$Date,log10(abondance_calidris$Nombre),type="o",col="blue",xlab="Date",ylab="Abondance of all Calidris",xlim=c(minD,maxD),main='Log of abundance of Calidris in the dataset',cex.lab=1.2,cex.main=1.5,xaxt="n")
axis(1, at=ticks, labels=nameTicks)
plot(abondance_anas$Date,log10(abondance_anas$Nombre),type="o",col="green",xlab="Date",ylab="Abondance of all Anas",xlim=c(minD,maxD),main='Log of abundance of Anas in the dataset',cex.lab=1.2,cex.main=1.5,xaxt="n")
axis(1, at=ticks, labels=nameTicks)


## smoothies average
# rather make three graphs!!!
par(mfrow=c(3,1)) 

plot(abondance_waders$Date,abondance_waders$Nombre,col="black",xlab="Date",ylab="Abondance of all Waders",xlim=c(minD,maxD),main='Abundance of waders in the dataset',cex.lab=1.2,cex.main=1.5,xaxt="n")
lines(abondance_waders$Date,movavg(abondance_waders$Nombre,3,type='t'),col="red")
axis(1, at=ticks, labels=nameTicks)

plot(abondance_calidris$Date,movavg(abondance_calidris$Nombre,3,type='s'),col="blue",xlab="Date",ylab="Abondance of all Calidris",xlim=c(minD,maxD),main='Abundance of Calidris in the dataset',cex.lab=1.2,cex.main=1.5,xaxt="n")
lines(abondance_calidris$Date,movavg(abondance_calidris$Nombre,3,type='t'),col="red")
axis(1, at=ticks, labels=nameTicks)


plot(abondance_anas$Date,movavg(abondance_anas$Nombre,3,type='s'),col="green",xlab="Date",ylab="Abondance of all Anas",xlim=c(minD,maxD),main='Abundance of Anas in the dataset',cex.lab=1.2,cex.main=1.5,xaxt="n")
axis(1, at=ticks, labels=nameTicks)


## smoothies average + log
# rather make three graphs!!!
par(mfrow=c(3,1)) 
minD=date("1981-01-01")
plot(abondance_waders$Date,log10(abondance_waders$Nombre),col="black",xlab="Date",ylab="Abondance of all Waders",xlim=c(minD,maxD),main='Log of abundance of waders in the dataset',cex.lab=1.2,cex.main=1.5,xaxt="n",pch=19)
lines(abondance_waders$Date,movavg(log10(abondance_waders$Nombre),3,type='t'),col="orange",lwd=1.5)
axis(1, at=ticks, labels=nameTicks)
#
plot(abondance_calidris$Date,log10(abondance_calidris$Nombre),col="blue",xlab="Date",ylab="Abondance of all Calidris",xlim=c(minD,maxD),main='Log of abundance of Calidris in the dataset',cex.lab=1.2,cex.main=1.5,xaxt="n",pch=19)
lines(abondance_calidris$Date,movavg(log10(abondance_calidris$Nombre),3,type='t'),col="orange",lwd=1.5)
axis(1, at=ticks, labels=nameTicks)
#
plot(abondance_anas$Date,log10(abondance_anas$Nombre),col="green",xlab="Date",ylab="Abondance of all Anas",xlim=c(minD,maxD),main='Log of abundance of Anas in the dataset',cex.lab=1.2,cex.main=1.5,xaxt="n",pch=19)
lines(abondance_anas$Date,movavg(log10(abondance_anas$Nombre),3,type='t'),col="orange",lwd=1.5)
axis(1, at=ticks, labels=nameTicks)


dev.off()

# ---------------------------------------------------------------------
#
#                              SCRIPT T-45
#
# ---------------------------------------------------------------------




abondance_anas = subset(DBt,grepl("^Anas",DBt$Nom_latin,ignore.case = TRUE))
abondance_calidris = subset(DBt,grepl("^Calidris",DBt$Nom_latin,ignore.case = TRUE))
abondance_waders = subset(DBt, DBt$Nom_latin %in% limicoles)
ticks = date(paste(unique(as.numeric(format(unique(DBt$Date), format = "%Y"))),"-01-01",sep=""))
ticks = ticks[seq(1,length(ticks),by = 2)]
nameTicks=unique(as.numeric(format(unique(DBt$Date), format = "%Y")))
nameTicks = nameTicks[seq(1,length(nameTicks),by = 2)]
listColor = c('green','red','blue','#F0C300','pink','lightgreen','violet','lightblue','darkred','#6600FF','orange','#BABABA','#40826D','#3A9D23','#C71585','#F88E55','#EDFF0C')
# start anas
temp_species = unique(as.character(abondance_anas$Nom_latin))
print (paste("Nombre d'espèces",length(temp_species)))

### avec le log
ymin = 0
ymax = log(max(abondance_anas$Nombre))
xmin = min(abondance_anas$Date)
xmax = max(abondance_anas$Date)
plot(abondance_anas$Date,log(abondance_anas$Nombre),col="black",xlab="Date",ylab="Log abundance",xlim=c(xmin,xmax),main='a)',cex.lab=1.2,cex.main=1.5,xaxt="n",pch=19,ylim=c(ymin,ymax),type="l")
for (s in 1:length(temp_species)){
  points(abondance_anas$Date[as.character(abondance_anas$Nom_latin) == temp_species[s]],log(abondance_anas$Nombre[abondance_anas$Nom_latin==temp_species[s]]),col=listColor[s],pch=19)
}
axis(1, at=ticks, labels=nameTicks)



# ---------------------------------------------------------------------
#
#                              SCRIPT T-45
#
# ---------------------------------------------------------------------
library(plotrix)
create_somme = function(colonneDate,colonneResultats){
  minAnnee = as.numeric(format(min(colonneDate), format = "%Y"))
  maxAnnee = as.numeric(format(max(colonneDate), format = "%Y"))
  nb = (maxAnnee-minAnnee)+1
  mat <- matrix(data=0:0,nrow=nb, ncol=12) #12 : 1col by month
  rownames(mat)=minAnnee:maxAnnee
  colnames(mat)=1:12
  for (i in 1:length(colonneDate)){
    annee = format(colonneDate[i], format = "%Y") #extraire un mois d'une date
    annee = as.numeric(annee)
    mois  = format(colonneDate[i], format = "%m") #extraire un mois d'une date
    mois = as.numeric(mois)
    #print (annee,mois)
    mat[as.character(annee),as.character(mois)]=mat[as.character(annee),as.character(mois)]+colonneResultats[i]
  }
  return (mat)
}



create_mat_percent=function(matrice,i_abondance_tot){
  species = unique(as.character(matrice$Nom_latin))
  mat <- matrix(data=0:0,nrow=length(species), ncol=3) 
  rownames(mat)=species
  colnames(mat)=c('somme','indice','percent')
  for (i in 1:length(species)){
    mat[species[i],'somme']=sum(matrice$Nombre[as.character(matrice$Nom_latin)==species[i]],na.rm=TRUE) # I crush seasonality
    mat[species[i],'indice']=sum(matrice$Nombre[as.character(matrice$Nom_latin)==species[i]],na.rm=TRUE)/i_abondance_tot # I crush seasonality
    mat[species[i],'percent']=(sum(matrice$Nombre[as.character(matrice$Nom_latin)==species[i]],na.rm=TRUE)/i_abondance_tot)*100 # I crush seasonality
  }
  return (mat)
}

d_ticks = date(paste(unique(as.numeric(format(unique(DBt$Date), format = "%Y"))),"-01-01",sep=""))
d_ticks = d_ticks[seq(1,length(d_ticks),by = 2)]
nameTicks=unique(as.numeric(format(unique(DBt$Date), format = "%Y")))
nameTicks = nameTicks[seq(1,length(nameTicks),by = 2)]
listColor = c('darkgreen','red','blue','#F0C300','pink','lightgreen','violet','lightblue','darkred','#6600FF','orange','#BABABA','#40826D','#3A9D23','#C71585','#F88E55','#00ffff','yellow','darkgrey','royalblue1','chocolate','antiquewhite1','black','deeppink','lightcoral','mediumorchid1','mediumaquamarine','olivedrab1','orangered','springgreen1','thistle1')


# The doughnut function permits to draw a donut plot (source : https://www.r-graph-gallery.com/130-ring-or-donut-chart/)
# I added option label.cex for change size of labels
doughnut <-
  function (x, labels = names(x), edges = 200, outer.radius = 0.8, 
            inner.radius=0.6, clockwise = FALSE,
            init.angle = if (clockwise) 90 else 0, density = NULL, 
            angle = 45, col = NULL, border = FALSE, lty = NULL, 
            main = NULL, label.cex=1,coord.x=c(1,1),...)
  {
    if (!is.numeric(x) || any(is.na(x) | x < 0))
      stop("'x' values must be positive.")
    if (is.null(labels))
      labels <- as.character(seq_along(x))
    else labels <- as.graphicsAnnot(labels)
    x <- c(0, cumsum(x)/sum(x))
    dx <- diff(x)
    nx <- length(dx)
    plot.new()
    pin <- par("pin")
    xlim <-coord.x
    ylim <- c(-1, 1)
    if (pin[1L] > pin[2L])
      xlim <- (pin[1L]/pin[2L]) * xlim
    else ylim <- (pin[2L]/pin[1L]) * ylim
    plot.window(xlim, ylim, "", asp = 1)
    if (is.null(col))
      col <- if (is.null(density))
        palette()
    else par("fg")
    col <- rep(col, length.out = nx)
    border <- rep(border, length.out = nx)
    lty <- rep(lty, length.out = nx)
    angle <- rep(angle, length.out = nx)
    density <- rep(density, length.out = nx)
    twopi <- if (clockwise)
      -2 * pi
    else 2 * pi
    t2xy <- function(t, radius) {
      t2p <- twopi * t + init.angle * pi/180
      list(x = radius * cos(t2p), 
           y = radius * sin(t2p))
    }
    for (i in 1L:nx) {
      n <- max(2, floor(edges * dx[i]))
      P <- t2xy(seq.int(x[i], x[i + 1], length.out = n),
                outer.radius)
      polygon(c(P$x, 0), c(P$y, 0), density = density[i], 
              angle = angle[i], border = border[i], 
              col = col[i], lty = lty[i])
      Pout <- t2xy(mean(x[i + 0:1]), outer.radius)
      lab <- as.character(labels[i])
      if (!is.na(lab) && nzchar(lab)) {
        lines(c(1, 1.05) * Pout$x, c(1, 1.05) * Pout$y)
        text(1.1 * Pout$x, 1.1 * Pout$y, labels[i], 
             xpd = TRUE, adj = ifelse(Pout$x < 0, 1, 0), cex=label.cex,
             #font=2,
             ...)
      }
      ## Add white disc          
      Pin <- t2xy(seq.int(0, 1, length.out = n*nx),
                  inner.radius)
      polygon(Pin$x, Pin$y, density = density[i], 
              angle = angle[i], border = border[i], 
              col = "white", lty = lty[i])
    }
    
    title(main = main, ...)
    invisible(NULL)
  }

# Let's use the function, it works like PiePlot !
# inner.radius controls the width of the ring!
#doughnut( c(3,5,9,12) , inner.radius=0.5, col=c(rgb(0.2,0.2,0.4,0.5), rgb(0.8,0.2,0.4,0.5), rgb(0.2,0.9,0.4,0.4) , rgb(0.0,0.9,0.8,0.4)) )

pdf("T45-Axe1-Teich-Abondance_of_waders_anas_and_calidris_vers2.pdf",width=15,height=10,family = 'Times')
# rather make three graphs!!!
par(mfrow=c(3,2))
layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = TRUE),heights = c(1,1,1),widths = c(3,1))
abondance_anas = subset(DBt,grepl("^Anas",DBt$Nom_latin,ignore.case = TRUE) & DBt$Annee>=1981)
abondance_calidris = subset(DBt,grepl("^Calidris",DBt$Nom_latin,ignore.case = TRUE)& DBt$Annee>=1981)
abondance_waders = subset(DBt, DBt$Nom_latin %in% limicoles & DBt$Annee>=1981)


## abondance totale
abondance_tot_anas = sum(abondance_anas$Nombre,na.rm=TRUE) # I crush seasonality
abondance_tot_calidris = sum(abondance_calidris$Nombre,na.rm=TRUE) 
abondance_tot_waders = sum(abondance_waders$Nombre,na.rm=TRUE) 

# matrice percent abondance without seasonality
mat_percent_abondance_tot_anas     = create_mat_percent(abondance_anas,abondance_tot_anas)
mat_percent_abondance_tot_calidris = create_mat_percent(abondance_calidris,abondance_tot_calidris)
mat_percent_abondance_tot_waders   = create_mat_percent(abondance_waders,abondance_tot_waders)

minimal_percent = 1.00001 # percent minimal for considere the species
# Anas First
mat=create_somme(abondance_anas$Date,abondance_anas$Nombre)
temp_date = expand.grid(15,colnames(mat),rownames(mat)) #15 because i need a day for build complete date.
temp_date = c(paste(temp_date[,1],temp_date[,2],temp_date[,3],sep='-'))
temp_date = as.Date(temp_date,format="%d-%m-%Y")
temp_value = c(t(mat))
temp_species = unique(as.character(abondance_anas$Nom_latin))
temp_species = temp_species[order(temp_species)]
ymin = 0
ymax = log(max(abondance_anas$Nombre))+0.5 #tocorrect
xmin = date("1981-01-01")
xmax = date("2016-04-01")
par(mar=c(1.5,4.3,0.1,1))
par(oma=c(0,0,0,0))

plot(temp_date,log(temp_value),col="black",xlab="",ylab="Log abundance",xlim=c(xmin,xmax),cex.lab=2,cex.axis=1.8,xaxt="n",pch=19,ylim=c(ymin,ymax),type="l",lwd =2)
lines(temp_date,log(temp_value+1),col="black",lwd =1,lty=3)
for (s in 1:length(temp_species)){
  points(abondance_anas$Date[as.character(abondance_anas$Nom_latin) == temp_species[s]],log(abondance_anas$Nombre[abondance_anas$Nom_latin==temp_species[s]]),col=listColor[s],pch=19)
}
#text(c(0,1), "a)",cex = 2)
mtext("\n (a)",side=3,cex = 1.5, adj = 0,padj=1)
axis(1, at=d_ticks, cex=1.8, ticks=TRUE, labels=FALSE) # pas de label

par(mar=c(0.1,0.1,0.1,0.1))
mat_percent_abondance_tot_anas = mat_percent_abondance_tot_anas[ order(row.names(mat_percent_abondance_tot_anas)), ]

doughnut( c(mat_percent_abondance_tot_anas[,'percent'][round(mat_percent_abondance_tot_anas[,'percent'])>=minimal_percent]) , inner.radius=0.15, col=listColor[which(round(mat_percent_abondance_tot_anas[,'percent'])>=minimal_percent)],outer.radius =0.45,label.cex=1.3,coord.x=c(-0.5,1) )

# next Calidris
par(mar=c(1.5,4.3,0.1,1))
mat=create_somme(abondance_calidris$Date,abondance_calidris$Nombre)
temp_date = expand.grid(15,colnames(mat),rownames(mat)) #15 because i need a day for build complete date.
temp_date = c(paste(temp_date[,1],temp_date[,2],temp_date[,3],sep='-'))
temp_date = as.Date(temp_date,format="%d-%m-%Y")
temp_value = c(t(mat))
temp_species = unique(as.character(abondance_calidris$Nom_latin))
temp_species = temp_species[order(temp_species)]
ymin = 0
ymax = log(max(abondance_calidris$Nombre))

plot(temp_date,log(temp_value),col="black",xlab="",ylab="Log abundance",xlim=c(xmin,xmax),cex.axis=1.8,cex.lab=2,xaxt="n",pch=19,ylim=c(ymin,ymax),type="l",lwd =2)
lines(temp_date,log(temp_value+1),col="black",lwd =1,lty=3)
for (s in 1:length(temp_species)){
  points(abondance_calidris$Date[as.character(abondance_calidris$Nom_latin) == temp_species[s]],log(abondance_calidris$Nombre[abondance_calidris$Nom_latin==temp_species[s]]),col=listColor[s],pch=19)
}
axis(1, at=d_ticks, cex=1.8, ticks=TRUE, labels=FALSE) # pas de label
mtext("\n (b)",side=3,cex = 1.5, adj = 0,padj = 1)
# plot.new()
par(mar=c(0.1,0.1,0.1,0.1))
mat_percent_abondance_tot_calidris = mat_percent_abondance_tot_calidris[ order(row.names(mat_percent_abondance_tot_calidris)), ]


doughnut( c(mat_percent_abondance_tot_calidris[,'percent'][round(mat_percent_abondance_tot_calidris[,'percent'])>=minimal_percent]) , inner.radius=0.145, col=listColor[which(round(mat_percent_abondance_tot_calidris[,'percent'])>=minimal_percent)],outer.radius =0.45,label.cex=1.28,coord.x=c(-0.7,1) )

# next waders
par(mar=c(4,4.3,0,1))
mat=create_somme(abondance_waders$Date,abondance_waders$Nombre)
temp_date = expand.grid(15,colnames(mat),rownames(mat)) #15 because i need a day for build complete date.
temp_date = c(paste(temp_date[,1],temp_date[,2],temp_date[,3],sep='-'))
temp_date = as.Date(temp_date,format="%d-%m-%Y")
temp_value = c(t(mat))
temp_species = unique(as.character(abondance_waders$Nom_latin))
temp_species = temp_species[order(temp_species)]
ymin = 0
ymax = log(max(abondance_waders$Nombre))
plot(temp_date,log(temp_value),col="black",xlab="Year",ylab="Log abundance",xlim=c(xmin,xmax),cex.axis=1.8,cex.lab=2,xaxt="n",pch=19,ylim=c(ymin,ymax),type="l",lwd =2)
lines(temp_date,log(temp_value+1),col="black",lwd =1,lty=3)
for (s in 1:length(temp_species)){
  points(abondance_waders$Date[as.character(abondance_waders$Nom_latin) == temp_species[s]],log(abondance_waders$Nombre[abondance_waders$Nom_latin==temp_species[s]]),col=listColor[s],pch=19)
}
mtext("\n (c)",side=3,cex = 1.5, adj = 0,padj = 1)
axis(1, at=d_ticks, labels=nameTicks,cex=1.8,cex.axis=1.8)
mat_percent_abondance_tot_waders = mat_percent_abondance_tot_waders[ order(row.names(mat_percent_abondance_tot_waders)), ]
par(mar=c(0.1,0.1,0.1,0.1))

doughnut( c(mat_percent_abondance_tot_waders[,'percent'][round(mat_percent_abondance_tot_waders[,'percent'])>=minimal_percent]) , inner.radius=0.15, col=listColor[which(round(mat_percent_abondance_tot_waders[,'percent'])>=minimal_percent)],outer.radius =0.45,label.cex=1.28,coord.x=c(-0.5,1),lty=1 )
dev.off()


pdf("T47-Axe1-Teich-les_3_especes.pdf",width=20,height=6,family = 'Times')
par(mar=c(4,5,3,1))
d_ticks = date(paste(unique(as.numeric(format(unique(DBt$Date), format = "%Y"))),"-01-01",sep=""))
d_ticks = d_ticks[seq(1,length(d_ticks),by = 2)]
nameTicks=unique(as.numeric(format(unique(DBt$Date), format = "%Y")))
nameTicks = nameTicks[seq(1,length(nameTicks),by = 2)]
abondance_egretta = subset(DBt,DBt$Nom_latin=='Egretta garzetta' & DBt$Annee>=1981)
abondance_phalacrocorax = subset(DBt,DBt$Nom_latin=='Phalacrocorax carbo' & DBt$Annee>=1981)
abondance_ardea = subset(DBt,DBt$Nom_latin=='Ardea cinerea' & DBt$Annee>=1981)
ymin = 0
ymax = max(abondance_egretta$Nombre,abondance_phalacrocorax$Nombre,abondance_ardea$Nombre,na.rm=TRUE)
ymax = log(max(abondance_egretta$Nombre,abondance_phalacrocorax$Nombre,abondance_ardea$Nombre,na.rm=TRUE))
xmin = date("1981-01-01")
xmax = date("2016-04-01")
plot(abondance_egretta$Date,log(abondance_egretta$Nombre),col="darkgreen",xlab="Year",ylab="Log abundance",xlim=c(xmin,xmax),cex.lab=1.8,cex.axis=1.5,xaxt="n",pch=19,ylim=c(ymin,ymax),type="o",lwd =1.8,cex=0.9)
lines(abondance_egretta$Date,log(abondance_egretta$Nombre+1),col="darkgreen",lwd =1,lty=3)
lines(abondance_phalacrocorax$Date,log(abondance_phalacrocorax$Nombre),col="darkred",type='o',lwd=1.8,pch=19,cex=0.9)
lines(abondance_phalacrocorax$Date,log(abondance_phalacrocorax$Nombre+1),col="darkred",lwd =1,lty=3)
lines(abondance_ardea$Date,log(abondance_ardea$Nombre),col="darkblue",type='o',lwd=1.8,pch=19,cex=0.9)
lines(abondance_ardea$Date,log(abondance_ardea$Nombre+1),col="darkblue",lwd =1,lty=3)
axis(1, at=d_ticks, labels=nameTicks,cex.axis=1.5)
legend("topright",legend = c("Egretta garzetta",'Phalacrocorax carbo','Ardea cinerea'),col=c('darkgreen','darkred','darkblue'),pch=19,lwd=1)
#dev.off()

#pdf("T47B-Axe1-Teich-les_3_especes_Heron+Aigrette.pdf",width=20,height=6,family = 'Times')
HeronEgret =subset(DBt,((DBt$Nom_latin=='Egretta garzetta' | DBt$Nom_latin=='Ardea cinerea') &DBt$Annee>=1981))
matHeronEgret = create_somme(HeronEgret$Date,HeronEgret$Nombre)
temp_date = expand.grid(15,colnames(matHeronEgret),rownames(matHeronEgret)) #15 because i need a day for build complete date.
temp_date = c(paste(temp_date[,1],temp_date[,2],temp_date[,3],sep='-'))
temp_date = as.Date(temp_date,format="%d-%m-%Y")
temp_value = c(t(matHeronEgret))
temp_species = rep("Heron+Egret",length(temp_date))
df_HeronEgret = cbind(as.character(temp_date),temp_value,temp_species)
colnames(df_HeronEgret)=c("Date","Nombre","Nom_latin")
df_HeronEgret=as.data.frame(df_HeronEgret)
df_HeronEgret$Nombre = as.numeric(as.character(df_HeronEgret$Nombre))
df_HeronEgret$Date = as.Date(as.character(df_HeronEgret$Date))
df_HeronEgret$Nom_latin = as.character(df_HeronEgret$Nom_latin)
par(mar=c(4,5,3,1))
d_ticks = date(paste(unique(as.numeric(format(unique(DBt$Date), format = "%Y"))),"-01-01",sep=""))
d_ticks = d_ticks[seq(1,length(d_ticks),by = 2)]
nameTicks=unique(as.numeric(format(unique(DBt$Date), format = "%Y")))
nameTicks = nameTicks[seq(1,length(nameTicks),by = 2)]
abondance_phalacrocorax = subset(DBt,DBt$Nom_latin=='Phalacrocorax carbo' & DBt$Annee>=1981)
ymin = 0
ymax = max(df_HeronEgret$Nombre,abondance_phalacrocorax$Nombre,na.rm=TRUE)
ymax = log(max(df_HeronEgret$Nombre,abondance_phalacrocorax$Nombre,na.rm=TRUE))
xmin = date("1981-01-01")
xmax = date("2016-04-01")
plot(df_HeronEgret$Date,log(df_HeronEgret$Nombre),col="darkgreen",xlab="Year",ylab="Log abundance",xlim=c(xmin,xmax),cex.lab=1.8,cex.axis=1.5,xaxt="n",pch=19,ylim=c(ymin,ymax),type="o",lwd =1.8,cex=0.9)
lines(df_HeronEgret$Date,log(df_HeronEgret$Nombre+1),col="darkgreen",lwd =1,lty=3)
lines(abondance_phalacrocorax$Date,log(abondance_phalacrocorax$Nombre),col="darkred",type='o',lwd=1.8,pch=19,cex=0.9)
lines(abondance_phalacrocorax$Date,log(abondance_phalacrocorax$Nombre+1),col="darkred",lwd =1,lty=3)

axis(1, at=d_ticks, labels=nameTicks,cex.axis=1.5)
legend("topright",legend = c("Sum Egretta garzetta + Ardea cinerea",'Phalacrocorax carbo'),col=c('darkgreen','darkred'),pch=19,lwd=1)
dev.off()

#----- T48
# test rapide 
# Calcul of synchrony for the 3 species, without modification.
temp_data = subset(DBt,((DBt$Nom_latin=='Egretta garzetta' | DBt$Nom_latin=='Phalacrocorax carbo' | DBt$Nom_latin=='Ardea cinerea') &DBt$Annee>=1981))
temp_data$Nombre=as.numeric(as.character(temp_data$Nombre))
SynchronyMonth(matrice = temp_data,file_out = "T48-Axe1-Teich_Synchrony_for_the_3_competing_species.pdf",titre = "Gross synchrony calculation between \nEgretta garzetta, Phalacrocorax carbo and Ardea cinerea\n",loreau = FALSE,ymin = -0.5,ymax=1)

## I put Egretta and Ardea in a same subset. I used the data of Coralie(summed_abundances)
HeronEgret = subset(summed_abundances,summed_abundances$Nom_latin=='HeronEgret')
temp_data = subset(DBt,(DBt$Nom_latin=='Phalacrocorax carbo' &DBt$Annee>=1981))
temp_data=subset(temp_data, select=c("Nom_latin", "Date","Nombre"))
HeronEgret=subset(HeronEgret,select=c("Nom_latin", "Date","Nombre")) # just to be sure
temp_data=rbind(temp_data,HeronEgret) # for use the same function.
temp_data$Nombre=as.numeric(as.character(temp_data$Nombre))
SynchronyMonth(matrice = temp_data,file_out = "T48b-Axe1-Teich_Synchrony_for_the_3_competing_species_with_HeronEgret_summed.pdf",titre = "Gross synchrony calculation between \nPhalacrocorax carbo and sum of Ardea cinerea and Egretta garzetta\n",loreau = FALSE,ymin = -0.5,ymax=1)

# -----
# Sum of values of Heron and Aigrette, I create my own dataset for check
# Egretta and Ardea have the same dates.
# After 2007, Phalacrocorax has differents dates.
HeronEgret =subset(DBt,((DBt$Nom_latin=='Egretta garzetta' | DBt$Nom_latin=='Ardea cinerea') &DBt$Annee>=1981))
matHeronEgret = create_somme(HeronEgret$Date,HeronEgret$Nombre)
temp_date = expand.grid(15,colnames(matHeronEgret),rownames(matHeronEgret)) #15 because i need a day for build complete date.
temp_date = c(paste(temp_date[,1],temp_date[,2],temp_date[,3],sep='-'))
temp_date = as.Date(temp_date,format="%d-%m-%Y")
temp_value = c(t(matHeronEgret))
temp_species = rep("Heron+Egret",length(temp_date))
Cormoran = subset(DBt,(DBt$Nom_latin=='Phalacrocorax carbo' &DBt$Annee>=1981))
Cormoran=subset(Cormoran, select=c("Nom_latin", "Date","Nombre"))
df_HeronEgret = cbind(as.character(temp_date),temp_value,temp_species)
colnames(df_HeronEgret)=c("Date","Nombre","Nom_latin")
temp_matrice = rbind(Cormoran,df_HeronEgret)
Annee = format(temp_matrice$Date,format = "%Y")
temp_matrice =cbind(temp_matrice,Annee)
temp_matrice$Annee=as.numeric(as.character(temp_matrice$Annee))
temp_matrice$Nom_latin=as.character(temp_matrice$Nom_latin)
temp_matrice$Nombre=as.numeric(as.character(temp_matrice$Nombre))

#----- T49

#SynchronyMonth(matrice = temp_matrice,file_out = "T48c-Axe1-Teich_Synchrony_for_the_3_competing_species_with_HeronEgret_summed.pdf",titre = "Gross synchrony calculation between \nPhalacrocorax carbo and sum of Ardea cinerea and Egretta garzetta\n",loreau = FALSE,ymin = -0.5,ymax=1)
SynchronySeason2(matrice = temp_matrice,"T49-Axe1-Teich_Synchrony_by_season_for_the_3_competing_species.pdf",titre ="Gross synchrony calculation between \nPhalacrocorax carbo and sum of Ardea cinerea and Egretta garzetta\n",Loreau = FALSE,Gross=TRUE,max_value=0.8,min_value=-0.4)

#-------- 
# Sum of values of Heron and Aigrette, I create my own dataset for check
# Egretta and Ardea have the same dates.
# After 2007, Phalacrocorax has differents dates. In this test, i changed the days of this data, all data is artificially made on the 15th of the month. it is to see the impact of the day on the calculation of the index of Gross.
Cormoran = subset(DBt,(DBt$Nom_latin=='Phalacrocorax carbo' &DBt$Annee>=1981))
matCormoran = create_somme(Cormoran$Date,Cormoran$Nombre)
temp_date_C = expand.grid(15,colnames(matCormoran),rownames(matCormoran)) #15 because i need a day for build complete date.
temp_date_C = c(paste(temp_date_C[,1],temp_date_C[,2],temp_date_C[,3],sep='-'))
temp_date_C = as.Date(temp_date_C,format="%d-%m-%Y")
temp_value_C = c(t(matCormoran))
temp_species_C = rep("Phalacrocorax carbo",length(temp_date))
df_Cormoran = cbind(as.character(temp_date_C),temp_value_C,temp_species_C)
colnames(df_Cormoran)=c("Date","Nombre","Nom_latin")
temp_matrice = rbind(df_Cormoran,df_HeronEgret)
temp_matrice= as.data.frame(temp_matrice)
temp_matrice$Date=as.Date(as.character(temp_matrice$Date),format="%Y-%m-%d")
temp_matrice$Nombre=as.numeric(as.character(temp_matrice$Nombre))
#SynchronyMonth(matrice = temp_matrice,file_out = "T48d-Axe1-Teich_Synchrony_for_the_3_competing_species_with_HeronEgret_summed.pdf",titre = "Gross synchrony calculation between \nPhalacrocorax carbo and sum of Ardea cinerea and Egretta garzetta\n",loreau = FALSE,ymin = -0.5,ymax=1)


# TEST BY SEASON FOR THE SAME SPECIES.
# type =c("winter-all","summer-all","winter>=2006","summer>=2006","winter<2006","summer<2006")
# SynchronySeason(matrice = temp_matrice,"T49-Axe1-Teich_Synchrony_by_season_for_the_3_competing_species.pdf",titre ="Gross synchrony calculation between \nPhalacrocorax carbo and sum of Ardea cinerea and Egretta garzetta\nSame Date",saison = c(c(11,12,1,2),c(6,7,8,9)),type=type )

Annee = format(temp_matrice$Date,format = "%Y")
temp_matrice =cbind(temp_matrice,Annee)
temp_matrice$Annee=as.numeric(as.character(temp_matrice$Annee))
temp_matrice$Nom_latin=as.character(temp_matrice$Nom_latin)
temp_matrice$Nombre=as.numeric(as.character(temp_matrice$Nombre))

SynchronySeason2(matrice = temp_matrice,"T49-Axe1-Teich_Synchrony_by_season_for_the_3_competing_species__same_date_for_Phalacrocorax.pdf",titre ="Gross synchrony calculation between \nPhalacrocorax carbo and sum of Ardea cinerea and Egretta garzetta\nSame Dates for Phalacrocorax",Loreau = FALSE,Gross=TRUE,max_value=0.8,min_value=-0.4)

# # with dataset of Coralie juste Heron Egret
# HeronEgret = subset(summed_abundances,summed_abundances$Nom_latin=='HeronEgret')
# temp_data = subset(DBt,(DBt$Nom_latin=='Phalacrocorax carbo' &DBt$Annee>=1981))
# temp_data=subset(temp_data, select=c("Nom_latin", "Date","Nombre"))
# HeronEgret=subset(HeronEgret,select=c("Nom_latin", "Date","Nombre")) # just to be sure
# temp_data=rbind(temp_data,HeronEgret) # for use the same function.
# Annee = format(temp_data$Date,format = "%Y")
# temp_data =cbind(temp_data,Annee)
# temp_data$Annee=as.numeric(as.character(temp_data$Annee))
# temp_data$Nom_latin=as.character(temp_data$Nom_latin)
# temp_data$Nombre=as.numeric(as.character(temp_data$Nombre))
# temp_data$Nombre=as.numeric(as.character(temp_data$Nombre))
# SynchronySeason2(matrice = temp_data,"T49b-Axe1-Teich_Synchrony_by_season_for_the_3_competing_species.pdf",titre ="Gross synchrony calculation between \nPhalacrocorax carbo and sum of Ardea cinerea and Egretta garzetta\n",Loreau = FALSE,Gross=TRUE)
# 
# # with dataset of Coralie juste Heron Egret et cormoran
# temp_data = subset(summed_abundances,summed_abundances$Nom_latin=='HeronEgret' | summed_abundances$Nom_latin=='Cormorant')
# temp_data=subset(temp_data, select=c("Nom_latin", "Date","Nombre"))
# temp_data$Date=as.Date(as.character(temp_data$Date),format="%Y-%m-%d")
# Annee = format(temp_data$Date,format = "%Y")
# temp_data =cbind(temp_data,Annee)
# temp_data$Annee=as.numeric(as.character(temp_data$Annee))
# temp_data$Nom_latin=as.character(temp_data$Nom_latin)
# temp_data$Nombre=as.numeric(as.character(temp_data$Nombre))
# SynchronySeason2(matrice = temp_matrice,"T49b-Axe1-Teich_Synchrony_by_season_for_the_3_competing_species__same_date_for_Phalacrocorax.pdf",titre ="Gross synchrony calculation between \nPhalacrocorax carbo and sum of Ardea cinerea and Egretta garzetta\nSame Dates for Phalacrocorax",Loreau = FALSE,Gross=TRUE)

# Mean abundance, for hot and cold seasons for 2015
#  sum(as.numeric(as.character(temp_matrice$Nombre[temp_matrice$Nom_latin=='Phalacrocorax carbo' & temp_matrice$Annee==2015 & as.numeric(format(temp_matrice$Date, format = "%m")) %in% c(5,6,7,8) ])))/length(temp_matrice$Nombre[temp_matrice$Nom_latin=='Phalacrocorax carbo' & temp_matrice$Annee==2015 & as.numeric(format(temp_matrice$Date, format = "%m")) %in% c(5,6,7,8) ])
# # 
#  sum(as.numeric(as.character(temp_matrice$Nombre[temp_matrice$Nom_latin=='Heron+Egret' & temp_matrice$Annee==2015 & as.numeric(format(temp_matrice$Date, format = "%m")) %in% c(5,6,7,8) ])))/length(temp_matrice$Nombre[temp_matrice$Nom_latin=='Heron+Egret' & temp_matrice$Annee==2015 & as.numeric(format(temp_matrice$Date, format = "%m")) %in% c(5,6,7,8) ])
# # 
# sum(as.numeric(as.character(temp_matrice$Nombre[( temp_matrice$Nom_latin=='Phalacrocorax carbo' & temp_matrice$Annee==2015 & as.numeric(format(temp_matrice$Date, format = "%m")) %in% c(11,12)) |  (temp_matrice$Nom_latin=='Phalacrocorax carbo' & temp_matrice$Annee==2016 & as.numeric(format(temp_matrice$Date, format = "%m")) %in% c(1,2)) ])))/length(temp_matrice$Nombre[( temp_matrice$Nom_latin=='Phalacrocorax carbo' & temp_matrice$Annee==2015 & as.numeric(format(temp_matrice$Date, format = "%m")) %in% c(11,12)) |  (temp_matrice$Nom_latin=='Phalacrocorax carbo' & temp_matrice$Annee==2016 & as.numeric(format(temp_matrice$Date, format = "%m")) %in% c(1,2)) ])
# # 
#  sum(as.numeric(as.character(temp_matrice$Nombre[( temp_matrice$Nom_latin=='Heron+Egret' & temp_matrice$Annee==2015 & as.numeric(format(temp_matrice$Date, format = "%m")) %in% c(11,12)) |  (temp_matrice$Nom_latin=='Heron+Egret' & temp_matrice$Annee==2016 & as.numeric(format(temp_matrice$Date, format = "%m")) %in% c(1,2)) ])))/length(temp_matrice$Nombre[( temp_matrice$Nom_latin=='Heron+Egret' & temp_matrice$Annee==2015 & as.numeric(format(temp_matrice$Date, format = "%m")) %in% c(11,12)) |  (temp_matrice$Nom_latin=='Heron+Egret' & temp_matrice$Annee==2015 & as.numeric(format(temp_matrice$Date, format = "%m")) %in% c(1,2)) ])
# # 
#  t = subset (DBt, (DBt$Nom_latin=='Egretta garzetta' | DBt$Nom_latin=='Ardea cinerea' ) & DBt$Annee==2015 & as.numeric(format(DBt$Date, format = "%m")) %in% c(5,6,7,8) )
#  t<-t[colnames(t) %in% c("Nom_latin","Date","Nombre")]


verification_abondance = function(temp_matrice){
  species = as.character(unique(temp_matrice$Nom_latin))
  min = min(temp_matrice$Annee)
  max = max(temp_matrice$Annee)
  for (s in 1:length(species)){
    print ("Saison froide")
    for (y in min:(max-1)){
      hiver = sum(as.numeric(as.character(temp_matrice$Nombre[( temp_matrice$Nom_latin==species[s] & temp_matrice$Annee==y & as.numeric(format(temp_matrice$Date, format = "%m")) %in% c(11,12)) |  (temp_matrice$Nom_latin==species[s] & temp_matrice$Annee==y+1 & as.numeric(format(temp_matrice$Date, format = "%m")) %in% c(1,2)) ])))/length(temp_matrice$Nombre[( temp_matrice$Nom_latin==species[s] & temp_matrice$Annee==y & as.numeric(format(temp_matrice$Date, format = "%m")) %in% c(11,12)) |  (temp_matrice$Nom_latin==species[s] & temp_matrice$Annee==y+1 & as.numeric(format(temp_matrice$Date, format = "%m")) %in% c(1,2)) ])
      print (paste(y,";",species[s],";",hiver))
    }
    print ("Saison Chaude")
    for (y in min:(max-1)){
      ete = sum(as.numeric(as.character(temp_matrice$Nombre[temp_matrice$Nom_latin==species[s]& temp_matrice$Annee==y & as.numeric(format(temp_matrice$Date, format = "%m")) %in% c(5,6,7,8) ])))/length(temp_matrice$Nombre[temp_matrice$Nom_latin==species[s] & temp_matrice$Annee==y & as.numeric(format(temp_matrice$Date, format = "%m")) %in% c(5,6,7,8) ])  
      print (paste(y,";",species[s],";",ete))
    }
  }
}
verification_abondance(subset(temp_matrice,temp_matrice$Annee>=2007))

###################################################################################################################
# Ratio abundance for each group  relative to total abundance !
# calidris

abondance_calidris = subset(DBt,grepl("^Calidris",DBt$Nom_latin,ignore.case = TRUE) & DBt$Annee>=1981)
print (paste("Ratio abundance Calidris : ",sum(abondance_calidris$Nombre)/sum(DBt$Nombre)*100," %",sep=""))
# anas
abondance_anas = subset(DBt,grepl("^Anas",DBt$Nom_latin,ignore.case = TRUE)& DBt$Annee>=1981)
print (paste("Ratio abundance Anas : ",sum(abondance_anas$Nombre)/sum(DBt$Nombre)*100," %",sep=""))
# waders
abondance_waders = subset(DBt,as.character(DBt$Nom_latin)%in% limicoles & DBt$Annee>=1981)
print (paste("Ratio abundance waders : ",sum(abondance_waders$Nombre)/sum(DBt$Nombre)*100," %",sep=""))
# frequent
abondance_frequent = subset(DBt,as.character(DBt$Nom_latin)%in% oiseaux_Frequents_t & DBt$Annee>=1981)
print (paste("Ratio abundance frequent : ",sum(abondance_frequent$Nombre)/sum(DBt$Nombre)*100," %",sep=""))
# list of sparow (passereaux in french), complete ?
sparrows = c("Alauda arvensis","Lullula arborea","Plectrophenax nivalis","Emberiza schoeniclus","Emberiza citrinella","Emberiza hortulana","Emberiza calandra","Carduelis carduelis","Sturnus vulgaris","Sylvia atricapilla","Sylvia borin","Sylvia communis","Sylvia undata","Garrulus glandarius","Muscicapa striata","Ficedula hypoleuca","Luscinia svecica","Luscinia svecica cyanecula / namnetum","Luscinia svecica cyanecula","Luscinia svecica namnetum","Delichon urbicum","Riparia riparia","Cecropis daurica","Hirundo rustica","Oriolus oriolus","Passer domesticus","Passer montanus","Lanius senator","Lanius collurio","Fringilla coelebs","Fringilla montifringilla","Anthus trivialis","Anthus pratensis","Anthus petrosus","Anthus campestris","Anthus spinoletta","Sitta europaea","Carduelis flammea")
abondance_sparrows = subset(DBt,as.character(DBt$Nom_latin)%in% sparrows & DBt$Annee>=1981)
print (paste("Ratio abundance sparrow : ",sum(abondance_sparrows$Nombre)/sum(DBt$Nombre)*100," %",sep=""))

###################################################################################################################
# synchrony with selected species T50
# calidris
abondance_calidris = subset(DBt,grepl("^Calidris",DBt$Nom_latin,ignore.case = TRUE) & DBt$Annee>=1981)
SynchronySeason2(matrice = abondance_calidris,"T50-Axe1-Teich_Synchrony_by_season_for_Calidris_without_species_with_low_abondance.pdf",titre ="Gross synchrony for Calidris \nwithout_species_with_low_abondance",Loreau = FALSE,Gross=TRUE,max_value=0.8,min_value=-0.4)

abondance_calidris = subset(DBt,DBt$Nom_latin %in% c("Calidris canutus","Calidris alpina", "Calidris ferruginea","Calidris minuta")& DBt$Annee>=1981)
SynchronySeason2(matrice = abondance_calidris,"T50-Axe1-Teich_Synchrony_by_season_for_Calidris_with_only_canutus_alpina_ferruginea_minuta.pdf",titre ="Gross synchrony for Calidris \nwith only canutus Alpina ferruginea et minuta",Loreau = FALSE,Gross=TRUE,max_value=0.8,min_value=-0.4)

#as.Date(dataF_ete$date,origin="1970-01-01") #pour passer du numéric
### coco matrice


#abondance_calidris_coco = subset(calidris_summed_abundances,calidris_summed_abundances$sp_all %in% c("Calidris canutus","Calidris alpina", "Calidris ferruginea","Calidris minuta"))
#colnames(abondance_calidris_coco)=c("X","Annee","Nom_latin","Nombre")
#Annee = format(abondance_calidris_coco$Date,format = "%Y")
#abondance_calidris_coco =cbind(abondance_calidris_coco,Annee)

# abondance_calidris_coco$Annee=as.numeric(as.character(abondance_calidris_coco$Annee))
# abondance_calidris_coco$Nom_latin=as.character(abondance_calidris_coco$Nom_latin)marcel-plenacoste-la-legende-des-filles-de-la-mer-d-irlande
# abondance_calidris_coco$Nombre=as.numeric(as.character(abondance_calidris_coco$Nombre))
#SynchronySeason2(matrice = abondance_calidris_coco,"T50-Axe1-Teich_Synchrony_by_season_for_Calidris_with_only_canutus_alpina_ferruginea_minuta_cocoVersion.pdf",titre ="Gross synchrony for Calidris \nwith only canutus Alpina ferruginea et minuta\nCoco version",Loreau = FALSE,Gross=TRUE)

# anas
abondance_anas = subset(DBt,grepl("^Anas",DBt$Nom_latin,ignore.case = TRUE)& DBt$Annee>=1981)
SynchronySeason2(matrice = abondance_anas,"T50-Axe1-Teich_Synchrony_by_season_for_anas_without_species_with_low_abondance.pdf",titre ="Gross synchrony for Anas\nwithout_species_with_low_abondance",Loreau = FALSE,Gross=TRUE,max_value=0.8,min_value=-0.4)

print ("Abondance waders")
# waders
abondance_waders = subset(DBt,as.character(DBt$Nom_latin)%in% limicoles & DBt$Annee>=1981)
SynchronySeason2(matrice = abondance_waders,"T50-Axe1-Teich_Synchrony_by_season_for_waders_without_species_with_low_abondance.pdf",titre ="Gross synchrony for Waders\nwithout_species_with_low_abondance",Loreau = FALSE,Gross=TRUE,max_value=0.8,min_value=-0.4)

print ("Abondance sparrows")
#sparrows
abondance_sparrows = subset(DBt,as.character(DBt$Nom_latin)%in% sparrows & DBt$Annee>=1981)
SynchronySeason2(matrice = abondance_sparrows,file_out ="T50-Axe1-Teich_Synchrony_by_season_for_sparrows_without_species_with_low_abondance.pdf",titre ="Gross synchrony for Sparrows",Loreau = FALSE,Gross=TRUE,max_value=0.8,min_value=-0.4)
