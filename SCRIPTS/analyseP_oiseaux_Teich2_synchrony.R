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
setwd(paste(DIRECTORY_ORIGIN,"IN/",sep=""))
options(nwarnings = 300) #nb de messages de warnings conservés



##################################################################################
# -- importation des données du Teich
#DBt<-read.csv(file="/home/caluome/Documents/DATA/DATA/le_Teich/DBWithMonthlyPhotoTeich.csv",header=TRUE,sep=",",dec=".")
#DBt<-read.csv(file="DBWithMonthlyPhotoTeich.csv",header=TRUE,sep=",",dec=".")
DBt<-read.csv(file="DBWithMonthlyPhotoTeich_maxsum.csv",header=TRUE,sep=",",dec=".")

DBt = subset(DBt,(DBt$Protocol==1 | DBt$Protocol==2) & DBt$Lieu_dit=="USN00-Réserve ornithologique (générique)")
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
              "Burrhinus oedicnemus","Charadrius dubius","Phalaropus lobatus","Pluvialis squatarola",
              "Pluvialis apricaria","Arenaria interpres","Vanellus vanellus")





setwd(paste(DIRECTORY_ORIGIN,"OUT/",sep=""))

# ---------------------------------------------------------------------
#
#                              SCRIPT T-36
#
# ---------------------------------------------------------------------
#### T36 : SYNCHRONY CALCUL by MONTH
SynchronyMonth=function(matrice,file_out,titre){ 
  pdf(file_out,width=14,height=8)
  par(mar=c(4,4.5,3,2.5))
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
  plot(1:12,vec_synchrony_Loreau,type="o",ylim=c(-0.05,1),pch=19, 
       main=titre,ylab="Value of synchrony",xlab="month",cex.main=2,cex.lab=2,
       cex=2,cex.axis=2)
  lines(1:12,vec_synchrony_Gross,col="orange",type = "o",pch=19,cex=2)
  #-----------------------
  # PERIODE >= 2006
	vec_synchrony_Loreau = rep(0,12) #Final results
  vec_synchrony_Gross = rep(0,12) #Final results
  for (m in 1:12){
    subdata   = subset(DBt,as.numeric(format(DBt$Date, format = "%m"))==m & as.numeric(format(DBt$Date, format = "%Y"))>=2006 & DBt$Nom_latin %in% limicoles)
    date      = subdata$Date
    date      = c(as.numeric(date))
    abundance = c(as.numeric(subdata$Nombre))
    species   = c(as.character(subdata$Nom_latin))
    dataF_DBt = data.frame(date,abundance,species)#Codyn need dataframe to work
    colnames(dataF_DBt)=c("date","abundance","species") 
    vec_synchrony_Loreau[m]= synchrony(dataF_DBt, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Loreau")
    vec_synchrony_Gross[m]=synchrony(dataF_DBt, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Gross") 
    
  }
  lines(1:12,vec_synchrony_Loreau,col="blue",type = "o",pch=19,cex=2)
  lines(1:12,vec_synchrony_Gross,col="red",type = "o",pch=19,cex=2)
  #-----------------------
  # PERIODE < 2006
  vec_synchrony_Loreau = rep(0,12) #Final results
  vec_synchrony_Gross = rep(0,12) #Final results
  for (m in 1:12){
    subdata   = subset(DBt,as.numeric(format(DBt$Date, format = "%m"))==m & as.numeric(format(DBt$Date, format = "%Y"))<2006 & DBt$Nom_latin %in% limicoles)
    date      = subdata$Date
    date      = c(as.numeric(date))
    abundance = c(as.numeric(subdata$Nombre))
    species   = c(as.character(subdata$Nom_latin))
    dataF_DBt = data.frame(date,abundance,species) #Codyn need dataframe to work
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
         inset = c(0,0),bty = "n",cex=0.7,pt.cex=1.5,lty=1,pt.lwd=1.5,lwd=1.5)
  
  dev.off()
}

## - Calls of the function n°T36
SynchronyMonth(DBt,"t36-Axe1-Teich-SynchronyComparison_byMonth_AllSpecies.pdf","Evolution of Synchrony from Loreau and Gross by month (all species)")
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
SynchronySeason(DBt,"t40-Axe1-Teich-Comparison_synchrony_by_season_all_species.pdf","Comparison Synchrony Season - All species",saison,type)
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
SynchronySeason2 = function(matrice,file_out,titre){
  saison=list(c(11,12,1,2,3),c(5,6,7,8,9))
  saison =  list(hiver_n=c(11, 12), hiver_n1=c(1,2,3),ete=c(5,6,7,8,9))
  
  vec_synchrony_Loreau=rep(0,6) #6 for 2season by 3 periods (all,before 2006 and from 2006)
  vec_synchrony_Gross=rep(0,6)
  abondance_ete=c()
  date_ete=c()
  species_ete=c()
  abondance_hiver=c()
  date_hiver=c()
  species_hiver=c()
  for(y in min(subdata$Annee):max(subdata$Annee)-1){
      subdata_ete = subset(matrice,as.numeric(format(matrice$Date, format = "%m")) %in% saison[['ete']] & as.numeric(format(matrice$Date, format = "%Y"))==y)
      subdata_hiver = subset(matrice,(as.numeric(format(matrice$Date, format = "%m")) %in% saison[['hiver_n']] & as.numeric(format(matrice$Date, format = "%Y"))==y) | (as.numeric(format(matrice$Date, format = "%m")) %in% saison[['hiver_n1']] & as.numeric(format(matrice$Date, format = "%Y"))==y+1))
      
      if(dim(subdata_ete)[1]>0){
        species_en_cours = c(unique(as.character(subdata_ete$Nom_latin)))
        for (e in 1:length(species_en_cours)){
          print (paste("sepcies",species_en_cours[e]))
          print (paste("mean",mean(subdata_ete$Nombre[subdata_ete$Nom_latin==species_en_cours[e]])))
          abondance_ete=c(abondance_ete,c(mean(subdata_ete$Nombre[subdata_ete$Nom_latin==species_en_cours[e]])))
          date_ete=c(date_ete,c(as.Date(as.character(subdata_ete$Date[1]),format="%Y-%m-%d"))) 
          species_ete =c(species_ete,c(species_en_cours[e]))
        }
      }
      if(dim(subdata_hiver)[1]>0){
        species_en_cours = c(unique(as.character(subdata_hiver$Nom_latin)))
        for (e in 1:length(species_en_cours)){
          print (paste("sepcies",species_en_cours[e]))
          print (paste("mean",mean(subdata_hiver$Nombre[subdata_hiver$Nom_latin==species_en_cours[e]])))
          abondance_hiver=c(abondance_hiver,c(mean(subdata_hiver$Nombre[subdata_hiver$Nom_latin==species_en_cours[e]])))
          date_hiver=c(date_hiver,c(as.Date(as.character(subdata_hiver$Date[1]),format="%Y-%m-%d"))) 
          species_hiver =c(species_hiver,c(species_en_cours[e]))
        }
      }
  }
  date      = c(as.numeric(date_ete))
  abundance = c(as.numeric(abondance_ete))
  species   = c(as.character(species_ete))
  dataF_ete = data.frame(date,abundance,species) #Codyn need dataframe
  colnames(dataF_ete)=c("date","abundance","species") 
  vec_synchrony_Loreau[1]= synchrony(dataF_ete, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Loreau")
  vec_synchrony_Gross[1]=synchrony(dataF_ete, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Gross")
  date      = c(as.numeric(date_hiver))
  abundance = c(as.numeric(abondance_hiver))
  species   = c(as.character(species_hiver))
  dataF_hiver = data.frame(date,abundance,species) 
  colnames(dataF_hiver)=c("date","abundance","species") 
  vec_synchrony_Loreau[2]= synchrony(dataF_hiver, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Loreau")
  vec_synchrony_Gross[2]=synchrony(dataF_hiver, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Gross")
  
  #-----------------------
  # PERIODE >= 2006

  for(y in 2006:max(subdata$Annee)-1){
    subdata_ete = subset(matrice,as.numeric(format(matrice$Date, format = "%m")) %in% saison[['ete']] & as.numeric(format(matrice$Date, format = "%Y"))==y)
    subdata_hiver = subset(matrice,(as.numeric(format(matrice$Date, format = "%m")) %in% saison[['hiver_n']] & as.numeric(format(matrice$Date, format = "%Y"))==y) | (as.numeric(format(matrice$Date, format = "%m")) %in% saison[['hiver_n1']] & as.numeric(format(matrice$Date, format = "%Y"))==y+1))
    
    if(dim(subdata_ete)[1]>0){
      species_en_cours = c(unique(as.character(subdata_ete$Nom_latin)))
      for (e in 1:length(species_en_cours)){
        print (paste("sepcies",species_en_cours[e]))
        print (paste("mean",mean(subdata_ete$Nombre[subdata_ete$Nom_latin==species_en_cours[e]])))
        abondance_ete=c(abondance_ete,c(mean(subdata_ete$Nombre[subdata_ete$Nom_latin==species_en_cours[e]])))
        date_ete=c(date_ete,c(as.Date(as.character(subdata_ete$Date[1]),format="%Y-%m-%d"))) 
        species_ete =c(species_ete,c(species_en_cours[e]))
      }
    }
    if(dim(subdata_hiver)[1]>0){
      species_en_cours = c(unique(as.character(subdata_hiver$Nom_latin)))
      for (e in 1:length(species_en_cours)){
        print (paste("sepcies",species_en_cours[e]))
        print (paste("mean",mean(subdata_hiver$Nombre[subdata_hiver$Nom_latin==species_en_cours[e]])))
        abondance_hiver=c(abondance_hiver,c(mean(subdata_hiver$Nombre[subdata_hiver$Nom_latin==species_en_cours[e]])))
        date_hiver=c(date_hiver,c(as.Date(as.character(subdata_hiver$Date[1]),format="%Y-%m-%d"))) 
        species_hiver =c(species_hiver,c(species_en_cours[e]))
      }
    }
  }
  date      = c(as.numeric(date_ete))
  abundance = c(as.numeric(abondance_ete))
  species   = c(as.character(species_ete))
  dataF_ete = data.frame(date,abundance,species) #Codyn need dataframe
  colnames(dataF_ete)=c("date","abundance","species") 
  vec_synchrony_Loreau[3]= synchrony(dataF_ete, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Loreau")
  vec_synchrony_Gross[3]=synchrony(dataF_ete, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Gross")
  date      = c(as.numeric(date_hiver))
  abundance = c(as.numeric(abondance_hiver))
  species   = c(as.character(species_hiver))
  dataF_hiver = data.frame(date,abundance,species)
  colnames(dataF_hiver)=c("date","abundance","species") 
  vec_synchrony_Loreau[4]= synchrony(dataF_hiver, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Loreau")
  vec_synchrony_Gross[4]=synchrony(dataF_hiver, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Gross")
  #-----------------------
  # PERIODE < 2006
  for(y in min(subdata$Annee):2006){
    subdata_ete = subset(matrice,as.numeric(format(matrice$Date, format = "%m")) %in% saison[['ete']] & as.numeric(format(matrice$Date, format = "%Y"))==y)
    subdata_hiver = subset(matrice,(as.numeric(format(matrice$Date, format = "%m")) %in% saison[['hiver_n']] & as.numeric(format(matrice$Date, format = "%Y"))==y) | (as.numeric(format(matrice$Date, format = "%m")) %in% saison[['hiver_n1']] & as.numeric(format(matrice$Date, format = "%Y"))==y+1))
    
    if(dim(subdata_ete)[1]>0){
      species_en_cours = c(unique(as.character(subdata_ete$Nom_latin)))
      for (e in 1:length(species_en_cours)){
        print (paste("sepcies",species_en_cours[e]))
        print (paste("mean",mean(subdata_ete$Nombre[subdata_ete$Nom_latin==species_en_cours[e]])))
        abondance_ete=c(abondance_ete,c(mean(subdata_ete$Nombre[subdata_ete$Nom_latin==species_en_cours[e]])))
        date_ete=c(date_ete,c(as.Date(as.character(subdata_ete$Date[1]),format="%Y-%m-%d"))) 
        species_ete =c(species_ete,c(species_en_cours[e]))
      }
    }
    if(dim(subdata_hiver)[1]>0){
      species_en_cours = c(unique(as.character(subdata_hiver$Nom_latin)))
      for (e in 1:length(species_en_cours)){
        print (paste("sepcies",species_en_cours[e]))
        print (paste("mean",mean(subdata_hiver$Nombre[subdata_hiver$Nom_latin==species_en_cours[e]])))
        abondance_hiver=c(abondance_hiver,c(mean(subdata_hiver$Nombre[subdata_hiver$Nom_latin==species_en_cours[e]])))
        date_hiver=c(date_hiver,c(as.Date(as.character(subdata_hiver$Date[1]),format="%Y-%m-%d"))) 
        species_hiver =c(species_hiver,c(species_en_cours[e]))
      }
    }
  }
  date      = c(as.numeric(date_ete))
  abundance = c(as.numeric(abondance_ete))
  species   = c(as.character(species_ete))
  dataF_ete = data.frame(date,abundance,species) #Codyn need dataframe
  colnames(dataF_ete)=c("date","abundance","species") 
  vec_synchrony_Loreau[5]= synchrony(dataF_ete, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Loreau")
  vec_synchrony_Gross[5]=synchrony(dataF_ete, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Gross")
  date      = c(as.numeric(date_hiver))
  abundance = c(as.numeric(abondance_hiver))
  species   = c(as.character(species_hiver))
  dataF_hiver = data.frame(date,abundance,species)
  colnames(dataF_hiver)=c("date","abundance","species") 
  vec_synchrony_Loreau[6]= synchrony(dataF_hiver, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Loreau")
  vec_synchrony_Gross[6]=synchrony(dataF_hiver, abundance.var = "abundance",species.var = "species",time.var = "date",metric="Gross")
  

  data_synchro =c(c(vec_synchrony_Loreau),c(vec_synchrony_Gross))
  max_value = max(data_synchro,na.rm=TRUE)+0.002
  data_synchro = matrix(data_synchro,nc=6, byrow=T)
  type =c("winter-all","summer-all","winter>=2006","summer>=2006","winter<2006","summer<2006")
  colnames(data_synchro) = type
  pdf(file_out,width=14,height=8)
  par(mar=c(4,4.5,3,2.5))
  par(oma = c(4.2, 0.5, 0.5, 0.5))
  barplot(data_synchro, xlab="X",ylab="Y", beside=T, col=c("#FFFFAA","#AAFFAA"), ylim=c(0,max_value), main=titre,cex.lab=1,cex.axis = 1,cex.main=2) # Tracer les barres
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom",c("Synchony by Loreau","Synchrony by Gross"), col=c("#FFFFAA","#AAFFAA"),pch=c(16,16), xpd = TRUE, 
         horiz=TRUE,
         inset = c(0,0),bty = "n",cex=0.8,pt.cex=1.5,lty=1,pt.lwd=1.5,lwd=1.5)
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

plot(abondance_calidris$Date,log10(abondance_calidris$Nombre),col="blue",xlab="Date",ylab="Abondance of all Calidris",xlim=c(minD,maxD),main='Log of abundance of Calidris in the dataset',cex.lab=1.2,cex.main=1.5,xaxt="n",pch=19)
lines(abondance_calidris$Date,movavg(log10(abondance_calidris$Nombre),3,type='t'),col="orange",lwd=1.5)
axis(1, at=ticks, labels=nameTicks)


plot(abondance_anas$Date,log10(abondance_anas$Nombre),col="green",xlab="Date",ylab="Abondance of all Anas",xlim=c(minD,maxD),main='Log of abundance of Anas in the dataset',cex.lab=1.2,cex.main=1.5,xaxt="n",pch=19)
lines(abondance_anas$Date,movavg(log10(abondance_anas$Nombre),3,type='t'),col="orange",lwd=1.5)
axis(1, at=ticks, labels=nameTicks)


dev.off()
