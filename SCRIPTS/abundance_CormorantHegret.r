rm(list=ls())
graphics.off()

library("lubridate")
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


# -- importation des données du Teich
DBt<-read.csv(file="IN/DBWithMonthlyPhotoTeich_completed.csv",header=TRUE,sep=",",dec=".")

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
              "Tringa nebularia","Tringa erythropus","Tringa ochropus","Tringa totanus",
              "Tringa glareola","Actitis hypoleucos","Philomachus pugnax","Numenius arquata","Numenius phaeopus",
              "Himantopus himantopus","Charadrius hiaticula","Charadrius alexandrinus","Haematopus ostralegus",
              "Burhinus oedicnemus","Charadrius dubius","Phalaropus lobatus","Pluvialis squatarola",
              "Pluvialis apricaria","Arenaria interpres","Vanellus vanellus","Calidris ferruginea")
#Tringa flavipes never appers with Protocol>0

birds_to_remove=c("Anas discors","Anas americana","Calidris melanotos","Calidris pusilla","Calidris ruficollis", "Calidris fuscicollis", "Calidris himantopus", "Burhinus oedicnemus","Phalaropus lobatus","Charadrius alexandrinus","Haematopus ostralegus","Calidris maritima","Aythya nyroca","Bucephala clangula","Melanitta nigra","Mergus serrator","Clangula hyemalis","Alopochen aegyptiaca", "Aix galericulata","Cygnus atratus","Tadorna ferruginea","Branta leucopsis","Anser fabalis","Anser albifrons","Cygnus cygnus","Mergus merganser","Anser brachyrhynchus")

DBt$Nom_latin=as.character(DBt$Nom_latin)
DBt=subset(DBt,!DBt$Nom_latin %in% birds_to_remove)


filename_pdf="abundance_3species_nolog.pdf"
pdf(paste("OUT/",filename_pdf,sep=""),width=20,height=6)
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
par(mar=c(4,6,3,1))
d_ticks = date(paste(unique(as.numeric(format(unique(DBt$Date), format = "%Y"))),"-01-01",sep=""))
d_ticks = d_ticks[seq(1,length(d_ticks),by = 2)]
nameTicks=unique(as.numeric(format(unique(DBt$Date), format = "%Y")))
nameTicks = nameTicks[seq(1,length(nameTicks),by = 2)]
abondance_phalacrocorax = subset(DBt,DBt$Nom_latin=='Phalacrocorax carbo' & DBt$Annee>=1981)
ymin = 0
ymax = max(df_HeronEgret$Nombre,abondance_phalacrocorax$Nombre,na.rm=TRUE)
#ymax = log10(max(df_HeronEgret$Nombre,abondance_phalacrocorax$Nombre,na.rm=TRUE))
xmin = date("1981-01-01")
xmax = date("2015-01-01")
plot(df_HeronEgret$Date,df_HeronEgret$Nombre,col="grey",xlab="Year",ylab="",xlim=c(xmin,xmax),cex.lab=1.8,cex.axis=1.5,xaxt="n",pch=19,ylim=c(ymin,ymax),type="o",cex=0.9,las=1,lwd=2)
lines(abondance_phalacrocorax$Date,abondance_phalacrocorax$Nombre,col="black",type='o',lwd=2,pch=19,cex=0.9,lty=2)
mtext("abundance",2,line=4.2,cex=1.8)
axis(1, at=d_ticks, labels=nameTicks,cex.axis=1.5)
legend("topright",legend = c("Grey heron + little egret",'Great cormorant'),col=c('grey','black'),pch=19,lwd=1,cex=1.5,bty="n",lty=c(1,2))
dev.off()

