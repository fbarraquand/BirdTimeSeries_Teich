####CP 06/09/2019: This is a cleaned up version of analyseP_oiseaux_Teich2_synchrony.R (entirely written by CAluome) to keep only what we need and use in the article (the rest might be put in exploratory analyses) 

rm(list=ls())
graphics.off()
library("lattice")
library ("RColorBrewer") #pour la génération automatique de couleur
library("corrplot")
library("lubridate")
library("stringr")
library("codyn")
library ("pracma")
DIRECTORY_ORIGIN = "./"
#setwd(paste(DIRECTORY_ORIGIN,"IN/",sep=""))
#options(nwarnings = 300) #nb de messages de warnings conservés

#update.packages(ask = F) #update all packages

##################################################################################
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
#summed_abundances<-read.csv(file="summed_abundances.csv",header=TRUE,sep=";",dec=".") 
#calidris_summed_abundances<-read.csv(file="coldseason_calidris_detailed.txt",header=TRUE,sep=",",dec=".") 
#setwd(paste(DIRECTORY_ORIGIN,"OUT/",sep=""))

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

filename_pdf="abundance_not_averaged.pdf"
pdf(paste("OUT/",filename_pdf,sep=""),width=15,height=10)
# rather make three graphs!!!
par(mfrow=c(3,2))
layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = TRUE),heights = c(1,1,1),widths = c(3,1))
abondance_anas = subset(DBt,grepl("^Anas",DBt$Nom_latin,ignore.case = TRUE) & DBt$Annee>=1981)
abondance_calidris = subset(DBt,grepl("^Calidris|Philomachus pugnax",DBt$Nom_latin,ignore.case = TRUE)& DBt$Annee>=1981)
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
ymax = log10(max(abondance_anas$Nombre))+0.5 #tocorrect
xmin = date("1981-01-01")
xmax = date("2015-01-01")
par(mar=c(1.5,4.3,0.1,1))
par(oma=c(0,0,0,0))

plot(temp_date,log10(temp_value),col="black",xlab="",ylab="Log10(abundance)",xlim=c(xmin,xmax),cex.lab=2,cex.axis=1.8,xaxt="n",pch=19,ylim=c(ymin,ymax),type="l",lwd =2,las=1)
lines(temp_date,log10(temp_value+1),col="black",lwd =1,lty=3)
for (s in 1:length(temp_species)){
  points(abondance_anas$Date[as.character(abondance_anas$Nom_latin) == temp_species[s]],log10(abondance_anas$Nombre[abondance_anas$Nom_latin==temp_species[s]]),col=listColor[s],pch=19)
}
#text(c(0,1), "a)",cex = 2)
mtext("\n (a)",side=3,cex = 1.5, adj = 0,padj=1)
axis(1, at=d_ticks, cex=1.8, ticks=TRUE, labels=FALSE) # pas de label

par(mar=c(0.1,0.1,0.1,0.1))
mat_percent_abondance_tot_anas = mat_percent_abondance_tot_anas[ order(row.names(mat_percent_abondance_tot_anas)), ]

id_1=which(rownames(mat_percent_abondance_tot_anas)=="Anas clypeata")
rownames(mat_percent_abondance_tot_anas)[id_1]="Spatula clypeata"
id_1=which(rownames(mat_percent_abondance_tot_anas)=="Anas penelope")
rownames(mat_percent_abondance_tot_anas)[id_1]="Mareca penelope"
id_1=which(rownames(mat_percent_abondance_tot_anas)=="Anas querquedula")
rownames(mat_percent_abondance_tot_anas)[id_1]="Spatula querquedula"
id_1=which(rownames(mat_percent_abondance_tot_anas)=="Anas strepera")
rownames(mat_percent_abondance_tot_anas)[id_1]="Mareca strepera"
doughnut( c(mat_percent_abondance_tot_anas[,'percent'][round(mat_percent_abondance_tot_anas[,'percent'])>=minimal_percent]) , inner.radius=0.15, col=listColor[which(round(mat_percent_abondance_tot_anas[,'percent'])>=minimal_percent)],outer.radius =0.45,label.cex=1.3,coord.x=c(-0.5,1),init.angle=50)
mtext(paste(nrow(mat_percent_abondance_tot_anas)," species",sep=""),1,line=-5)

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
ymax = log10(max(abondance_calidris$Nombre))

plot(temp_date,log10(temp_value),col="black",xlab="",ylab="Log10(abundance)",xlim=c(xmin,xmax),cex.axis=1.8,cex.lab=2,xaxt="n",pch=19,ylim=c(ymin,ymax),type="l",lwd =2,las=1)
lines(temp_date,log10(temp_value+1),col="black",lwd =1,lty=3)
for (s in 1:length(temp_species)){
  points(abondance_calidris$Date[as.character(abondance_calidris$Nom_latin) == temp_species[s]],log10(abondance_calidris$Nombre[abondance_calidris$Nom_latin==temp_species[s]]),col=listColor[s],pch=19)
}
axis(1, at=d_ticks, cex=1.8, ticks=TRUE, labels=FALSE) # pas de label
mtext("\n (b)",side=3,cex = 1.5, adj = 0,padj = 1)
# plot.new()
par(mar=c(0.1,0.1,0.1,0.1))
mat_percent_abondance_tot_calidris = mat_percent_abondance_tot_calidris[ order(row.names(mat_percent_abondance_tot_calidris)), ]


doughnut( c(mat_percent_abondance_tot_calidris[,'percent'][round(mat_percent_abondance_tot_calidris[,'percent'])>=minimal_percent]) , inner.radius=0.145, col=listColor[which(round(mat_percent_abondance_tot_calidris[,'percent'])>=minimal_percent)],outer.radius =0.45,label.cex=1.28,coord.x=c(-0.7,1) ,init.angle=25)
mtext(paste(nrow(mat_percent_abondance_tot_calidris)," species",sep=""),1,line=-5)

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
ymax = log10(max(abondance_waders$Nombre))
plot(temp_date,log10(temp_value),col="black",xlab="Year",ylab="Log10(abundance)",xlim=c(xmin,xmax),cex.axis=1.8,cex.lab=2,xaxt="n",pch=19,ylim=c(ymin,ymax),type="l",lwd =2,las=1)
lines(temp_date,log10(temp_value+1),col="black",lwd =1,lty=3)
for (s in 1:length(temp_species)){
  points(abondance_waders$Date[as.character(abondance_waders$Nom_latin) == temp_species[s]],log10(abondance_waders$Nombre[abondance_waders$Nom_latin==temp_species[s]]),col=listColor[s],pch=19)
}
mtext("\n (c)",side=3,cex = 1.5, adj = 0,padj = 1)
axis(1, at=d_ticks, labels=nameTicks,cex=1.8,cex.axis=1.8)
mat_percent_abondance_tot_waders = mat_percent_abondance_tot_waders[ order(row.names(mat_percent_abondance_tot_waders)), ]
par(mar=c(0.1,0.1,0.1,0.1))

doughnut( c(mat_percent_abondance_tot_waders[,'percent'][round(mat_percent_abondance_tot_waders[,'percent'])>=minimal_percent]) , inner.radius=0.15, col=listColor[which(round(mat_percent_abondance_tot_waders[,'percent'])>=minimal_percent)],outer.radius =0.45,label.cex=1.28,coord.x=c(-0.5,1),lty=1,init.angle=-10)
mtext(paste(nrow(mat_percent_abondance_tot_waders)," species",sep=""),1,line=-5)

dev.off()
system(paste("cp OUT/",filename_pdf," Submission_JAE/Revisions/",filename_pdf,sep=""))

#pdf("T47-Axe1-Teich-les_3_especes.pdf",width=20,height=6,family = 'Times')
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

filename_pdf="abundance_3species.pdf"
pdf(paste("OUT/",filename_pdf,sep=""),width=20,height=6)
#pdf("T47B-Axe1-Teich-les_3_especes_Heron+Aigrette.pdf",width=20,height=6)#,family = 'Times') #CP changed this
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
ymax = log10(max(df_HeronEgret$Nombre,abondance_phalacrocorax$Nombre,na.rm=TRUE))
xmin = date("1981-01-01")
xmax = date("2015-01-01")
plot(df_HeronEgret$Date,log10(df_HeronEgret$Nombre),col="grey",xlab="Year",ylab="Log10(abundance)",xlim=c(xmin,xmax),cex.lab=1.8,cex.axis=1.5,xaxt="n",pch=19,ylim=c(ymin,ymax),type="o",cex=0.9,las=1,lwd=2)
lines(df_HeronEgret$Date,log10(df_HeronEgret$Nombre+1),col="grey",lwd =1,lty=3)
lines(abondance_phalacrocorax$Date,log10(abondance_phalacrocorax$Nombre),col="black",type='o',lwd=2,pch=19,cex=0.9,lty=2)
lines(abondance_phalacrocorax$Date,log10(abondance_phalacrocorax$Nombre+1),col="black",lwd =1,lty=3)

axis(1, at=d_ticks, labels=nameTicks,cex.axis=1.5)
legend("topright",legend = c("Grey heron + little egret",'Great cormorant'),col=c('grey','black'),pch=19,lwd=1,cex=1.5,bty="n",lty=c(1,2))
dev.off()
system(paste("cp OUT/",filename_pdf," Submission_JAE/Revisions/",filename_pdf,sep=""))

###################################################################################################################
# Ratio abundance for each group  relative to total abundance !
# calidris

abondance_calidris = subset(DBt,grepl("^Calidris|Philomachus pugnax",DBt$Nom_latin,ignore.case = TRUE) & DBt$Annee>=1981)
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

#### T30
minAnnee = 1981
vec_abondance_waders      = rep(0,(maxAnnee-minAnnee+1)) # abondances des limicoles
vec_abondance_ducks       = rep(0,(maxAnnee-minAnnee+1)) # abondances des limicoles
vec_abondance_all_others  = rep(0,(maxAnnee-minAnnee+1)) # abondances des autres oiseaux.
vec_abondance_all         = rep(0,(maxAnnee-minAnnee+1)) # abondances de tous les oiseaux
vec_ratio_limicoles       = rep(0,(maxAnnee-minAnnee+1)) # vecteur qui stocke l'information ratio entre les limicoles et les non limicoles.

#Define Ducks : All Anatidae + Fulica atra
tab_tmp=read.csv(file="IN/Initial_files/data_ROT20160324.csv",header=TRUE,sep="\t")
sp=tab_tmp$Nom_latin
fam=tab_tmp$Famille

ducks=c(as.character(unique(sp[fam=="Anatidae"])),"Fulica atra") #Other Rallidae do not dive to feed, I considered this feature as a definition of ducks


# on fait les calculs pour toutes les années

for (y in minAnnee:maxAnnee){
  vec_abondance_all[y-minAnnee+1] = sum(DBt$Nombre[DBt$Annee==y],na.rm = TRUE)
  vec_abondance_waders[y-minAnnee+1] = sum(DBt$Nombre[DBt$Annee == y & DBt$Nom_latin %in% limicoles])
  vec_abondance_ducks[y-minAnnee+1] = sum(DBt$Nombre[DBt$Annee == y & DBt$Nom_latin %in% ducks])
  vec_abondance_all_others[y-minAnnee+1]=vec_abondance_all[y-minAnnee+1]-( vec_abondance_waders[y-minAnnee+1] +vec_abondance_ducks[y-minAnnee+1]  )
  vec_ratio_limicoles[y-minAnnee+1]  = vec_abondance_waders[y-minAnnee+1]/vec_abondance_all[y-minAnnee+1]
}

# on superpose les 3 courbes sur le graphe : abondance des limicoles, des non-limicoles et ratio.
filename_pdf="ProportionsLimicoles.pdf"
pdf(paste("OUT/",filename_pdf,sep=""),width=12,height=8)
par(mar=c(4,4,3,3))
par(oma = c(3,1,0,1))
plot(minAnnee:maxAnnee,vec_ratio_limicoles,type="o",xlab="Year",ylab="Ratio",
     col="black",lwd=2,cex=1.5,pch=15,lty=2,main="",las=1) #main =evolution of abundance of shorebirds, ducks and other species
abline(v=2006,col="black",lwd=2,lty=3)
par(new=TRUE)
plot(minAnnee:maxAnnee,log10(vec_abondance_all),type="o",
     xlab="Year",ylab="",xaxt="n",yaxt="n",pch=19,col="deepskyblue",ylim=c(1,log10(max(vec_abondance_all))))
lines(minAnnee:maxAnnee,log10(vec_abondance_waders),type="o",col="darkorchid",pch=19)
lines(minAnnee:maxAnnee,log10(vec_abondance_ducks),type="o",col="olivedrab",pch=19)
lines(minAnnee:maxAnnee,log10(vec_abondance_all_others),type="o",col="salmon",pch=19)
axis(4,col="black",las=1)
mtext("Log10(abundance)",side=4,outer = TRUE)
#légende
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom",c("Ratio waders",'All birds','Waders',"Waterfowl",'All other birds'),
       col=c("black","deepskyblue","darkorchid",'olivedrab','salmon'),xpd = TRUE, 
       ncol=5 , bty = "n",cex=1,pch=c(15,19,19,19,19),lty=c(2,1,1,1,1),lwd=c(2,1,1,1,1))
dev.off()
system(paste("cp OUT/",filename_pdf," Submission_JAE/Revisions/",filename_pdf,sep=""))

#### Just checking Cormorant and HeronEgret when they are averaged
w1=read.table("IN/warmseason_abundances_summed_v2_wtoutrarespecies.csv",header=T,sep=";")
c2=read.table("IN/coldseason_abundances_summed_v2_wtoutrarespecies.csv",header=T,sep=";")

filename_pdf="averaged_CorHerEgr_log.pdf"
pdf(paste("OUT/",filename_pdf,sep=""),width=10)
par(mfrow=c(2,2),mar=c(3,4,2,1))
plot(minAnnee:maxAnnee,w1$HeronEgret,col="grey",t="l",xlab="Year",ylab="Abundance",lwd=2,main="Warm",ylim=c(0,max(w1$HeronEgret,na.rm=T)))
lines(minAnnee:maxAnnee,w1$Cormorant,col="black",lty=2,lwd=2)
plot(minAnnee:maxAnnee,c2$Cormorant,col="black",lty=2,t="l",ylim=c(0,max(c2$HeronEgret,na.rm=T)),ylab="",xlab="Year",lwd=2,main="Cold")
lines(minAnnee:maxAnnee,c2$HeronEgret,col="grey",lwd=2)

plot(minAnnee:maxAnnee,log10(w1$HeronEgret),col="grey",t="l",xlab="Year",ylab="log10(abundance)",lwd=2,main="",ylim=c(0,log10(max(w1$HeronEgret,na.rm=T))))
lines(minAnnee:maxAnnee,log10(w1$Cormorant),col="black",lty=2,lwd=2)
plot(minAnnee:maxAnnee,log10(c2$Cormorant),col="black",lty=2,t="l",ylim=c(0,log10(max(c2$HeronEgret,na.rm=T))),ylab="",xlab="Year",lwd=2,main="")
lines(minAnnee:maxAnnee,log10(c2$HeronEgret),col="grey",lwd=2)
legend('bottomright',c('Grey heron + little egret','Great cormorant'),pch=16,lty=c(1,2),col=c("grey","black"),bty="n")
dev.off()
system(paste("cp OUT/",filename_pdf," Submission_JAE/Revisions/",filename_pdf,sep=""))
