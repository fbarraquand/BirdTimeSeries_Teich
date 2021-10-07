###CP Just extracting part of the code of analyseP_oiseaux... entirely written by CAluome, just for efficiency
###CP 04/09/2019: Adding las=1 to have horizontal y-axis tick labels and updating taxonomy
###CP 07/03/2020: Modified shapes of the points to make the figure easier to read, that's why there are "v2" versions

rm(list=ls())
graphics.off()

#First, define functions (these functions were written by CA)

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

#Load Database
DBt<-read.csv(file="IN/DBWithMonthlyPhotoTeich_completed.csv",header=TRUE,sep=",",dec=".")
DBt = subset(DBt,(((DBt$Protocol==1 | DBt$Protocol==2) & DBt$Lieu_dit=="USN00-Réserve ornithologique (générique)") | ((DBt$Protocol==1 | DBt$Protocol==2)  & DBt$Lieu_dit=="USN01-Artigues-Réserve ornithologique")))

DBt$Date=as.Date(as.character(DBt$Date))
minAnnee = as.numeric(format(min(DBt$Date), format = "%Y"))
maxAnnee = as.numeric(format(max(DBt$Date), format = "%Y"))

DBt<-DBt[!colnames(DBt) %in% c("X","Ref")] #edit du 14/04/2017 : there are 3 duplicates
DBt = unique.matrix(DBt) 
limicoles = c("Recurvirostra avosetta","Limosa limosa","Limosa lapponica","Calidris temminckii","Calidris canutus",
              "Calidris alba","Calidris alpina","Calidris minuta","Calidris maritima" ,"Gallinago gallinago",
              "Tringa nebularia","Tringa erythropus","Tringa ochropus","Tringa totanus",
              "Tringa glareola","Actitis hypoleucos","Philomachus pugnax","Numenius arquata","Numenius phaeopus",
              "Himantopus himantopus","Charadrius hiaticula","Charadrius alexandrinus","Haematopus ostralegus",
              "Burhinus oedicnemus","Charadrius dubius","Phalaropus lobatus","Pluvialis squatarola",
              "Pluvialis apricaria","Arenaria interpres","Vanellus vanellus","Calidris ferruginea")
#Tringa flavipes never appers with Protocol>0
listColor = c('darkorchid','#F0C300','pink','lightgreen','darkgreen','lightblue','darkred','#6600FF','orange','#BABABA','#40826D','#3A9D23','#C71585','#F88E55','#00ffff','yellow','darkgrey','royalblue1','chocolate','antiquewhite1','black','deeppink','lightcoral','mediumorchid1','mediumaquamarine','olivedrab1','orangered','springgreen1','thistle1')
#listColor = c('darkgreen','red','blue','#F0C300','pink','lightgreen','violet','lightblue','darkred','#6600FF','orange','#BABABA','#40826D','#3A9D23','#C71585','#F88E55','#00ffff','yellow','darkgrey','royalblue1','chocolate','antiquewhite1','black','deeppink','lightcoral','mediumorchid1','mediumaquamarine','olivedrab1','orangered','springgreen1','thistle1')
minimal_percent = 1.00001

filename_pdf="average_abundance_timeseries_v2.pdf"
#pdf(paste("OUT/",filename_pdf,sep=""),width=15,height=10)
pdf(paste("Submission_JAE/Revisions_R2/",filename_pdf,sep=""),width=15,height=10)
par(mfrow=c(3,2))
layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = TRUE),heights = c(1,1,1),widths = c(3,1))
#Prepare doughnuts
abondance_anas = subset(DBt,grepl("^Anas",DBt$Nom_latin,ignore.case = TRUE) & DBt$Annee>=1981)
abondance_calidris = subset(DBt,grepl("^Calidris",DBt$Nom_latin,ignore.case = TRUE)& DBt$Annee>=1981)
abondance_waders = subset(DBt, DBt$Nom_latin %in% limicoles & DBt$Annee>=1981)


## abondance totale
abondance_tot_anas = sum(abondance_anas$Nombre,na.rm=TRUE) # I crush seasonality
abondance_tot_calidris = sum(abondance_calidris$Nombre,na.rm=TRUE)
abondance_tot_waders = sum(abondance_waders$Nombre,na.rm=TRUE)

# matrice percent abondance without seasonality
mat_percent_abondance_tot_anas     = create_mat_percent(abondance_anas,abondance_tot_anas)
mat_percent_abondance_tot_anas=mat_percent_abondance_tot_anas[order(mat_percent_abondance_tot_anas[,'somme'],decreasing=T),]
mat_percent_abondance_tot_calidris = create_mat_percent(abondance_calidris,abondance_tot_calidris)
mat_percent_abondance_tot_calidris=mat_percent_abondance_tot_calidris[order(mat_percent_abondance_tot_calidris[,'somme'],decreasing=T),]
mat_percent_abondance_tot_waders   = create_mat_percent(abondance_waders,abondance_tot_waders)
mat_percent_abondance_tot_waders=mat_percent_abondance_tot_waders[order(mat_percent_abondance_tot_waders[,'somme'],decreasing=T),]


#db_warm_tot=read.csv("IN/warmseason_abundances_asdataframe_summed.csv",sep=";",header=T)
#db_cold_tot=read.csv("IN/coldseason_abundances_asdataframe_summed.csv",sep=";",header=T)
db_warm_tot=read.csv("IN/warmseason_abundances_asdataframe_summed_v2_wrarespecies.csv",sep=";",header=T)
db_cold_tot=read.csv("IN/coldseason_abundances_asdataframe_summed_v2_wrarespecies.csv",sep=";",header=T)


#ANAS
#db_warm=read.csv("IN/warmseason_anas_detailed.txt",sep=",",header=T)
#db_cold=read.csv("IN/coldseason_anas_detailed.txt",sep=",",header=T)
db_warm=read.csv("IN/warmseason_anas_detailed_v2_abundances_wrarespecies.txt",sep=",",header=T)
db_cold=read.csv("IN/coldseason_anas_detailed_v2_abundances_wrarespecies.txt",sep=",",header=T)

temp_species=rownames(mat_percent_abondance_tot_anas)

dd=unique(db_cold_tot$Date)
dd=dd[1:(length(dd)-1)]
plou=c()
temp=c()
for(d in dd){
	temp=c(temp,d,d+0.5)
	plou=c(plou,db_cold_tot$Abundance[db_cold_tot$Species=="Anas"&db_cold_tot$Date==d],db_warm_tot$Abundance[db_warm_tot$Species=="Anas"&db_warm_tot$Date==d])
}
par(mar=c(3.9,4.75,0.1,0.1))
plot(temp,log10(plou+1),col="black",xlab="",ylab="log10(average abundance)",cex.lab=1.5,cex.axis=1.5,xaxt="n",pch=c(21,24),type="o",lwd =2,ylim=c(-0.1,max(log10(plou),na.rm=T)),bg=c("blue","red"),las=1,cex=1.25)
#points(db_warm_tot$Date[idw]+0.5,db_warm_tot$Abundance[idw],col="black",pch=17,ylim=c(0,max(db_warm_tot$Abundance[idw],na.rm=T)),type="p",lwd =2)
for (s in 1:length(temp_species)){
	temp=c()
	plou=c()
	for(d in dd){
        	temp=c(temp,d,d+0.5)
        	plou=c(plou,db_cold$abundance_cold[db_cold$sp_all == temp_species[s]&db_cold$dates==d],db_warm$abundance_warm[db_warm$sp_all == temp_species[s]&db_warm$dates==d])
	}
	points(temp,log10(plou+1),col=listColor[s],bg=listColor[s],pch=c(21,24),lty=2,t="o",cex=1.25)
#  points(db_cold$dates[as.character(db_cold$sp_all) == temp_species[s]],log10(db_cold$abundance_cold[db_cold$sp_all==temp_species[s]]),col=listColor[s],bg=listColor[s],pch=21,t="o",lty=2)
#  points(db_warm$dates[as.character(db_warm$sp_all) == temp_species[s]]+0.5,log10(db_warm$abundance_warm[db_warm$sp_all==temp_species[s]]),col=listColor[s],bg=listColor[s],pch=23)
}
#text(c(0,1), "a)",cex = 2)
mtext("\n (a)",side=3,cex = 1.5, adj = 0,padj=1)

par(mar=c(0.1,0.1,0.1,0.1))
id_1=which(rownames(mat_percent_abondance_tot_anas)=="Anas clypeata")
rownames(mat_percent_abondance_tot_anas)[id_1]="Spatula clypeata"
id_1=which(rownames(mat_percent_abondance_tot_anas)=="Anas penelope")
rownames(mat_percent_abondance_tot_anas)[id_1]="Mareca penelope"
id_1=which(rownames(mat_percent_abondance_tot_anas)=="Anas querquedula")
rownames(mat_percent_abondance_tot_anas)[id_1]="Spatula querquedula"
id_1=which(rownames(mat_percent_abondance_tot_anas)=="Anas strepera")
rownames(mat_percent_abondance_tot_anas)[id_1]="Mareca strepera"
doughnut( c(mat_percent_abondance_tot_anas[,'percent'][round(mat_percent_abondance_tot_anas[,'percent'])>=minimal_percent]) , inner.radius=0.15, col=listColor[which(round(mat_percent_abondance_tot_anas[,'percent'])>=minimal_percent)],outer.radius =0.45,label.cex=1.3,coord.x=c(-0.5,1), init.angle=50)


#CALIDRIS
#db_warm=read.csv("IN/warmseason_calidris_detailed.txt",sep=",",header=T)
#db_cold=read.csv("IN/coldseason_calidris_detailed.txt",sep=",",header=T)
db_warm=read.csv("IN/warmseason_calidris_detailed_v2_abundances_wrarespecies.txt",sep=",",header=T)
db_cold=read.csv("IN/coldseason_calidris_detailed_v2_abundances_wrarespecies.txt",sep=",",header=T)
temp_species=rownames(mat_percent_abondance_tot_calidris)

plou=c()
temp=c()
for(d in dd){
        temp=c(temp,d,d+0.5)
        plou=c(plou,db_cold_tot$Abundance[db_cold_tot$Species=="Calidris"&db_cold_tot$Date==d],db_warm_tot$Abundance[db_warm_tot$Species=="Calidris"&db_warm_tot$Date==d])
}
par(mar=c(3.9,4.8,0.1,0.1))
plot(temp,log10(plou+1),col="black",xlab="",ylab="log10(average abundance)",cex.lab=1.5,cex.axis=1.5,xaxt="n",pch=c(21,24),type="o",lwd =2,ylim=c(-0.1,max(log10(plou),na.rm=T)),bg=c("blue","red"),las=1,cex=1.25)
for (s in 1:length(temp_species)){
        temp=c()
        plou=c()
        for(d in dd){
                temp=c(temp,d,d+0.5)
                plou=c(plou,db_cold$abundance_cold[db_cold$sp_all == temp_species[s]&db_cold$dates==d],db_warm$abundance_warm[db_warm$sp_all == temp_species[s]&db_warm$dates==d])
        }
        points(temp,log10(plou+1),col=listColor[s],bg=listColor[s],pch=c(21,24),lty=2,t="o",cex=1.25)
}
#text(c(0,1), "a)",cex = 2)
mtext("\n (b)",side=3,cex = 1.5, adj = 0,padj=1)


par(mar=c(0.1,0.1,0.1,0.1))

doughnut( c(mat_percent_abondance_tot_calidris[,'percent'][round(mat_percent_abondance_tot_calidris[,'percent'])>=minimal_percent]) , inner.radius=0.145, col=listColor[which(round(mat_percent_abondance_tot_calidris[,'percent'])>=minimal_percent)],outer.radius =0.45,label.cex=1.28,coord.x=c(-0.7,1) )



#WADERS
#db_warm=read.csv("IN/warmseason_waders_detailed.txt",sep=",",header=T)
#db_cold=read.csv("IN/coldseason_waders_detailed.txt",sep=",",header=T)
db_warm=read.csv("IN/warmseason_waders_detailed_v2_abundances_wrarespecies.txt",sep=",",header=T)
db_cold=read.csv("IN/coldseason_waders_detailed_v2_abundances_wrarespecies.txt",sep=",",header=T)
temp_species=rownames(mat_percent_abondance_tot_waders)


plou=c()
temp=c()
for(d in dd){
        temp=c(temp,d,d+0.5)
        plou=c(plou,db_cold_tot$Abundance[db_cold_tot$Species=="Waders"&db_cold_tot$Date==d],db_warm_tot$Abundance[db_warm_tot$Species=="Waders"&db_warm_tot$Date==d])
}
par(mar=c(3.9,4.75,0.1,0.1))
plot(temp,log10(plou+1),col="black",xlab="Year",ylab="log10(average abundance)",cex.lab=1.5,cex.axis=1.5,pch=c(21,24),type="o",lwd =2,ylim=c(-0.1,max(log10(plou),na.rm=T)),bg=c("blue","red"),las=1,cex=1.25)
for (s in 1:length(temp_species)){
        temp=c()
        plou=c()
        for(d in dd){
                temp=c(temp,d,d+0.5)
                plou=c(plou,db_cold$abundance_cold[db_cold$sp_all == temp_species[s]&db_cold$dates==d],db_warm$abundance_warm[db_warm$sp_all == temp_species[s]&db_warm$dates==d])
        }
        points(temp,log10(plou+1),col=listColor[s],bg=listColor[s],pch=c(21,24),lty=2,t="o",cex=1.25)
}
#text(c(0,1), "a)",cex = 2)
mtext("\n (c)",side=3,cex = 1.5, adj = 0,padj=1)

par(mar=c(0.1,0.1,0.1,0.1))
rownames(mat_percent_abondance_tot_waders)[1]="C. alpina"
doughnut( c(mat_percent_abondance_tot_waders[,'percent'][round(mat_percent_abondance_tot_waders[,'percent'])>=minimal_percent]) , inner.radius=0.15, col=listColor[which(round(mat_percent_abondance_tot_waders[,'percent'])>=minimal_percent)],outer.radius =0.45,label.cex=1.28,coord.x=c(-0.5,1),lty=1)
dev.off()

system(paste("cp OUT/",filename_pdf," Submission_JAE/Revisions/",filename_pdf,sep=""))
