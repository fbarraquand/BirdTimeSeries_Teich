## CP 08/07/2020 Script to plot the two wavelet analyses of waders and waterfowl together
## Modified to use iaaft surrogates instead of phase shifts in the wavelet domain

rm(list=ls())
graphics.off()
source("SCRIPTS/iaaft.R")
library('mvcwt')
source("SCRIPTS/image_mvcwt_for_colormaps.r")
source("SCRIPTS/image_mvcwt_for_pvalues.r")
library("RColorBrewer")
library("lubridate")

set.seed(42)

anrands=1000 #Number of surrogates

biomass=T #If biomass==False, we only look at abundances. If biomass==TRUE, we weigh the abundances with the mass of each bird

if(biomass){
end_bio="biomasses"
}else{
end_bio="abundances"
}

#normalize_seq=c(T,F)
normalize_seq=F
doyouload=F #True if analyses are taken from previous files, False if we want to relaunch an analysis

for(normalize in normalize_seq){

if(normalize){
end_nor="scaled"
}else{
end_nor="NOTscaled"
}

##################################################################################
# -- importation des données du Teich
DBt<-read.csv(file="IN/DBWithMonthlyPhotoTeich_completed.csv",header=TRUE,sep=",",dec=".")
DBt = subset(DBt,(DBt$Protocol==1 | DBt$Protocol==2) & (DBt$Lieu_dit=="USN00-Réserve ornithologique (générique)" | DBt$Lieu_dit=="USN01-Artigues-Réserve ornithologique")  )#& !(DBt$Annee %in% c(2007,2008)))
DBt$Date=as.Date(as.character(DBt$Date))
minAnnee = as.numeric(format(min(DBt$Date), format = "%Y"))
maxAnnee = as.numeric(format(max(DBt$Date), format = "%Y"))

DBt<-DBt[!colnames(DBt) %in% c("X","Ref")] #edit du 14/04/2017 : vu qu'il y a 3 doublons
DBt = unique.matrix(DBt) # élimination de cette manière, mais il faudrait aller regarder dans la photo mensuelle 
# les doublons sont présents dans le fichier d'origine, 
# pour les virer avec unique.matrix je suis obligée de sortir les deux premières colonnes

#-----------------------------------------------------------------------------------------
#Combien d'espèces? 
limicoles = sort(c("Recurvirostra avosetta","Limosa limosa","Limosa lapponica","Calidris temminckii","Calidris canutus",
              "Calidris alba","Calidris alpina","Calidris minuta","Calidris maritima" ,"Gallinago gallinago",
              "Tringa nebularia","Tringa erythropus","Tringa ochropus","Tringa totanus",
              "Tringa glareola","Actitis hypoleucos","Philomachus pugnax","Numenius arquata","Numenius phaeopus",
              "Himantopus himantopus","Charadrius hiaticula","Charadrius alexandrinus","Haematopus ostralegus",
              "Burhinus oedicnemus","Charadrius dubius","Phalaropus lobatus","Pluvialis squatarola",
              "Pluvialis apricaria","Arenaria interpres","Vanellus vanellus","Calidris ferruginea"))

sp_to_ignore=c("Anas discors","Anas americana","Calidris melanotos","Calidris pusilla","Calidris ruficollis", "Calidris fuscicollis", "Calidris himantopus", "Burhinus oedicnemus","Phalaropus lobatus","Charadrius alexandrinus","Haematopus ostralegus","Calidris maritima","Aythya nyroca","Bucephala clangula","Melanitta nigra","Mergus serrator","Clangula hyemalis","Alopochen aegyptiaca", "Aix galericulata","Cygnus atratus","Tadorna ferruginea","Branta leucopsis","Anser fabalis","Anser albifrons","Cygnus cygnus","Mergus merganser","Anser brachyrhynchus")

#Waders are predefined, we just remove the species for which there were problems (observations only after or before 2006, which means that they would be 0 for one of the dataset)
limicoles=limicoles[!(limicoles %in% sp_to_ignore)]

tab_tmp=read.csv(file="IN/Initial_files/data_ROT20160324.csv",header=TRUE,sep="\t")
sp=tab_tmp$Nom_latin
fam=tab_tmp$Famille

#Waterfowl are defined as species from the Anatidae family + Fulica atra which is functionnally part of this group
ducks=c(as.character(unique(sp[fam=="Anatidae"])),"Fulica atra")
ducks=ducks[!(ducks %in% sp_to_ignore)]
ducks=intersect(ducks,DBt$Nom_latin)

year_min=1980
DBt=subset(DBt,year(Date)>year_min) #Before 1980, waders were not registered

dates=sort(unique(DBt$Date))

tab_limicoles=matrix(0,nrow=length(dates),ncol=length(limicoles))
colnames(tab_limicoles)=limicoles
tab_ducks=matrix(0,nrow=length(dates),ncol=length(ducks))
colnames(tab_ducks)=ducks

#We extract the two different tables of waders and waterfowl from the complete database
######If we want to analyse the biomasses, we need to weigh the abundances
tab_bm=read.table("IN/Information_Trait_Oiseaux_20190903_allimportantbirds_meanmass.csv",sep=";",dec=".",header=T)

for(id in 1:length(dates)){
        for (s in limicoles){
                id_d=which(DBt$Date==dates[id]&DBt$Nom_latin==s)
                if(length(id_d)>0){
                        tab_limicoles[id,s]=DBt$Nombre[id_d]
                }
                if(biomass){
                  tab_limicoles[id,s]=tab_limicoles[id,s]*tab_bm$Mean_mass[as.character(tab_bm$Species)==s]
		  if(length(tab_bm$Mean_mass[as.character(tab_bm$Species)==s])==0){stop("limicoles pas bm")}
                }
        }
        for (s in ducks){
                id_d=which(DBt$Date==dates[id]&DBt$Nom_latin==s)
                if(length(id_d)>0){
                        tab_ducks[id,s]=DBt$Nombre[id_d]
                }
                if(biomass){
                  tab_ducks[id,s]=tab_ducks[id,s]*tab_bm$Mean_mass[as.character(tab_bm$Species)==s]
		  if(length(tab_bm$Mean_mass[as.character(tab_bm$Species)==s])==0){stop("ducks pas bm")}
                }
        }

}

for(s in limicoles){
if(normalize){
        tab_limicoles[,s]=scale(tab_limicoles[,s])
}
}

for(s in ducks){
if(normalize){
        tab_ducks[,s]=scale(tab_ducks[,s])
}
}

x=(dates-dates[1])/365.25

#This function computes the Morlet wavelet transform for each bird species separately
print("Limicoles")
print(Sys.time())
#We need to regularize the dates as species were not always registered at the same day of the month after 2006
seq_x=seq(1,365.25*35+3*30.5,length.out=423)/365.25

mm=mvcwt(x,tab_limicoles,min.scale=mean(diff(seq_x))*3,max.scale=10.0,nscales=100,loc=seq_x)


year_min=1981

#This function computes the wavelet ratio of the whole community (see Keitt's paper in 2008)
if(!doyouload){
ref_wmr_wader=wmr(mm) #First, we compute the observed wavelet modulus ratio, from which we extract a table of the reference vaues
ref_val=ref_wmr_wader$z[,,1]
ref_wmr_wader$x=ref_wmr_wader$x+year_min #Change the dates to be "human-readable"

#Now, we store the values that could be observed under H0 with the iaaft
tab_values_iaaft=array(NA,dim=c(length(mm$x),length(mm$y),anrands+1))
tab_values_iaaft[,,anrands+1]=ref_val #We store the reference value as it also used to compute the final p-value
prog.bar = txtProgressBar(min = 0, max = anrands,style = 3)
for(i in 1:anrands){
        setTxtProgressBar(prog.bar, i)
        tab_tmp=tab_limicoles
	for(s in limicoles){
        tab_tmp[,s]=iaaft_surrogate(tab_limicoles[,s]) #Compute a surrogate with the iaaft method, for each species in the wader group
	}
	mmtmp=mvcwt(x,tab_tmp,min.scale=mean(diff(seq_x))*3,max.scale=10.0,nscales=100,loc=seq_x) #Compute the wavelet decompositon for each time series
        wmr_tmp=wmr(mmtmp) #Compute the wavelet ratio based on the surrogate time series
        tab_values_iaaft[,,i]=wmr_tmp$z[,,1] #Store the corresponding wavelet modulus ratio
}

#Now compute Pr(X<=x) (we will compute the corresponding p-values in the image. We stick to what is stored in z.boot in normal wmr.boot function
tab_pval=array(NA,dim=c(length(mm$x),length(mm$y),1))
for(i in 1:length(mm$x)){
        for(j in 1:length(mm$y)){
                #tab_pval[i,j,1]= 2*min(sum(tab_values_iaaft[i,j,] >= ref_val[i,j]),sum(tab_values_iaaft[i,j,] < ref_val[i,j]))/(anrands+1) #With this line, we could extract directly the p-value
                tab_pval[i,j,1]= sum(tab_values_iaaft[i,j,] <= ref_val[i,j])/(anrands+1)
                if(tab_pval[i,j,1]>1){stop()}

        }
}
ref_wmr_wader$z.boot=tab_pval

#locations and scales do not have the same length, we need to make it homogeneous before storing everything
if(length(ref_wmr_wader$x)>length(ref_wmr_wader$y)){
        yy=c(ref_wmr_wader$y,rep(NA,length(ref_wmr_wader$x)-length(ref_wmr_wader$y)))
        xx=ref_wmr_wader$x
}else{
        xx=c(ref_wmr_wader$x,rep(NA,length(ref_wmr_wader$y)-length(ref_wmr_wader$x)))
        yy=ref_wmr_wader$y
}
tab_xy=cbind(xx,yy)
colnames(tab_xy)=c("x","y")
write.table(tab_xy,paste("../Teich_resultsLFS/IAAFT_analyses_Gross100-1000surrogates/tab_xy_mr_waders",end_bio,"_",end_nor,"_with",anrands,"_IAAFT.csv",sep=""),sep=";",dec=".",col.names=T,row.names=F)

tab_z=ref_wmr_wader$z
write.table(as.matrix(tab_z[,,1]),paste("../Teich_resultsLFS/IAAFT_analyses_Gross100-1000surrogates/tab_z_mr_waders",end_bio,"_",end_nor,"_with",anrands,"_IAAFT.csv",sep=""),sep=";",dec=".",col.names=F,row.names=F)

tab_z.boot=ref_wmr_wader$z.boot
write.table(as.matrix(tab_z.boot[,,1]),paste("../Teich_resultsLFS/IAAFT_analyses_Gross100-1000surrogates/tab_zboot_mr_waders",end_bio,"_",end_nor,"_with",anrands,"_IAAFT.csv",sep=""),sep=";",dec=".",col.names=F,row.names=F)

}else{

#Load all data
ref_wmr_wader = wmr(mm)

tmp_xy=read.csv(paste("../Teich_resultsLFS/IAAFT_analyses_Gross100-1000surrogates/tab_xy_mr_waders",end_bio,"_",end_nor,"_with",anrands,"_IAAFT.csv",sep=""),header=T,sep=";",dec=".")
ref_wmr_wader$x=tmp_xy[!is.na(tmp_xy[,"x"]),"x"]
ref_wmr_wader$y=tmp_xy[!is.na(tmp_xy[,"y"]),"y"]

tmp_z=as.matrix(read.csv(paste("../Teich_resultsLFS/IAAFT_analyses_Gross100-1000surrogates/tab_z_mr_waders",end_bio,"_",end_nor,"_with",anrands,"_IAAFT.csv",sep=""),header=F,sep=";",dec="."))
tmp_array_z=array(0,dim=c(dim(tmp_z),1))
tmp_array_z[,,1]=tmp_z
ref_wmr_wader$z=tmp_array_z

tmp_z.boot=as.matrix(read.csv(paste("../Teich_resultsLFS/IAAFT_analyses_Gross100-1000surrogates/tab_zboot_mr_waders",end_bio,"_",end_nor,"_with",anrands,"_IAAFT.csv",sep=""),header=F,sep=";",dec="."))
tmp_array_z.boot=array(0,dim=c(dim(tmp_z.boot),1))
tmp_array_z.boot[,,1]=tmp_z.boot
ref_wmr_wader$z.boot=tmp_array_z.boot
}


#Same thing as above, but for waterfowl
print("Canards")
print(Sys.time())
#We regularize the sequence of dates
seq_x=seq(1,365.25*35+3*30.5,length.out=423)/365.25

#We compute the wavelet decomposition for each species of the waterfowl group
mm=mvcwt(x,tab_ducks,min.scale=mean(diff(seq_x))*3,max.scale=10.0,nscales=100,loc=seq_x)

if(!doyouload){
#We compute the observed wavelet modulus ratio
ref_wmr_waterfowl=wmr(mm)
ref_val=ref_wmr_waterfowl$z[,,1]
ref_wmr_waterfowl$x=ref_wmr_waterfowl$x+year_min #Change the dates to be "human-readable"

#We compute the wavelet modulus ratio for the surrogates built under H0
tab_values_iaaft=array(NA,dim=c(length(mm$x),length(mm$y),anrands+1))
tab_values_iaaft[,,anrands+1]=ref_val
prog.bar = txtProgressBar(min = 0, max = anrands,style = 3)
for(i in 1:anrands){
        setTxtProgressBar(prog.bar, i)
        tab_tmp=tab_ducks
	for(s in ducks){
        tab_tmp[,s]=iaaft_surrogate(tab_ducks[,s])
	}
        #mmtmp=mvcwt(x,tab_tmp,min.scale=0.2,max.scale=10.0,nscales=100,loc=regularize(x,nsteps=length(x)/2))
	mmtmp=mvcwt(x,tab_tmp,min.scale=mean(diff(seq_x))*3,max.scale=10.0,nscales=100,loc=seq_x)
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
ref_wmr_waterfowl$z.boot=tab_pval

#We store all Pr(X<=x) as well as scales, locations and the observed values
if(length(ref_wmr_waterfowl$x)>length(ref_wmr_waterfowl$y)){
        yy=c(ref_wmr_waterfowl$y,rep(NA,length(ref_wmr_waterfowl$x)-length(ref_wmr_waterfowl$y)))
        xx=ref_wmr_waterfowl$x
}else{
        xx=c(ref_wmr_waterfowl$x,rep(NA,length(ref_wmr_waterfowl$y)-length(ref_wmr_waterfowl$x)))
        yy=ref_wmr_waterfowl$y
}

tab_xy=cbind(xx,yy)
colnames(tab_xy)=c("x","y")
write.table(tab_xy,paste("../Teich_resultsLFS/IAAFT_analyses_Gross100-1000surrogates/tab_xy_mr_waterfowl",end_bio,"_",end_nor,"_with",anrands,"_IAAFT.csv",sep=""),sep=";",dec=".",col.names=T,row.names=F)

tab_z=ref_wmr_waterfowl$z
write.table(as.matrix(tab_z[,,1]),paste("../Teich_resultsLFS/IAAFT_analyses_Gross100-1000surrogates/tab_z_mr_waterfowl",end_bio,"_",end_nor,"_with",anrands,"_IAAFT.csv",sep=""),sep=";",dec=".",col.names=F,row.names=F)

tab_z.boot=ref_wmr_waterfowl$z.boot
write.table(as.matrix(tab_z.boot[,,1]),paste("../Teich_resultsLFS/IAAFT_analyses_Gross100-1000surrogates/tab_zboot_mr_waterfowl",end_bio,"_",end_nor,"_with",anrands,"_IAAFT.csv",sep=""),sep=";",dec=".",col.names=F,row.names=F)


}else{

#We need to create a wmr object to store the values, so that it is recognized by the function image_mvcwt. We then store the values from previous analyses in it.
ref_wmr_waterfowl = wmr(mm)

tmp_xy=read.csv(paste("../Teich_resultsLFS/IAAFT_analyses_Gross100-1000surrogates/tab_xy_mr_waterfowl",end_bio,"_",end_nor,"_with",anrands,"_IAAFT.csv",sep=""),header=T,sep=";",dec=".")
ref_wmr_waterfowl$x=tmp_xy[!is.na(tmp_xy[,"x"]),"x"]
ref_wmr_waterfowl$y=tmp_xy[!is.na(tmp_xy[,"y"]),"y"]
ref_wmr_waterfowl$x=ref_wmr_waterfowl$x-year_min 
print("Here, we subtract 1981 because there was a mistake in the .csv file (we added 1981 twice)")

tmp_z=as.matrix(read.csv(paste("../Teich_resultsLFS/IAAFT_analyses_Gross100-1000surrogates/tab_z_mr_waterfowl",end_bio,"_",end_nor,"_with",anrands,"_IAAFT.csv",sep=""),header=F,sep=";",dec="."))
tmp_array_z=array(0,dim=c(dim(tmp_z),1))
tmp_array_z[,,1]=tmp_z
ref_wmr_waterfowl$z=tmp_array_z

tmp_z.boot=as.matrix(read.csv(paste("../Teich_resultsLFS/IAAFT_analyses_Gross100-1000surrogates/tab_zboot_mr_waterfowl",end_bio,"_",end_nor,"_with",anrands,"_IAAFT.csv",sep=""),header=F,sep=";",dec="."))
tmp_array_z.boot=array(0,dim=c(dim(tmp_z.boot),1))
tmp_array_z.boot[,,1]=tmp_z.boot
ref_wmr_waterfowl$z.boot=tmp_array_z.boot
}

#Finally, we plot both the waterfowl and wader wavelet modulus ratios
pdf(paste("Submission_JAE/Revisions_R2/wavelet_wader_waterfowl",end_bio,"_",end_nor,"_with",anrands,"nocorrection_smallergrid_IAAFT.pdf",sep=""),height=15,width=12)
layout(matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=T),widths=c(10,2))

#par(mar=c(4,5,3,3))
par(mar=c(6,6,3,3))


  image_mvcwt_for_colormaps(ref_wmr_wader,reset.par=F,cex.axis=4,z.fun="Mod",amain="Wader",adj="None")
#  image_mvcwt_for_pvalues(mr_object,reset.par=F,cex.axis=4,z.fun="Mod",amain="Wader",adj="None")
mtext("a)",side=2,line=-2,at=0.98,cex=1.5,outer=T,las=1)
mtext("Scale (years)", side=2, line=-2,at=0.775, outer = TRUE,cex=1.5)


#par(mar=c(4,5,3,3))
par(mar=c(6,6,3,3))


mtext("b)",side=2,line=-2,at=0.49,cex=1.5,outer=T,las=1)
mtext("Scale (years)", side=2, line=-2,at=0.28, outer = TRUE,cex=1.5)
mtext("Years",side=1, line=-2, at=0.425,outer = TRUE,cex=1.5)

image_mvcwt_for_colormaps(ref_wmr_waterfowl,reset.par=F,cex.axis=4,z.fun="Mod",amain="Waterfowl",adj="None")
dev.off()
}

