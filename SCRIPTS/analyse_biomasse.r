#########
##CP 2019/09/03: This script checks that we know the biomass of every species that is used in the analyses
########

rm(list=ls())
graphics.off()
library(stringr)

sp_to_ignore=c("Anas discors","Anas americana","Calidris melanotos","Calidris pusilla","Calidris ruficollis", "Calidris fuscicollis", "Calidris himantopus", "Burhinus oedicnemus","Phalaropus lobatus","Charadrius alexandrinus","Haematopus ostralegus","Calidris maritima","Aythya nyroca","Bucephala clangula","Melanitta nigra","Mergus serrator","Clangula hyemalis","Alopochen aegyptiaca", "Aix galericulata","Cygnus atratus","Tadorna ferruginea","Branta leucopsis","Anser fabalis","Anser albifrons","Cygnus cygnus","Mergus merganser","Anser brachyrhynchus")

####Waders
db_warm=read.csv("IN/warmseason_waders_detailed.txt",sep=",",header=T)
db_warm=db_warm[,c(2,3,4)]
names(db_warm)=c("dates","sp_data_frame","abundance")
db_warm$sp_data_frame=as.character(db_warm$sp_data_frame)
limicoles=unique(db_warm$sp_data_frame)
limicoles=limicoles[!(limicoles %in%sp_to_ignore)]

#Ducks
db_warm=read.csv("IN/warmseason_duck_detailed.txt",sep=",",header=T)
db_warm=db_warm[,c(2,3,4)]
names(db_warm)=c("dates","sp_data_frame","abundance")
db_warm$sp_data_frame=as.character(db_warm$sp_data_frame)
ducks=unique(db_warm$sp_data_frame)
ducks=ducks[!(ducks %in%sp_to_ignore)]

#Freq
db_warm=read.csv("IN/warmseason_freq_detailed.txt",sep=",",header=T)
db_warm=db_warm[,c(2,3,4)]
names(db_warm)=c("dates","sp_data_frame","abundance")
db_warm$sp_data_frame=as.character(db_warm$sp_data_frame)
freq=unique(db_warm$sp_data_frame)
freq=freq[!(freq %in%sp_to_ignore)]

sp_ok=unique(c(limicoles,ducks,freq))

###Open list of species for which we know the biomass
if(1==0){
tab_biom=read.table("IN/Information_Trait_Oiseaux_20190903.csv",header=T,sep=";")
sp_biom=as.character(tab_biom$Species)
sp_diff=setdiff(sp_ok,sp_biom)

vec=cbind(sp_diff,rep(NA,length(sp_diff)))
colnames(vec)=c("Species","Mass")
tab_biom_bis=rbind(tab_biom,vec)
write.table(tab_biom_bis,"IN/Information_Trait_Oiseaux_20190903_allimportantbirds.csv",sep=";",na = "",row.names=F,col.names=T)
} #This part to be used only once to be able to identify the species we were missing.

if(1==0){
tab_biom=read.table("IN/Information_Trait_Oiseaux_20190903_allimportantbirds.csv",header=T,sep=";")
for(i in 1:nrow(tab_biom)){
  tab_biom[i,3]=mean(as.numeric(str_split(tab_biom[i,2],"-")[[1]]))
}
write.table(tab_biom,"IN/Information_Trait_Oiseaux_20190903_allimportantbirds_meanmass.csv",sep=";",na = "",row.names=F,col.names=T)
} #This part to be used only once to compute average masses

