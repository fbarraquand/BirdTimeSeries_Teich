### CP 25/07/2019 Taking the code from CA to make it simpler to have the names of all species for Appendices

graphics.off()
rm(list=ls())


DBt<-read.csv(file="IN/DBWithMonthlyPhotoTeich_completed.csv",header=TRUE,sep=",",dec=".")
DBt = subset(DBt,(((DBt$Protocol==1 | DBt$Protocol==2) & DBt$Lieu_dit=="USN00-Réserve ornithologique (générique)") | ((DBt$Protocol==1 | DBt$Protocol==2)  & DBt$Lieu_dit=="USN01-Artigues-Réserve ornithologique")))

DBt$Date=as.Date(as.character(DBt$Date))
minAnnee = as.numeric(format(min(DBt$Date), format = "%Y"))
maxAnnee = as.numeric(format(max(DBt$Date), format = "%Y"))

sp_to_ignore=c("Anas discors","Anas americana","Calidris melanotos","Calidris pusilla","Calidris ruficollis", "Calidris fuscicollis", "Calidris himantopus", "Burhinus oedicnemus","Phalaropus lobatus","Charadrius alexandrinus","Haematopus ostralegus","Calidris maritima","Aythya nyroca","Bucephala clangula","Melanitta nigra","Mergus serrator","Clangula hyemalis","Alopochen aegyptiaca", "Aix galericulata","Cygnus atratus","Tadorna ferruginea","Branta leucopsis","Anser fabalis","Anser albifrons","Cygnus cygnus","Mergus merganser","Anser brachyrhynchus")


total_total_bird=sum(DBt$Nombre,na.rm=TRUE)
DBt=subset(DBt,!DBt$Nom_latin %in% sp_to_ignore)


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

total_bird = sum(DBt$Nombre,na.rm=TRUE)
abondance_ducks = subset(DBt, DBt$Nom_latin %in% ducks)
abondance_waders = subset(DBt, DBt$Nom_latin %in% limicoles)
abondance_frequents_birds = subset(DBt, DBt$Nom_latin %in% freq)

# Ducks, 17 species
unique_ducks =unique(as.character(abondance_ducks$Nom_latin))
unique_ducks= unique_ducks[order(unique_ducks)]
tab_abondance_ducks= matrix(data = NA,ncol = 1,nrow=length(unique_ducks))
rownames(tab_abondance_ducks)=c(unique_ducks)
for (i in 1:length(unique_ducks)){
  tab_abondance_ducks[as.character(unique_ducks[i]),] = round((sum(abondance_ducks$Nombre[abondance_ducks$Nom_latin==unique_ducks[i]],na.rm=TRUE)/total_bird)*100,3)
}

#WADERS, 26 species
unique_waders =unique(as.character(abondance_waders$Nom_latin))
unique_waders = unique_waders[order(unique_waders)]
tab_abondance_waders= matrix(data = NA,ncol = 1,nrow=length(unique_waders))
rownames(tab_abondance_waders)=c(unique_waders)
for (i in 1:length(unique_waders)){
  tab_abondance_waders[as.character(unique_waders[i]),] = round((sum(abondance_waders$Nombre[abondance_waders$Nom_latin==unique_waders[i]],na.rm=TRUE)/total_bird)*100,3)
}

# frequents birds
unique_frequents_birds =unique(as.character(abondance_frequents_birds$Nom_latin))
unique_frequents_birds=unique_frequents_birds[!(unique_frequents_birds %in% unique_ducks) & !(unique_frequents_birds %in% unique_waders)]
unique_frequents_birds = unique_frequents_birds[order(unique_frequents_birds)]
tab_abondance_frequents_birds= matrix(data = NA,ncol = 1,nrow=length(unique_frequents_birds))
rownames(tab_abondance_frequents_birds)=c(unique_frequents_birds)
for (i in 1:length(unique_frequents_birds)){
  tab_abondance_frequents_birds[as.character(unique_frequents_birds[i]),] = round((sum(abondance_frequents_birds$Nombre[abondance_frequents_birds$Nom_latin==unique_frequents_birds[i]],na.rm=TRUE)/total_bird)*100,3)
}

tmp=total_total_bird-sum(abondance_ducks$Nombre)-sum(abondance_waders$Nombre)-sum(abondance_frequents_birds$Nombre[abondance_frequents_birds$Nom_latin %in% unique_frequents_birds])
print(tmp/total_total_bird)

SpeciesL = as.character(unique(DBt$Nom_latin)) #-> 280
SpeciesL_order =SpeciesL[order(SpeciesL)]
mat_n_obs_oiseaux_t=matrix(data=NA,nrow=length(SpeciesL),ncol=1)
rownames(mat_n_obs_oiseaux_t)=c(SpeciesL_order)

for (i in 1:length(SpeciesL_order)){
  mat_n_obs_oiseaux_t[as.character(SpeciesL_order[i]),]=sum(as.character(DBt$Nom_latin)==SpeciesL_order[i],na.rm = TRUE)
}

# anas
subset(mat_n_obs_oiseaux_t,rownames(mat_n_obs_oiseaux_t) %in% ducks)
subset(mat_n_obs_oiseaux_t,rownames(mat_n_obs_oiseaux_t) %in% limicoles)
subset(mat_n_obs_oiseaux_t,rownames(mat_n_obs_oiseaux_t) %in% unique_frequents_birds)

######## Small modif by CP (11/07/19) to have relative number of occurrences and print abundance and occurrence in a tex formatted table
#For relative frequency
occ_ducks=subset(mat_n_obs_oiseaux_t/length(unique(DBt$Date)),rownames(mat_n_obs_oiseaux_t) %in% ducks)
occ_waders=subset(mat_n_obs_oiseaux_t/length(unique(DBt$Date)),rownames(mat_n_obs_oiseaux_t) %in% limicoles)
occ_freq=subset(mat_n_obs_oiseaux_t/length(unique(DBt$Date)),rownames(mat_n_obs_oiseaux_t) %in% unique_frequents_birds)

library(xtable)
tab=cbind(occ_freq*100,tab_abondance_frequents_birds)
print.xtable(xtable(tab,digits=c(1,1,3)),"OUT/abundance_occurrence_freq.tex" ,type="latex")

tab=cbind(occ_ducks*100,tab_abondance_ducks)
print.xtable(xtable(tab,digits=c(1,1,3)),"OUT/abundance_occurrence_ducks.tex" ,type="latex")

tab=cbind(occ_waders*100,tab_abondance_waders)
print.xtable(xtable(tab,digits=c(1,1,3)),"OUT/abundance_occurrence_waders.tex" ,type="latex")

