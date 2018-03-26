rm(list=ls())
graphics.off()

DBt<-read.csv(file="/home/cpicoche/Documents/Data_to_treat/TRANSFERT_LIMICOLES/IN/DBWithMonthlyPhotoTeich.csv",header=TRUE,sep=",",dec=".")
DBt = subset(DBt,(DBt$Protocol==1 | DBt$Protocol==2) & DBt$Lieu_dit=="USN00-Réserve ornithologique (générique)")
DBt<-DBt[!colnames(DBt) %in% c("X","Ref")] #edit du 14/04/2017 : vu qu'il y a 3 doublons
DBt = unique.matrix(DBt) # élimination de cette manière, mais il faudrait aller regarder dans la photo mensuelle 
# les doublons sont présents dans le fichier d'origine, 
# pour les virer avec unique.matrix je suis obligée de sortir les deux premières colonnes

#-----------------------------------------------------------------------------------------
#Combien d'espèces? 
SpeciesL = as.character(unique(DBt$Nom_latin))
vec_n_obs_oiseaux_t=rep(0,length(SpeciesL))
for (i in 1:length(SpeciesL)){
  vec_n_obs_oiseaux_t[i]=sum(as.character(DBt$Nom_latin)==SpeciesL[i],na.rm = TRUE)
}
oiseaux_Frequents_t=SpeciesL[vec_n_obs_oiseaux_t>75]


data_migration=read.table('/home/cpicoche/Documents/Data_to_treat/TRANSFERT_LIMICOLES/IN/oiseaux.csv',header=TRUE,sep=";")
Species_migration=as.character(unique(data_migration$Nom_Latin))
data_migration$TYPE=gsub("Nicheur","Breeding",data_migration$TYPE)
data_migration$TYPE=gsub("Hivernant","Wintering",data_migration$TYPE)
data_migration$TYPE=gsub("Migrateur","Migratory",data_migration$TYPE)
data_migration$TYPE=gsub("occasionnel","occasional",data_migration$TYPE)
data_migration$TYPE=gsub("sédentaire","sedentary",data_migration$TYPE)


cross_species=intersect(SpeciesL,Species_migration)
#cross_species=intersect(oiseaux_Frequents_t,Species_migration)

species_not_recognized=length(SpeciesL)-length(cross_species)

data_migration=subset(data_migration,Nom_Latin %in% cross_species)

ty=unique(data_migration$TYPE)

plop=table(data_migration$TYPE)
plop=sort(plop[plop>0])

pdf("OUT/known_migratory_birds.pdf")
par(mar=c(2,18,1,1))
barplot(plop,horiz=TRUE,las=1,cex.names=0.7)
dev.off()


