#CP 2019/09/04 Just checking the autocorrelation in the yearly data

rm(list=ls())
graphics.off()

sp_to_ignore=c("Anas discors","Anas americana","Calidris melanotos","Calidris pusilla","Calidris ruficollis", "Calidris fuscicollis", "Calidris himantopus", "Burhinus oedicnemus","Phalaropus lobatus","Charadrius alexandrinus","Haematopus ostralegus","Calidris maritima","Aythya nyroca","Bucephala clangula","Melanitta nigra","Mergus serrator","Clangula hyemalis","Alopochen aegyptiaca", "Aix galericulata","Cygnus atratus","Tadorna ferruginea","Branta leucopsis","Anser fabalis","Anser albifrons","Cygnus cygnus","Mergus merganser","Anser brachyrhynchus")

mat_acf=array(NA,dim=c(2,2,67),dimnames=list(c("Warm","Cold"),c("Lag1","Lag2"),rep("NA",67)))

####Waders
db_warm=read.csv("IN/warmseason_waders_detailed.txt",sep=",",header=T)
db_warm=db_warm[,c(2,3,4)]
names(db_warm)=c("dates","sp_data_frame","abundance")
db_warm$sp_data_frame=as.character(db_warm$sp_data_frame)
limicoles=unique(db_warm$sp_data_frame)
limicoles=limicoles[!(limicoles %in%sp_to_ignore)]
limicoles=sort(limicoles)
i=0
for(s in limicoles){
i=i+1
dimnames(mat_acf)[[3]][i]=s
 mat_acf[1,,i]=acf(db_warm$abundance[db_warm$sp_data_frame==s],plot=F)$acf[2:3]
}

#Ducks
db_warm=read.csv("IN/warmseason_duck_detailed.txt",sep=",",header=T)
db_warm=db_warm[,c(2,3,4)]
names(db_warm)=c("dates","sp_data_frame","abundance")
db_warm$sp_data_frame=as.character(db_warm$sp_data_frame)
ducks=unique(db_warm$sp_data_frame)
ducks=ducks[!(ducks %in%sp_to_ignore)]
ducks=sort(ducks)
for(s in ducks){
if(!(s %in% dimnames(mat_acf)[[3]])){
i=i+1
dimnames(mat_acf)[[3]][i]=s
 mat_acf[1,,i]=acf(db_warm$abundance[db_warm$sp_data_frame==s],plot=F)$acf[2:3]
}else{
 mat_acf[1,,s]=acf(db_warm$abundance[db_warm$sp_data_frame==s],plot=F)$acf[2:3]
}
}


#Freq
db_warm=read.csv("IN/warmseason_freq_detailed.txt",sep=",",header=T)
db_warm=db_warm[,c(2,3,4)]
names(db_warm)=c("dates","sp_data_frame","abundance")
db_warm$sp_data_frame=as.character(db_warm$sp_data_frame)
freq=unique(db_warm$sp_data_frame)
freq=freq[!(freq %in%sp_to_ignore)]
freq=sort(freq)
for(s in freq){
if(!(s %in% dimnames(mat_acf)[[3]])){
i=i+1
dimnames(mat_acf)[[3]][i]=s
 mat_acf[1,,i]=acf(db_warm$abundance[db_warm$sp_data_frame==s],plot=F)$acf[2:3]
}else{
 mat_acf[1,,s]=acf(db_warm$abundance[db_warm$sp_data_frame==s],plot=F)$acf[2:3]
}
}

#######Cold
####Waders
db_warm=read.csv("IN/coldseason_waders_detailed.txt",sep=",",header=T)
db_warm=db_warm[,c(2,3,4)]
names(db_warm)=c("dates","sp_data_frame","abundance")
db_warm$sp_data_frame=as.character(db_warm$sp_data_frame)
limicoles=unique(db_warm$sp_data_frame)
limicoles=limicoles[!(limicoles %in%sp_to_ignore)]
limicoles=sort(limicoles)
for(s in limicoles){
 mat_acf[2,,s]=acf(db_warm$abundance[db_warm$sp_data_frame==s],plot=F)$acf[2:3]
}

#Ducks
db_warm=read.csv("IN/coldseason_duck_detailed.txt",sep=",",header=T)
db_warm=db_warm[,c(2,3,4)]
names(db_warm)=c("dates","sp_data_frame","abundance")
db_warm$sp_data_frame=as.character(db_warm$sp_data_frame)
ducks=unique(db_warm$sp_data_frame)
ducks=ducks[!(ducks %in%sp_to_ignore)]
for(s in ducks){
 mat_acf[2,,s]=acf(db_warm$abundance[db_warm$sp_data_frame==s],plot=F)$acf[2:3]
}


#Freq
db_warm=read.csv("IN/coldseason_freq_detailed.txt",sep=",",header=T)
db_warm=db_warm[,c(2,3,4)]
names(db_warm)=c("dates","sp_data_frame","abundance")
db_warm$sp_data_frame=as.character(db_warm$sp_data_frame)
freq=unique(db_warm$sp_data_frame)
freq=freq[!(freq %in%sp_to_ignore)]
for(s in freq){
 mat_acf[2,,s]=acf(db_warm$abundance[db_warm$sp_data_frame==s],plot=F)$acf[2:3]
}

pdf("OUT/acf_birds.pdf")
par(mfrow=c(2,2))
for(l in 1:2){
	hist(mat_acf["Warm",l,],col="orange",xlim=c(-0.3,1.0),ylim=c(0,15),main="Warm")
	mtext(paste("lag",l),2,line=2)
	hist(mat_acf["Cold",l,],col="blue",xlim=c(-0.3,1.0),ylim=c(0,15),main="Cold")
}
dev.off()
