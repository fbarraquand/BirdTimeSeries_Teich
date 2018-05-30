# Script by CPicoche 2018
#The aim is to compare synchrony indices at different taxonomic levels
rm(list=ls())
graphics.off()
library("R.utils")
source("SCRIPTS/test_synchrony_Gross.r")

thresh=0.05

# Grand Cormoran Phalacrocorax carbo
# Héron cendré Ardea cinerea
# Aigrette garzette Egretta garzetta

freq=loadToEnv("OUT/data_freq.RData")

#Cut years before 1981
ddmin=as.numeric(as.Date("1981-01-01"))
freq$wintering_all=subset(freq$wintering_all,dates>ddmin)
freq$wintering_pre_2006=subset(freq$wintering_pre_2006,dates>ddmin)
freq$breeding_all=subset(freq$breeding_all,dates>ddmin)
freq$breeding_pre_2006=subset(freq$breeding_pre_2006,dates>ddmin)

plou=list()
#Compute synchrony
val=c("wintering_all","wintering_pre_2006","wintering_post_2006","breeding_all","breeding_pre_2006","breeding_post_2006")
color=rep(c("Black","Lightblue","Darkblue"),2)
for (nid in 1:length(val)){
        n=val[nid]
        cormoran_all=subset(freq[[n]],sp_data_frame=="Phalacrocorax carbo")
        names(cormoran_all)=c("sp_data_frame","dates","abundance")
        heron_aigrette_tmp=subset(freq[[n]],sp_data_frame %in% c("Ardea cinerea","Egretta garzetta"))
        heron_aigrette_all=aggregate(heron_aigrette_tmp[,3],list(heron_aigrette_tmp$dates),sum)
        heron_aigrette_all=cbind(rep("HeronAigrette",length(heron_aigrette_all[[1]])),heron_aigrette_all)
        names(heron_aigrette_all)=c("sp_data_frame","dates","abundance")
                new_data_frame=merge(cormoran_all,heron_aigrette_all,all=T)
                plou[[nid]]=community_sync_Gross(new_data_frame,nrands=100)
}

