#Deprecated, C. Aluome and CP just wanted to compare their Databases
rm(list=ls())
graphics.off()
source("SCRIPTS/test_synchrony_Gross.r")

print("Calidris")
db_warm=read.csv("IN/cricri_matrix_warm.csv",sep=";",header=T)
db_cold=read.csv("IN/cricri_matrix_cold.csv",sep=";",header=T)

#db_warm=db_warm[,c(2,3,4)]
#db_cold=db_cold[,c(2,3,4)]
#Right format for synchrony scripts
names(db_warm)=c("dates","sp_data_frame","abundance")
names(db_cold)=c("dates","sp_data_frame","abundance")

list1=unique(db_warm$sp_data_frame)
sp_ignore=c("Calidris melanotos","Calidris pusilla","Calidris pusilla","Calidris ruficollis","Calidris fuscicollis","Calidris himantopus","Calidris maritima")
list2=list1[!(list1 %in% sp_ignore)]
list3=c("Calidris canutus", "Calidris alpina", "Calidris ferruginea","Calidris minuta")

#List1
db_warm_all=subset(db_warm)
db_warm_pre_2006=subset(db_warm,dates<13649)
db_warm_post_2006=subset(db_warm,dates>=13649)
db_cold_all=subset(db_cold)
db_cold_pre_2006=subset(db_cold,dates<13649)
db_cold_post_2006=subset(db_cold,dates>=13649)

#Compute synchrony values
synch_warm_all=community_sync_Gross(db_warm_all,nrands=100)
synch_warm_pre_2006=community_sync_Gross(db_warm_pre_2006,nrands=100)
synch_warm_post_2006=community_sync_Gross(db_warm_post_2006,nrands=100)

synch_cold_all=community_sync_Gross(db_cold_all,nrands=100)
synch_cold_pre_2006=community_sync_Gross(db_cold_pre_2006,nrands=100)
synch_cold_post_2006=community_sync_Gross(db_cold_post_2006,nrands=100)

list_calidris_1=list(synch_cold_all,synch_cold_pre_2006,synch_cold_post_2006,synch_warm_all,synch_warm_pre_2006,synch_warm_post_2006)

#List2
db_warm_all=subset(db_warm,sp_data_frame %in% list2 )
db_warm_pre_2006=subset(db_warm,sp_data_frame %in% list2  &dates<13532)
db_warm_post_2006=subset(db_warm,sp_data_frame %in% list2  &dates>=13532)
db_cold_all=subset(db_cold,sp_data_frame %in% list2 )
db_cold_pre_2006=subset(db_cold,sp_data_frame %in% list2 &dates<13532)
db_cold_post_2006=subset(db_cold, sp_data_frame %in% list2 &dates>=13532)

synch_warm_all=community_sync_Gross(db_warm_all,nrands=100)
synch_warm_pre_2006=community_sync_Gross(db_warm_pre_2006,nrands=100)
synch_warm_post_2006=community_sync_Gross(db_warm_post_2006,nrands=100)

synch_cold_all=community_sync_Gross(db_cold_all,nrands=100)
synch_cold_pre_2006=community_sync_Gross(db_cold_pre_2006,nrands=100)
synch_cold_post_2006=community_sync_Gross(db_cold_post_2006,nrands=100)

list_calidris_2=list(synch_cold_all,synch_cold_pre_2006,synch_cold_post_2006,synch_warm_all,synch_warm_pre_2006,synch_warm_post_2006)

#List3
db_warm_all=subset(db_warm,sp_data_frame %in% list3  &dates<2016)
db_warm_pre_2006=subset(db_warm,sp_data_frame %in% list3  &dates<=2006)
db_warm_post_2006=subset(db_warm,sp_data_frame %in% list3  &dates>2006&dates<2016)
db_cold_all=subset(db_cold,sp_data_frame %in% list3  &dates<2016)
db_cold_pre_2006=subset(db_cold,sp_data_frame %in% list3 &dates<=2006)
db_cold_post_2006=subset(db_cold, sp_data_frame %in% list3 &dates>2006&dates<2016)

synch_warm_all=community_sync_Gross(db_warm_all,nrands=100)
synch_warm_pre_2006=community_sync_Gross(db_warm_pre_2006,nrands=100)
synch_warm_post_2006=community_sync_Gross(db_warm_post_2006,nrands=100)

synch_cold_all=community_sync_Gross(db_cold_all,nrands=100)
synch_cold_pre_2006=community_sync_Gross(db_cold_pre_2006,nrands=100)
synch_cold_post_2006=community_sync_Gross(db_cold_post_2006,nrands=100)

list_calidris_3=list(synch_cold_all,synch_cold_pre_2006,synch_cold_post_2006,synch_warm_all,synch_warm_pre_2006,synch_warm_post_2006)

