rm(list=ls())
graphics.off()
library('codyn')
source("SCRIPTS/test_synchrony_Gross.r")

mm=c('','_using_max_values')
ss=c('freq','Anas','Calidris','limicoles')

for (s in ss){
	for (m in mm){
filename=paste('OUT/data_',s,m,'.RData',sep="")
print(filename)
load(filename)

swinter_all=community_sync_Gross(winter_all,nrands=100)
sspring_all=community_sync_Gross(spring_all,nrands=100)
ssummer_all=community_sync_Gross(summer_all,nrands=100)
sautumn_all=community_sync_Gross(autumn_all,nrands=100)
swintering_all=community_sync_Gross(wintering_all,nrands=100)
sbreeding_all=community_sync_Gross(breeding_all,nrands=100)

swinter_pre_2006=community_sync_Gross(winter_pre_2006,nrands=100)
sspring_pre_2006=community_sync_Gross(spring_pre_2006,nrands=100)
ssummer_pre_2006=community_sync_Gross(summer_pre_2006,nrands=100)
sautumn_pre_2006=community_sync_Gross(autumn_pre_2006,nrands=100)
swintering_pre_2006=community_sync_Gross(wintering_pre_2006,nrands=100)
sbreeding_pre_2006=community_sync_Gross(breeding_pre_2006,nrands=100)

swinter_post_2006=community_sync_Gross(winter_post_2006,nrands=100)
sspring_post_2006=community_sync_Gross(spring_post_2006,nrands=100)
ssummer_post_2006=community_sync_Gross(summer_post_2006,nrands=100)
sautumn_post_2006=community_sync_Gross(autumn_post_2006,nrands=100)
swintering_post_2006=community_sync_Gross(wintering_post_2006,nrands=100)
sbreeding_post_2006=community_sync_Gross(breeding_post_2006,nrands=100)

data=list()
data$winter=c(swinter_all,swinter_pre_2006,swinter_post_2006)
data$spring=c(sspring_all,sspring_pre_2006,sspring_post_2006)
data$summer=c(ssummer_all,ssummer_pre_2006,ssummer_post_2006)
data$autumn=c(sautumn_all,sautumn_pre_2006,sautumn_post_2006)
data_freq_season4_Gross=data

data=list()
data$wintering=c(swintering_all,swintering_pre_2006,swintering_post_2006)
data$breeding=c(sbreeding_all,sbreeding_pre_2006,sbreeding_post_2006)
data_freq_season2_Gross=data

save(data_freq_season2_Gross,data_freq_season4_Gross,file=paste("indices_",s,m,"_with_p_values.RData",sep=""))
}
}

