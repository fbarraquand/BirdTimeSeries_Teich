
rm(list=ls())
graphics.off()
source("SCRIPTS/analyse_oiseaux_saisonnier.R")

ssp=c('Anas','Calidris','freq','limicoles')
mmeanmax=c('mean','max')

for (s in ssp){
	for (m in mmeanmax){
		print(paste(s,m))
		prepare_data(s,m)
	}
}

