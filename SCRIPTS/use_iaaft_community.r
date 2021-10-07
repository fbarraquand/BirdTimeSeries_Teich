# This script is based on Fig5.R which is now in the explo folder and is used to compare Gross and Keitt indices for the most frequent birds with IAAFT surrogates instead of Keitt's classical phase shifts for the wavelet

rm(list=ls())
graphics.off()
source("SCRIPTS/test_synchrony_Gross.r")
library('mvcwt')
source("SCRIPTS/image_mvcwt_for_colormaps.r") 
library("RColorBrewer")
library("lubridate")

set.seed(42)

thresh=0.1
type_correct="BH" #was Bonferroni before. Only used for the Gross index as we do no correction for the Keitt index
amethod_w="shift" #Method used for the surrogates for the Gross method
anrands=1000
biomass=F
if(biomass){
end_bio="biomasses"
}else{
end_bio="abundances"
}

doyouload=T #False if you want to launch the analyses again ; TRUE if you just want to do the plots

normalize_seq=c(TRUE,FALSE)

for(normalize in normalize_seq){
print(paste("normalize",normalize))

if(normalize){
end_nor="scaled"
}else{
end_nor="NOTscaled"
}

#Species to ignore as their frequency is too low or their observations might be concentrated during one phase of the analysis (before or after 2006), which can cause problems in the analyses
sp_to_ignore=c("Anas discors","Anas americana","Calidris melanotos","Calidris pusilla","Calidris ruficollis", "Calidris fuscicollis", "Calidris himantopus", "Burhinus oedicnemus","Phalaropus lobatus","Charadrius alexandrinus","Haematopus ostralegus","Calidris maritima","Aythya nyroca","Bucephala clangula","Melanitta nigra","Mergus serrator","Clangula hyemalis","Alopochen aegyptiaca", "Aix galericulata","Cygnus atratus","Tadorna ferruginea","Branta leucopsis","Anser fabalis","Anser albifrons","Cygnus cygnus","Mergus merganser","Anser brachyrhynchus")

print("Gross index")
print(Sys.time())
db_warm=read.csv(paste("IN/warmseason_freq_detailed_v2_",end_bio,"_wtoutrarespecies.txt",sep=""))
db_cold=read.csv(paste("IN/coldseason_freq_detailed_v2_",end_bio,"_wtoutrarespecies.txt",sep=""))
db_warm=db_warm[,c(2,3,4)]
db_cold=db_cold[,c(2,3,4)]

#Right format for synchrony scripts
names(db_warm)=c("dates","sp_data_frame","abundance")
names(db_cold)=c("dates","sp_data_frame","abundance")

#Subset to separate before and after 2006, remove species that we should ignore, and ignore dates before 2016 which can cause problems
db_warm_all=subset(db_warm,!(as.character(sp_data_frame) %in%sp_to_ignore) &dates<2016)
db_warm_pre_2006=subset(db_warm,!(as.character(sp_data_frame) %in%sp_to_ignore)  &dates<=2006)
db_warm_post_2006=subset(db_warm,!(as.character(sp_data_frame) %in%sp_to_ignore)  &dates>2006&dates<2016)
db_cold_all=subset(db_cold,!(as.character(sp_data_frame) %in%sp_to_ignore) &dates<2016)
db_cold_pre_2006=subset(db_cold,!(as.character(sp_data_frame) %in%sp_to_ignore) &dates<=2006)
db_cold_post_2006=subset(db_cold,!(as.character(sp_data_frame) %in%sp_to_ignore)&dates>2006&dates<2016)


if(normalize){
	list_db=list(db_cold_all,db_cold_pre_2006,db_cold_post_2006,db_warm_all,db_warm_pre_2006,db_warm_post_2006)
	for(d in 1:length(list_db)){
		species=unique(list_db[[d]]$sp_data_frame)
		for(s in species){
			list_db[[d]]$abundance[list_db[[d]]$sp_data_frame==s]=scale(list_db[[d]]$abundance[list_db[[d]]$sp_data_frame==s])
		}
	}
	db_cold_all=list_db[[1]]
	db_cold_pre_2006=list_db[[2]]
	db_cold_post_2006=list_db[[3]]
	db_warm_all=list_db[[4]]
	db_warm_pre_2006=list_db[[5]]
	db_warm_post_2006=list_db[[6]]
}



if(doyouload){
	mat_save=read.table(paste("../Teich_resultsLFS/IAAFT_analyses_Gross100-1000surrogates/tab_data_frame_Gross_community_",end_bio,"_",end_nor,"_with",anrands,".csv",sep=""),sep=";",dec=".",header=T)
	list_freq=list()
	for(v in 1:nrow(mat_save)){
		list_freq[[v]]=list(obs=as.numeric(mat_save[v,"obs"]),pval=as.numeric(mat_save[v,"pval"]),alternative=as.character(mat_save[v,"alternative"]),rands=as.numeric(c(mat_save[v,grep("rands",colnames(mat_save))])))
	}
	

}else{

#Compute synchrony values
synch_warm_all=community_sync_Gross(db_warm_all,nrands=anrands,method=amethod_w)
synch_warm_pre_2006=community_sync_Gross(db_warm_pre_2006,nrands=anrands,method=amethod_w)
synch_warm_post_2006=community_sync_Gross(db_warm_post_2006,nrands=anrands,method=amethod_w)

synch_cold_all=community_sync_Gross(db_cold_all,nrands=anrands,method=amethod_w)
synch_cold_pre_2006=community_sync_Gross(db_cold_pre_2006,nrands=anrands,method=amethod_w)
synch_cold_post_2006=community_sync_Gross(db_cold_post_2006,nrands=anrands,method=amethod_w)

list_freq=list(synch_cold_all,synch_cold_pre_2006,synch_cold_post_2006,synch_warm_all,synch_warm_pre_2006,synch_warm_post_2006)

#Save analyses in a separate file
mat_save=matrix(NA,nrow=length(list_freq),ncol=3+anrands+1) #3 for obs, pval, alternative ; nrands for all the values... of rands. We add 1 to nrands because we also use the observed value in the computation of the pvalues
colnames(mat_save)=c("obs",paste("rands",1:(anrands+1),sep=""),"pval","alternative")
for(v in 1:length(list_freq)){
	mat_save[v,"obs"]=list_freq[[v]]$obs
	mat_save[v,"pval"]=list_freq[[v]]$pval
	mat_save[v,"alternative"]=list_freq[[v]]$alternative
	for(r in 1:anrands){
		mat_save[v,paste("rands",r,sep="")]=list_freq[[v]]$rands[r]
	}
	mat_save[v,paste("rands",anrands+1,sep="")]=list_freq[[v]]$rands[anrands+1]
}
write.table(mat_save,paste("../Teich_resultsLFS/IAAFT_analyses_Gross100-1000surrogates/tab_data_frame_Gross_community_",end_bio,"_",end_nor,"_with",anrands,".csv",sep=""),sep=";",col.names=TRUE,row.names=F,dec=".")
write.table(mat_save,paste("../Teich_resultsLFS/IAAFT_analyses_Gross100-1000surrogates/tab_data_frame_Gross_community_",end_bio,"_",end_nor,"_with",anrands,".csv",sep=""),sep=";",col.names=TRUE,row.names=F,dec=".")

}


print(Sys.time())

#Adjust p-values for FDR
mat=rep(NA,length(list_freq))
for(v in 1:length(list_freq)){
        mat[v]=list_freq[[v]]$pval
}
mat_adj=p.adjust(mat,method=type_correct)
for(v in 1:length(list_freq)){
        list_freq[[v]]$pval=mat_adj[v]
}

essai_taxo=list_freq

#Plot Gross indices
upsi=0.02
color=rep(c("Black","Lightblue","Dodgerblue2"),2)
filename_pdf=paste("synchrony_indices_frequent_",end_bio,"_",end_nor,"_with",anrands,"rand_nocorrection_smallgrid_IAAFT_corrected_loc.pdf",sep="")
pdf(paste("Submission_JAE/Revisions_R2/",filename_pdf,sep=""),width=9,height=9)

layout(matrix(c(1,1,2,3),nrow=2,ncol=2,byrow=T),widths=c(10,2))
par(mar=c(3,5,2,3))

plot(0,0,t="n",ylim=c(-1.,1.0),xlim=c(0,7.5),xaxt="n",xlab="",ylab=expression(eta),cex.lab=1.5,cex.axis=1.5,las=1)
mtext("a)",side=2,line=-2,at=0.96,cex=1.5,outer=T,las=1)
axis(1,at=c(2,6),lab=c("Cold","Warm"),cex.axis=1.5)
for (v in 1:6){
        obs=essai_taxo[[v]]$obs
        print(obs)
        p_s=essai_taxo[[v]]$pval
        x=v+(v>3)
        points(x,obs,pch=21,col=color[v],bg=color[v],cex=2)
        boxplot(essai_taxo[[v]]$rands[1:(anrands)],at=x,add=T,boxwex=0.25,range=0,yaxt="n",xaxt="n")
        if (p_s<thresh){
                points(x,as.numeric(obs),pch='*',col="red",cex=2)
                }

        }
lines(c(0,7.5),c(0,0),lty=2,lwd=2)

legend("bottomleft",c("All","Pre-2006","Post-2006"),pch=NA,fill=c("black","Lightblue","Dodgerblue2"),pt.cex=2,bty="n",cex=1.5)


############################################ WAVELETS ########################################

# -- importation des données du Teich
DBt<-read.csv(file="IN/DBWithMonthlyPhotoTeich_completed.csv",header=TRUE,sep=",",dec=".")

DBt = subset(DBt,(((DBt$Protocol==1 | DBt$Protocol==2) & DBt$Lieu_dit=="USN00-Réserve ornithologique (générique)") | ((DBt$Protocol==1 | DBt$Protocol==2)  & DBt$Lieu_dit=="USN01-Artigues-Réserve ornithologique")))

DBt$Date=as.Date(as.character(DBt$Date))
minAnnee = as.numeric(format(min(DBt$Date), format = "%Y"))
maxAnnee = as.numeric(format(max(DBt$Date), format = "%Y"))

DBt<-DBt[!colnames(DBt) %in% c("X","Ref")] #edit du 14/04/2017 : there are 3 duplicates
DBt = unique.matrix(DBt) # delete with this line,  but I must check the monthly picture
# duplicates are present in the original file, it's not a code pb
# to delete with unique.matrix I must delete the first two columns.

#minimum of observation to keep a species = Frequent Species
SpeciesL = as.character(unique(DBt$Nom_latin)) #-> 280
vec_n_obs_oiseaux_t=rep(0,length(SpeciesL))
for (i in 1:length(SpeciesL)){
  vec_n_obs_oiseaux_t[i]=sum(as.character(DBt$Nom_latin)==SpeciesL[i],na.rm = TRUE)
}

freq_birds = SpeciesL[vec_n_obs_oiseaux_t>75] 

#Warning : I'm cutting the beginning of the series here, because there are no waders before 1980
year_min=1980
DBt=subset(DBt,year(Date)>year_min & Nom_latin%in%freq_birds)

dates=sort(unique(DBt$Date))

tab_freq=matrix(0,nrow=length(dates),ncol=length(freq_birds))
colnames(tab_freq)=freq_birds
#Build the table we want to analyse from the database with all species
for(id in 1:length(dates)){
        for (s in freq_birds){
                id_d=which(DBt$Date==dates[id]&DBt$Nom_latin==s)
                if(length(id_d)>0){
                        tab_freq[id,s]=DBt$Nombre[id_d]
                }
        }
}

#Normalize if necessary
for(s in freq_birds){
if(normalize){
        tab_freq[,s]=scale(tab_freq[,s])
}

}

x=(dates-dates[1])/365.25

year_min=1981
#This function computes the Morlet wavelet transform for each bird species separately
print("Keitt index")
print(Sys.time())

#####Regularize data because some of them have been collected at different dates during the month after 2006
seq_x=seq(1,365.25*35+3*30.5,length.out=423)/365.25
#Computes the wavelet decomposition for time series in the table
mm=mvcwt(x,tab_freq,min.scale=mean(diff(seq_x))*3,max.scale=10.0,nscales=100,loc=seq_x)
print(paste(Sys.time(),"after mvcwt"))

#This function computes the wavelet ratio of the whole community (see Keitt's paper in 2008)
if(!doyouload){
#Computes the observed wavelet modulus ratio from the wavelet decomposition
ref_wmr=wmr(mm)
ref_val=ref_wmr$z[,,1]
ref_wmr$x=ref_wmr$x+year_min #Change the dates to be "human-readable"

#Computes IAAFT surrogates under H0 and the corresponding wavelet modulus ratio
tab_values_iaaft=array(NA,dim=c(length(mm$x),length(mm$y),anrands+1))
tab_values_iaaft[,,anrands+1]=ref_val
prog.bar = txtProgressBar(min = 0, max = anrands,style = 3)
for(i in 1:anrands){
        setTxtProgressBar(prog.bar, i)
        tab_tmp=tab_freq
        for(s in freq_birds){
        tab_tmp[,s]=iaaft_surrogate(tab_freq[,s])
        }
	mmtmp=mvcwt(x,tab_tmp,min.scale=mean(diff(seq_x))*3,max.scale=10.0,nscales=100,loc=seq_x)
        wmr_tmp=wmr(mmtmp)
        tab_values_iaaft[,,i]=wmr_tmp$z[,,1]
}

#Now compute the one-tailed p-value Pr(X<=x_obs) where X is the test statistic, for all pixels in the image. The switch to two-tailed p-values is done when calling image_mvcwt_for_colormaps.r
tab_pval=array(NA,dim=c(length(mm$x),length(mm$y),1))
for(i in 1:length(mm$x)){
        for(j in 1:length(mm$y)){
                tab_pval[i,j,1]= sum(tab_values_iaaft[i,j,] <= ref_val[i,j])/(anrands+1)
                if(tab_pval[i,j,1]>1){stop()}

        }
}
ref_wmr$z.boot=tab_pval

#there is not the same number of locations and scales, so we need to complete the vectors before being able to bind them
if(length(ref_wmr$x)>length(ref_wmr$y)){
	yy=c(ref_wmr$y,rep(NA,length(ref_wmr$x)-length(ref_wmr$y)))
	xx=ref_wmr$x
}else{
	xx=c(ref_wmr$x,rep(NA,length(ref_wmr$y)-length(ref_wmr$x)))
	yy=ref_wmr$y
}
tab_xy=cbind(xx,yy)
colnames(tab_xy)=c("x","y")
write.table(tab_xy,paste("../Teich_resultsLFS/IAAFT_analyses_Gross100-1000surrogates/tab_xy_mr_freq",end_bio,"_",end_nor,"_with",anrands,"_IAAFT_corrected_loc.csv",sep=""),sep=";",dec=".",col.names=T,row.names=F)

tab_z=ref_wmr$z
write.table(as.matrix(tab_z[,,1]),paste("../Teich_resultsLFS/IAAFT_analyses_Gross100-1000surrogates/tab_z_mr_freq",end_bio,"_",end_nor,"_with",anrands,"_IAAFT_corrected_loc.csv",sep=""),sep=";",dec=".",col.names=F,row.names=F)

tab_z.boot=ref_wmr$z.boot
write.table(as.matrix(tab_z.boot[,,1]),paste("../Teich_resultsLFS/IAAFT_analyses_Gross100-1000surrogates/tab_zboot_mr_freq",end_bio,"_",end_nor,"_with",anrands,"_IAAFT_corrected_loc.csv",sep=""),sep=";",dec=".",col.names=F,row.names=F)

}else{
#We need to create a wmr object to store the values, so that it is recognized by the function image_mvcwt. We then store the values from previous analyses in it.
ref_wmr = wmr(mm)

tmp_xy=read.csv(paste("../Teich_resultsLFS/IAAFT_analyses_Gross100-1000surrogates/tab_xy_mr_freq",end_bio,"_",end_nor,"_with",anrands,"_IAAFT_corrected_loc.csv",sep=""),header=T,sep=";",dec=".")
ref_wmr$x=tmp_xy[!is.na(tmp_xy[,"x"]),"x"]
ref_wmr$y=tmp_xy[!is.na(tmp_xy[,"y"]),"y"]

tmp_z=as.matrix(read.csv(paste("../Teich_resultsLFS/IAAFT_analyses_Gross100-1000surrogates/tab_z_mr_freq",end_bio,"_",end_nor,"_with",anrands,"_IAAFT_corrected_loc.csv",sep=""),header=F,sep=";",dec="."))
tmp_array_z=array(0,dim=c(dim(tmp_z),1))
tmp_array_z[,,1]=tmp_z
ref_wmr$z=tmp_array_z

tmp_z.boot=as.matrix(read.csv(paste("../Teich_resultsLFS/IAAFT_analyses_Gross100-1000surrogates/tab_zboot_mr_freq",end_bio,"_",end_nor,"_with",anrands,"_IAAFT_corrected_loc.csv",sep=""),header=F,sep=";",dec="."))
tmp_array_z.boot=array(0,dim=c(dim(tmp_z.boot),1))
tmp_array_z.boot[,,1]=tmp_z.boot
ref_wmr$z.boot=tmp_array_z.boot
}

print(paste(Sys.time(),"before image"))
par(mar=c(5,5,3,3))
image_mvcwt_for_colormaps(ref_wmr,reset.par=F,cex.axis=4,z.fun="Mod",adj="None")
mtext("b)",side=2,line=-2,at=0.48,cex=1.5,outer=T,las=1)
mtext("Scale (years)", side=2, line=-2,at=0.28, outer = TRUE,cex=1.3)
mtext("Years",side=1, line=-2, at=0.425,outer = TRUE,cex=1.3)
print(paste(Sys.time(),"after image"))

#abline(v=2006,lwd=3,col="black") #This is supposed to change in 2006 with water management
print("After wavelet")
print(Sys.time())
dev.off()
}
