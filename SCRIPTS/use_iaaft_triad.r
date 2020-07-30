# Script by CPicoche 2018
#The aim is to compare synchrony indices at different taxonomic levels for Cormoran/HeronEgret. This is based on Fig5.R, then adapted to use iaaft instead of Keitt's classical randomization.

rm(list=ls())
graphics.off()
source("SCRIPTS/test_synchrony_Gross.r")
library('mvcwt')
source("SCRIPTS/image_mvcwt_for_colormaps.r")

library("RColorBrewer")

set.seed(42)

thresh=0.1
type_correct="BH" #was Bonferroni before, used only for the Gross index
log_b=F #We checked the effect of using the logarithm of data before
amethod="iaaft" #method used for the Gross index. We do not use the shift method here because there are only two groups and the combination of shifts which can be performed with only two time series is limited
anrands=1000 #Number of surrogates
normalize_seq=c(T,F)
doyouload=T #True if we want to take analyses that have been performed previously, False if we want to relaunch analyses

biomass=F
if(biomass){
end_bio="biomasses"
}else{
end_bio="abundances"
}

for(normalize in normalize_seq){
if(normalize){
end_nor="scaled"
}else{
end_nor="NOTscaled"
}

# Grand Cormoran Phalacrocorax carbo
# Héron cendré Ardea cinerea
# Aigrette garzette Egretta garzetta

#Load data
db_warm=read.csv(paste("IN/warmseason_",end_bio,"_asdataframe_summed_v2_wtoutrarespecies.csv",sep=""),sep=";",header=T)
db_cold=read.csv(paste("IN/coldseason_",end_bio,"_asdataframe_summed_v2_wtoutrarespecies.csv",sep=""),sep=";",header=T)

#Right format for synchrony scripts
names(db_warm)=c("dates","sp_data_frame","abundance")
names(db_cold)=c("dates","sp_data_frame","abundance")

if(log_b){
        db_warm$abundance=log(db_warm$abundance+1)
        db_cold$abundance=log(db_cold$abundance+1)
	end_log="log"
}else{
	end_log="NOlog"
}

#Subset to separate database before and after 2006 + remove problematic data
db_warm_all=subset(db_warm,sp_data_frame %in% c("Cormorant","HeronEgret")&dates<2016)
db_warm_pre_2006=subset(db_warm,sp_data_frame %in% c("Cormorant","HeronEgret")&dates<=2006)
db_warm_post_2006=subset(db_warm,sp_data_frame %in% c("Cormorant","HeronEgret")&dates>2006&dates<2016)#Because using the NA for 2016 makes results wrong
db_cold_all=subset(db_cold,sp_data_frame %in% c("Cormorant","HeronEgret")&dates<2016)
db_cold_pre_2006=subset(db_cold,sp_data_frame %in% c("Cormorant","HeronEgret")&dates<=2006)
db_cold_post_2006=subset(db_cold,sp_data_frame %in% c("Cormorant","HeronEgret")&dates>2006&dates<2016)

#Normalize species by species
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

#Load synchrony values
if(doyouload){

        mat_save=read.table(paste("../Teich_resultsLFS/IAAFT_analyses_Gross100-1000surrogates/tab_data_frame_Gross_triad_",end_bio,"_",end_nor,"_with",anrands,".csv",sep=""),sep=";",dec=".",header=T)
        essai_taxo=list()
        for(v in 1:nrow(mat_save)){
                essai_taxo[[v]]=list(obs=as.numeric(mat_save[v,"obs"]),pval=as.numeric(mat_save[v,"pval"]),alternative=as.character(mat_save[v,"alternative"]),rands=as.numeric(c(mat_save[v,grep("rands",colnames(mat_save))])))
        }


}else{

#Compute synchrony values

print("Before Gross")
print(Sys.time())
synch_warm_all=community_sync_Gross(db_warm_all,nrands=anrands,method=amethod)
synch_warm_pre_2006=community_sync_Gross(db_warm_pre_2006,nrands=anrands,method=amethod)
synch_warm_post_2006=community_sync_Gross(db_warm_post_2006,nrands=anrands,method=amethod)

synch_cold_all=community_sync_Gross(db_cold_all,nrands=anrands,method=amethod)
synch_cold_pre_2006=community_sync_Gross(db_cold_pre_2006,nrands=anrands,method=amethod)
synch_cold_post_2006=community_sync_Gross(db_cold_post_2006,nrands=anrands,method=amethod)

essai_taxo=list(synch_cold_all,synch_cold_pre_2006,synch_cold_post_2006,synch_warm_all,synch_warm_pre_2006,synch_warm_post_2006)

#Save synchrony values
mat_save=matrix(NA,nrow=length(essai_taxo),ncol=3+anrands+1) #3 for obs, pval, alternative ; nrands for all the values... of rands. We add 1 to nrands because we also use the observed value in the computation of the pvalues
colnames(mat_save)=c("obs",paste("rands",1:(anrands+1),sep=""),"pval","alternative")
for(v in 1:length(essai_taxo)){
        mat_save[v,"obs"]=essai_taxo[[v]]$obs
        mat_save[v,"pval"]=essai_taxo[[v]]$pval
        mat_save[v,"alternative"]=essai_taxo[[v]]$alternative
        for(r in 1:anrands){
                mat_save[v,paste("rands",r,sep="")]=essai_taxo[[v]]$rands[r]
        }
        mat_save[v,paste("rands",anrands+1,sep="")]=essai_taxo[[v]]$rands[anrands+1]
}
write.table(mat_save,paste("../Teich_resultsLFS/IAAFT_analyses_Gross100-1000surrogates/tab_data_frame_Gross_triad_",end_bio,"_",end_nor,"_with",anrands,".csv",sep=""),sep=";",col.names=TRUE,row.names=F,dec=".")

}

print("After Gross")
print(Sys.time())

#Adjust p-value for FDR
mat=rep(NA,length(essai_taxo))
for(v in 1:length(essai_taxo)){
        mat[v]=essai_taxo[[v]]$pval
}
mat_adj=p.adjust(mat,method=type_correct)
for(v in 1:length(essai_taxo)){
        essai_taxo[[v]]$pval=mat_adj[v]
}


#Plot everything
color=rep(c("Black","Lightblue","Dodgerblue2"),2)
pdf(paste("Submission_JAE/Revisions_R2/triad_synchrony_2panels_",end_bio,"_",end_nor,"_",end_log,"_test_with",anrands,"rand_nocorrection_smallgrid_IAAFT.pdf",sep=""),width=9,height=9) 
layout(matrix(c(1,1,2,3),nrow=2,ncol=2,byrow=T),widths=c(10,2))

#par(mfrow=c(1,1),mar=c(3,3.5,2,.25),oma=c(1,2,1,.25),mgp=c(3,1,0),xpd=NA)
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
        boxplot(essai_taxo[[v]]$rands[1:anrands],at=x,add=T,boxwex=0.25,range=0,yaxt="n",xaxt="n")
        if (p_s<thresh){
                points(x,as.numeric(obs),pch='*',col="red",cex=2)
                }

        }
lines(c(0,7.5),c(0,0),lty=2,lwd=2)

legend("bottomleft",c("All","Pre-2006","Post-2006"),pch=NA,fill=c("black","Lightblue","Dodgerblue2"),pt.cex=2,bty="n",cex=1.5)


############################################ WAVELETS ########################################

#Load data
db=read.csv(paste("IN/summed_",end_bio,"_v2_wtoutrarespecies.csv",sep=""),sep=";",header=T)
db$Date=as.Date(db$Date)

db=subset(db,Nom_latin %in% c("Cormorant","HeronEgret"))

#Build a table with the right format to be analysed
dates=unique(db$Date)
tab=matrix(0,nrow=length(dates),ncol=2)
colnames(tab)=c("Cormorant","HeronEgret")
rownames(tab)=dates
for(id in 1:length(dates)){
        for (s in c('Cormorant',"HeronEgret")){
                id_d=which(db$Date==dates[id]&db$Nom_latin==s)
                if(length(id_d)>0){
			if(log_b){
                        	tab[id,s]=log(db$Nombre[id_d]+1)
			}else{
                        	tab[id,s]=db$Nombre[id_d]
			}
                }
        }
}

if(normalize){
	tab[,"Cormorant"]=scale(tab[,'Cormorant'])
	tab[,"HeronEgret"]=scale(tab[,'HeronEgret'])
}

x=(dates-dates[1])/365.25
year_min=1981
#This function computes the Morlet wavelet transform for each bird species separately
print("Before wavelet")
print(Sys.time())

seq_x=seq(1,365.25*35+3*30.5,length.out=423)/365.25

mm=mvcwt(x,tab,min.scale=mean(diff(seq_x))*3,max.scale=10.0,nscales=100,loc=seq_x)

if(!doyouload){
#This function computes the wavelet ratio of the whole community (see Keitt's paper in 2008)
ref_wmr=wmr(mm)
ref_wmr$x=ref_wmr$x+year_min #Change the dates to be "human-readable"
ref_val=ref_wmr$z[,,1]

#Computes iaaft surrogates and the corresponding wavelet modulus ratio for each of them
tab_values_iaaft=array(NA,dim=c(length(mm$x),length(mm$y),anrands+1))
tab_values_iaaft[,,anrands+1]=ref_val
prog.bar = txtProgressBar(min = 0, max = anrands,style = 3)
for(i in 1:anrands){
        setTxtProgressBar(prog.bar, i)
        tab_tmp=tab
        tab_tmp[,"Cormorant"]=iaaft_surrogate(tab[,'Cormorant'])
        tab_tmp[,"HeronEgret"]=iaaft_surrogate(tab[,'HeronEgret'])
        #mmtmp=mvcwt(x,tab_tmp,min.scale=0.2,max.scale=10.0,nscales=100,loc=regularize(x,nsteps=length(x)/2))
	mmtmp=mvcwt(x,tab_tmp,min.scale=mean(diff(seq_x))*3,max.scale=10.0,nscales=100,loc=seq_x)
        wmr_tmp=wmr(mmtmp)
        tab_values_iaaft[,,i]=wmr_tmp$z[,,1]
}

#Computes Pr(X<=x) for each observed wavelet modulus ratio
tab_pval=array(NA,dim=c(length(mm$x),length(mm$y),1))
for(i in 1:length(mm$x)){
        for(j in 1:length(mm$y)){
                #tab_pval[i,j,1]= 2*min(sum(tab_values_iaaft[i,j,] >= ref_val[i,j]),sum(tab_values_iaaft[i,j,] < ref_val[i,j]))/(nrep+1) #This line could be used if we wanted to output the p-value directly. Here we stick to the philosophy of the mvcwt package and output Pr(X<=x)
                tab_pval[i,j,1]= sum(tab_values_iaaft[i,j,] <= ref_val[i,j])/(anrands+1)
#                if(tab_pval[i,j,1]>1){stop()}

        }
}


ref_wmr$z.boot=tab_pval

#There is not the same number of locations and scales so we need to transform a bit the vectors in order to bind them
if(length(ref_wmr$x)>length(ref_wmr$y)){
        yy=c(ref_wmr$y,rep(NA,length(ref_wmr$x)-length(ref_wmr$y)))
        xx=ref_wmr$x
}else{
        xx=c(ref_wmr$x,rep(NA,length(ref_wmr$y)-length(ref_wmr$x)))
        yy=ref_wmr$y
}
tab_xy=cbind(xx,yy)
colnames(tab_xy)=c("x","y")
write.table(tab_xy,paste("../Teich_resultsLFS/IAAFT_analyses_Gross100-1000surrogates/tab_xy_mr_triad",end_bio,"_",end_nor,"_with",anrands,"_IAAFT.csv",sep=""),sep=";",dec=".",col.names=T,row.names=F)

tab_z=ref_wmr$z
write.table(as.matrix(tab_z[,,1]),paste("../Teich_resultsLFS/IAAFT_analyses_Gross100-1000surrogates/tab_z_mr_triad",end_bio,"_",end_nor,"_with",anrands,"_IAAFT.csv",sep=""),sep=";",dec=".",col.names=F,row.names=F)

tab_z.boot=ref_wmr$z.boot
write.table(as.matrix(tab_z.boot[,,1]),paste("../Teich_resultsLFS/IAAFT_analyses_Gross100-1000surrogates/tab_zboot_mr_triad",end_bio,"_",end_nor,"_with",anrands,"_IAAFT.csv",sep=""),sep=";",dec=".",col.names=F,row.names=F)


}else{
#We need to create a wmr object to store the values, so that it is recognized by the function image_mvcwt. We then store the values from previous analyses in it.
ref_wmr = wmr(mm)

tmp_xy=read.csv(paste("../Teich_resultsLFS/IAAFT_analyses_Gross100-1000surrogates/IAAFT_analyses_Gross100-1000surrogates/tab_xy_mr_triad",end_bio,"_",end_nor,"_with",anrands,"_IAAFT.csv",sep=""),header=T,sep=";",dec=".")
ref_wmr$x=tmp_xy[!is.na(tmp_xy[,"x"]),"x"]
ref_wmr$y=tmp_xy[!is.na(tmp_xy[,"y"]),"y"]

tmp_z=as.matrix(read.csv(paste("../Teich_resultsLFS/IAAFT_analyses_Gross100-1000surrogates/tab_z_mr_triad",end_bio,"_",end_nor,"_with",anrands,"_IAAFT.csv",sep=""),header=F,sep=";",dec="."))
tmp_array_z=array(0,dim=c(dim(tmp_z),1))
tmp_array_z[,,1]=tmp_z
ref_wmr$z=tmp_array_z

tmp_z.boot=as.matrix(read.csv(paste("../Teich_resultsLFS/IAAFT_analyses_Gross100-1000surrogates/tab_zboot_mr_triad",end_bio,"_",end_nor,"_with",anrands,"_IAAFT.csv",sep=""),header=F,sep=";",dec="."))
tmp_array_z.boot=array(0,dim=c(dim(tmp_z.boot),1))
tmp_array_z.boot[,,1]=tmp_z.boot
ref_wmr$z.boot=tmp_array_z.boot

}

par(mar=c(5,5,2,3))
image_mvcwt_for_colormaps(ref_wmr,reset.par=F,cex.axis=4,z.fun="Mod",adj="None")
mtext("b)",side=2,line=-2,at=0.48,cex=1.5,outer=T,las=1)
mtext("Scale (years)", side=2, line=-2,at=0.28, outer = TRUE,cex=1.3)
mtext("Years",side=1, line=-2, at=0.425,outer = TRUE,cex=1.3)

#abline(v=2006,lwd=3,col="black") #This is supposed to change in 2006 with water management
print("After wavelet")
print(Sys.time())
dev.off()
} 
