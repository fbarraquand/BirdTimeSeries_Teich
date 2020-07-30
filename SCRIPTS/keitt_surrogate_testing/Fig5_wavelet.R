# Script by CPicoche 2018
#The aim is to compare synchrony indices at different taxonomic levels
#06/03/2020: Just checking the time necessary for 100 and 1000 randomization
#01/07/2020: Checking the effect of log and standardization

rm(list=ls())
graphics.off()
source("SCRIPTS/test_synchrony_Gross.r")
library('mvcwt')
source("Submission_JAE/Revisions_R2/Simulations_response/image_mvcwt_for_colormaps.r")

library("RColorBrewer")

set.seed(42)

thresh=0.1
type_correct="BH" #was Bonferroni before
log_b=F
amethod="iaaft"
anrands=1000
#normalize_seq=c(T,F)
normalize_seq=c(T)
doyouload=F

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

db_warm_all=subset(db_warm,sp_data_frame %in% c("Cormorant","HeronEgret")&dates<2016)
db_warm_pre_2006=subset(db_warm,sp_data_frame %in% c("Cormorant","HeronEgret")&dates<=2006)
db_warm_post_2006=subset(db_warm,sp_data_frame %in% c("Cormorant","HeronEgret")&dates>2006&dates<2016)#Because using the NA for 2016 makes results wrong
db_cold_all=subset(db_cold,sp_data_frame %in% c("Cormorant","HeronEgret")&dates<2016)
db_cold_pre_2006=subset(db_cold,sp_data_frame %in% c("Cormorant","HeronEgret")&dates<=2006)
db_cold_post_2006=subset(db_cold,sp_data_frame %in% c("Cormorant","HeronEgret")&dates>2006&dates<2016)

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

#Compute synchrony values
if(doyouload){

        mat_save=read.table(paste("OUT/tab_data_frame_Gross_triad_",end_bio,"_",end_nor,"_with",anrands,".csv",sep=""),sep=";",dec=".",header=T)
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
write.table(mat_save,paste("OUT/tab_data_frame_Gross_triad_",end_bio,"_",end_nor,"_with",anrands,".csv",sep=""),sep=";",col.names=TRUE,row.names=F,dec=".")

}

print("After Gross")
print(Sys.time())

#if(1==0){ ##remove that after
#Plot everything
mat=rep(NA,length(essai_taxo))
for(v in 1:length(essai_taxo)){
        mat[v]=essai_taxo[[v]]$pval
}
mat_adj=p.adjust(mat,method=type_correct)
for(v in 1:length(essai_taxo)){
        essai_taxo[[v]]$pval=mat_adj[v]
}


color=rep(c("Black","Lightblue","Dodgerblue2"),2)
pdf(paste("Submission_JAE/Revisions_R2/triad_synchrony_2panels_",end_bio,"_",end_nor,"_",end_log,"_test_with",anrands,"rand_nocorrection_smallgrid.pdf",sep=""),width=9,height=9) 
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

db=read.csv(paste("IN/summed_",end_bio,"_v2_wtoutrarespecies.csv",sep=""),sep=";",header=T)
db$Date=as.Date(db$Date)

db=subset(db,Nom_latin %in% c("Cormorant","HeronEgret"))
db_tmp=db
names(db_tmp)=c("sp_data_frame","dates","abundance")
db_av=subset(db_tmp,dates<as.Date("2007-01-01"))
db_ap=subset(db_tmp,dates>=as.Date("2007-01-01"))
db_tmp$dates=as.numeric(db_tmp$dates)
db_av$dates=as.numeric(db_av$dates)
db_ap$dates=as.numeric(db_ap$dates)

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
mm=mvcwt(x,tab,min.scale=0.2,max.scale=10.0,nscales=100,loc=regularize(x,nsteps=ceiling(length(x)/2)))

if(!doyouload){
#This function computes the wavelet ratio of the whole community (see Keitt's paper in 2008)
ref_wmr = wmr.boot(mm, smoothing = 1,reps=anrands)
ref_wmr$x=ref_wmr$x+year_min #Change the dates to be "human-readable"

if(length(ref_wmr$x)>length(ref_wmr$y)){
        yy=c(ref_wmr$y,rep(NA,length(ref_wmr$x)-length(ref_wmr$y)))
        xx=ref_wmr$x
}else{
        xx=c(ref_wmr$x,rep(NA,length(ref_wmr$y)-length(ref_wmr$x)))
        yy=ref_wmr$y
}
tab_xy=cbind(xx,yy)
colnames(tab_xy)=c("x","y")
write.table(tab_xy,paste("OUT/tab_xy_mr_triad",end_bio,"_",end_nor,"_with",anrands,"rands.csv",sep=""),sep=";",dec=".",col.names=T,row.names=F)

tab_z=ref_wmr$z
write.table(as.matrix(tab_z[,,1]),paste("OUT/tab_z_mr_triad",end_bio,"_",end_nor,"_with",anrands,"rands.csv",sep=""),sep=";",dec=".",col.names=F,row.names=F)

tab_z.boot=ref_wmr$z.boot
write.table(as.matrix(tab_z.boot[,,1]),paste("OUT/tab_zboot_mr_triad",end_bio,"_",end_nor,"_with",anrands,"rands.csv",sep=""),sep=";",dec=".",col.names=F,row.names=F)

}else{
ref_wmr = wmr(mm)

tmp_xy=read.csv(paste("OUT/tab_xy_mr_triad",end_bio,"_",end_nor,"_with",anrands,"rands.csv",sep=""),header=T,sep=";",dec=".")
ref_wmr$x=tmp_xy[!is.na(tmp_xy[,"x"]),"x"]
ref_wmr$y=tmp_xy[!is.na(tmp_xy[,"y"]),"y"]

tmp_z=as.matrix(read.csv(paste("OUT/tab_z_mr_triad",end_bio,"_",end_nor,"_with",anrands,"rands.csv",sep=""),header=F,sep=";",dec="."))
tmp_array_z=array(0,dim=c(dim(tmp_z),1))
tmp_array_z[,,1]=tmp_z
ref_wmr$z=tmp_array_z

tmp_z.boot=as.matrix(read.csv(paste("OUT/tab_zboot_mr_triad",end_bio,"_",end_nor,"_with",anrands,"rands.csv",sep=""),header=F,sep=";",dec="."))
tmp_array_z.boot=array(0,dim=c(dim(tmp_z.boot),1))
tmp_array_z.boot[,,1]=tmp_z.boot
ref_wmr$z.boot=tmp_array_z.boot

}

par(mar=c(3,5,2,3))
image_mvcwt_for_colormaps(ref_wmr,reset.par=F,cex.axis=4,z.fun="Mod",adj="None")
mtext("b)",side=2,line=-2,at=0.48,cex=1.5,outer=T,las=1)

#abline(v=2006,lwd=3,col="black") #This is supposed to change in 2006 with water management
print("After wavelet")
print(Sys.time())
dev.off()
} 
#
