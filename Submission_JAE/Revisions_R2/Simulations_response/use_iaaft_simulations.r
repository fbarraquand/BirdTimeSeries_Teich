# Script by CPicoche 2018
#The aim is to compare synchrony indices at different taxonomic levels
#06/03/2020: Just checking the time necessary for 100 and 1000 randomization
#01/07/2020: Checking the effect of log and standardization

rm(list=ls())
graphics.off()
source("../../../SCRIPTS/test_synchrony_Gross.r")
library('mvcwt')
source("../../../SCRIPTS/image_mvcwt_for_colormaps.r") 
library("RColorBrewer")
library("lubridate")

set.seed(42)

thresh=0.1
type_correct="BH" #was Bonferroni before
amethod="shift"
anrands=1000
biomass=F
if(biomass){
end_bio="biomasses"
}else{
end_bio="abundances"
}

doyouload=F #False if you want to launch the analyses again ; TRUE if you just want to do the plots

normalize_seq=c(TRUE,FALSE)

end_of_file_seq=c("4sp_alpha2_beta4","4sp_alpha4_beta2")
explain=c("compensation","synchrony")

for(e in 1:length(end_of_file_seq)){
end_of_file=end_of_file_seq[e]
#############################
tab_bm=read.table(paste("MockData_SAD_",end_of_file,".csv",sep=""),sep=";",dec=".",header=T)

rep=1

tab_bm=subset(tab_bm,Rep==rep)

tab_bm$Species=as.character(tab_bm$Species)
tab_bm$Time=as.numeric(as.character(tab_bm$Time))
tab_bm$Abundance=as.numeric(as.character(tab_bm$Abundance))

dates_tmp=unique(tab_bm$Time)

#Convert 1,2,3 to month and year

dates=rep(NA,length(dates_tmp))
yy=0
for(d in 1:length(dates)){
        if(dates_tmp[d]%%12==1){
                yy=yy+1
        }
        mm=dates_tmp[d]%%12
        if(mm==0){mm=12}
        tmp=as.Date(paste(yy,mm,"15",sep="-"))
        dates[d]=tmp
}

dates=as.Date(dates,origin="01-01-01")
sp=unique(tab_bm$Species)

tab_species=matrix(0,nrow=length(dates),ncol=length(sp))
colnames(tab_species)=sp

for(id in 1:length(dates)){
        for (s in sp){
                id_d=which(tab_bm$Time==dates_tmp[id]&tab_bm$Species==s)
                if(length(id_d)>0){
                        tab_species[id,s]=tab_bm$Abundance[id_d]
                }
        }
}

yy_seq=unique(year(dates))
tab_annual_warm=as.data.frame(matrix(NA,nrow=length(yy_seq)*length(sp),ncol=3))
tab_annual_cold=as.data.frame(matrix(NA,nrow=length(yy_seq)*length(sp),ncol=3))
names(tab_annual_warm)=names(tab_annual_cold)=c("dates","sp_data_frame","abundance")

for(normalize in normalize_seq){
print(paste("normalize",normalize))
if(normalize){
end_nor="scaled"
}else{
end_nor="NOTscaled"
}
tab_species_tmp=tab_species
for(s in sp){
if(normalize){
        tab_species_tmp[,s]=scale(tab_species[,s])
}

}

x=0
for(s in sp){
for(y in 1:length(yy_seq)){
        x=x+1
        subset=tab_species_tmp[year(dates)==yy_seq[y],s]
        tab_annual_warm[x,"dates"]=yy_seq[y]
        tab_annual_warm[x,"sp_data_frame"]=s
        tab_annual_warm[x,"abundance"]=mean(subset[7:10])
        tab_annual_cold[x,"dates"]=yy_seq[y]
        tab_annual_cold[x,"sp_data_frame"]=s
        tab_annual_cold[x,"abundance"]=mean(subset[1:4])
}
}

if(doyouload){

        mat_save=read.table(paste("tab_data_frame_Gross_simulated_data_",end_nor,"_with",anrands,"_",end_of_file,".csv",sep=""),sep=";",dec=".",header=T)
        essai_taxo=list()
        for(v in 1:nrow(mat_save)){
                essai_taxo[[v]]=list(obs=as.numeric(mat_save[v,"obs"]),pval=as.numeric(mat_save[v,"pval"]),alternative=as.character(mat_save[v,"alternative"]),rands=as.numeric(c(mat_save[v,grep("rands",colnames(mat_save))])))
        }

}else{



print("Before Gross")
print(Sys.time())
synch_warm_all=community_sync_Gross(tab_annual_warm,nrands=anrands,method=amethod)

synch_cold_all=community_sync_Gross(tab_annual_cold,nrands=anrands,method=amethod)

essai_taxo=list(synch_cold_all,synch_warm_all)

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
write.table(mat_save,paste("tab_data_frame_Gross_simulated_data_",end_nor,"_with",anrands,"_",end_of_file,".csv",sep=""),sep=";",col.names=TRUE,row.names=F,dec=".")
}



print("After Gross")
print(Sys.time())

mat=rep(NA,length(essai_taxo))
for(v in 1:length(essai_taxo)){
        mat[v]=essai_taxo[[v]]$pval
}
mat_adj=p.adjust(mat,method=type_correct)
for(v in 1:length(essai_taxo)){
        essai_taxo[[v]]$pval=mat_adj[v]
}


color=rep(c("Black","Lightblue","Dodgerblue2"),2)

pdf(paste("wavelet_simu_",end_nor,"_",end_of_file,"_with",anrands,"_two_panels_no_correction_smallergrid_IAAFT.pdf",sep=""),width=10,height=10)

layout(matrix(c(1,1,2,3),nrow=2,ncol=2,byrow=T),widths=c(10,2))

#par(mfrow=c(1,1),mar=c(3,3.5,2,.25),oma=c(1,2,1,.25),mgp=c(3,1,0),xpd=NA)
par(mar=c(4,5,2,3))

plot(0,0,t="n",ylim=c(-1.,1.0),xlim=c(0,3),xaxt="n",xlab="",ylab=expression(eta),cex.lab=1.5,cex.axis=1.5,las=1,main=explain[e])
mtext("a)",side=2,line=-2,at=0.96,cex=1.5,outer=T,las=1)
axis(1,at=c(1,2),lab=c("Cold","Warm"),cex.axis=1.5)
for (v in 1:2){
        obs=essai_taxo[[v]]$obs
        print(obs)
        p_s=essai_taxo[[v]]$pval
        points(v,obs,pch=21,col=color[v],bg=color[v],cex=2)
        boxplot(essai_taxo[[v]]$rands[1:100],at=v,add=T,boxwex=0.25,range=0,yaxt="n",xaxt="n")
        if (p_s<thresh){
                points(v,as.numeric(obs),pch='*',col="red",cex=2)
                }

        }
lines(c(0,7.5),c(0,0),lty=2,lwd=2)

x=(dates-dates[1])/365.25
#Here, I don't need to regularize data as they are already cleaned in the lines before

#This function computes the Morlet wavelet transform for each bird species separately
print(Sys.time())
seq_x=seq(1,365.25*35,length.out=420)/365.25

mm=mvcwt(seq_x,tab_species_tmp,min.scale=mean(diff(seq_x))*3,max.scale=10.0,nscales=100,loc=seq_x) ####WARNING. THIS ONLY WORKS BECAUSE x IS ARTIFICIAL)

print(paste(Sys.time(),"after mvcwt"))

#This function computes the wavelet ratio of the whole community (see Keitt's paper in 2008)
if(!doyouload){

ref_wmr=wmr(mm)
ref_val=ref_wmr$z[,,1]

tab_values_iaaft=array(NA,dim=c(length(mm$x),length(mm$y),anrands+1))
tab_values_iaaft[,,anrands+1]=ref_val
prog.bar = txtProgressBar(min = 0, max = anrands,style = 3)
for(i in 1:anrands){
        setTxtProgressBar(prog.bar, i)
        tab_tmp=tab_species_tmp
        for(s in 1:ncol(tab_species)){
        tab_tmp[,s]=iaaft_surrogate(tab_species_tmp[,s])
        }
	mmtmp=mvcwt(x,tab_tmp,min.scale=mean(diff(seq_x))*3,max.scale=10.0,nscales=100,loc=seq_x)
        wmr_tmp=wmr(mmtmp)
        tab_values_iaaft[,,i]=wmr_tmp$z[,,1]
}

tab_pval=array(NA,dim=c(length(mm$x),length(mm$y),1))
for(i in 1:length(mm$x)){
        for(j in 1:length(mm$y)){
#                tab_pval[i,j,1]= 2*min(sum(tab_values_iaaft[i,j,] >= ref_val[i,j]),sum(tab_values_iaaft[i,j,] < ref_val[i,j]))/(anrands+1)
                tab_pval[i,j,1]= sum(tab_values_iaaft[i,j,] <= ref_val[i,j])/(anrands+1)
                if(tab_pval[i,j,1]>1){stop()}

        }
}
ref_wmr$z.boot=tab_pval

if(length(ref_wmr$x)>length(ref_wmr$y)){
	yy=c(ref_wmr$y,rep(NA,length(ref_wmr$x)-length(ref_wmr$y)))
	xx=ref_wmr$x
}else{
	xx=c(ref_wmr$x,rep(NA,length(ref_wmr$y)-length(ref_wmr$x)))
	yy=ref_wmr$y
}
tab_xy=cbind(xx,yy)
colnames(tab_xy)=c("x","y")
write.table(tab_xy,paste("tab_xy_mr_simulated_data_",end_nor,"_with",anrands,"_",end_of_file,"_IAAFT.csv",sep=""),sep=";",dec=".",col.names=T,row.names=F)

tab_z=ref_wmr$z
write.table(as.matrix(tab_z[,,1]),paste("tab_z_mr_simulated_data_",end_nor,"_with",anrands,"_",end_of_file,"_IAAFT.csv",sep=""),sep=";",dec=".",col.names=F,row.names=F)

tab_z.boot=ref_wmr$z.boot
write.table(as.matrix(tab_z.boot[,,1]),paste("tab_zboot_mr_simulated_data_",end_nor,"_with",anrands,"_",end_of_file,"_IAAFT.csv",sep=""),sep=";",dec=".",col.names=F,row.names=F)

}else{
ref_wmr=wmr(mm)

tmp_xy=read.csv(paste("tab_xy_mr_simulated_data_",end_nor,"_with",anrands,"_",end_of_file,"_IAAFT.csv",sep=""),header=T,sep=";",dec=".")
ref_wmr$x=tmp_xy[!is.na(tmp_xy[,"x"]),"x"]
ref_wmr$y=tmp_xy[!is.na(tmp_xy[,"y"]),"y"]


tmp_z=as.matrix(read.csv(paste("tab_z_mr_simulated_data_",end_nor,"_with",anrands,"_",end_of_file,"_IAAFT.csv",sep=""),header=F,sep=";",dec="."))
tmp_array_z=array(0,dim=c(dim(tmp_z),1))
tmp_array_z[,,1]=tmp_z
ref_wmr$z=tmp_array_z

tmp_z.boot=as.matrix(read.csv(paste("tab_zboot_mr_simulated_data_",end_nor,"_with",anrands,"_",end_of_file,"_IAAFT.csv",sep=""),header=F,sep=";",dec="."))
tmp_array_z.boot=array(0,dim=c(dim(tmp_z.boot),1))
tmp_array_z.boot[,,1]=tmp_z.boot
ref_wmr$z.boot=tmp_array_z.boot
}


print(paste(Sys.time(),"before image"))
par(mar=c(3,5,2,3))
image_mvcwt_for_colormaps(ref_wmr,reset.par=F,cex.axis=4,z.fun="Mod",adj="None")
mtext("b)",side=2,line=-2,at=0.48,cex=1.5,outer=T,las=1)
print(paste(Sys.time(),"after image"))

#abline(v=2006,lwd=3,col="black") #This is supposed to change in 2006 with water management
print("After wavelet")
print(Sys.time())
dev.off()
}
}
