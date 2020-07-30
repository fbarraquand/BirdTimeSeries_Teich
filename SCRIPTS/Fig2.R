####CP: original script was written in 2018 to compute the Gross index for synchrony (Gross et al. 2014) on taxonomic and functional groups of birds, and to assess their significance. 
#### 2019/06/25: Corrected an error on standard deviation. They were computed under H0 (on shifted time series) and were associated to the observed value (errobars were built around the observed value with the sd from rands)
### 2019/06/27: Added BH p-val correction. Tried IAAFT and Ebisuzaki surrogates for the "between" case, and proposed to also use them for the "within" case.
### 2019/07/03: Corrected the way rare species were ignored (not taken into account the right way before). This does not change anything, actually, as the Gross function already removed them, but this is cleaner (and removes the warnings).
### 2019/09/04: Changed the script to adapt to biomasses
### 2019/09/05: Presentation details (las=1, Anas-> Anatini, etc.)
### 2019/09/12: Ducks are actually waterfowl
### 2020/07/09: Added the normalization option and 1000 rands
### 2020/07: Added saving files to plot them later without having to relaunch the whole analysis 

rm(list=ls())
graphics.off()

box_index=T #If box_index is true, draw boxplots for the whole distribution for the Gross index with shifted time series. Otherwise, we just show a line for the 5%-95% percentiles
type_correct="BH" #Was bonferroni before
amethod_b="iaaft" #Could also be shift or ebisuzaki. IAAFT is used between groups as shifting with two groups only may lead to the same surrogates (there is a limited number of combinations of shift you can do with only two time series
amethod_w="shift" #Shift is used within groups
anrands=1000 #Number of surrogates we want to test the significance
biomass=F
if(biomass){
end_bio="biomasses"
}else{
end_bio="abundances"
}

doyouload=F #True if you use files that already exist, False if you want to relaunch the analyses

normalize_seq=c(TRUE,FALSE)

source("SCRIPTS/test_synchrony_Gross.r") #Warning: test_synchrony Gross uses SCRIPTS/iaaft.r It sometimes is switched to ../../SCRIPTS/iaaft.r as it also used in other folders that are not located at the root of this folder
set.seed(42)
thresh=0.1

sp_to_ignore=c("Anas discors","Anas americana","Calidris melanotos","Calidris pusilla","Calidris ruficollis", "Calidris fuscicollis", "Calidris himantopus", "Burhinus oedicnemus","Phalaropus lobatus","Charadrius alexandrinus","Haematopus ostralegus","Calidris maritima","Aythya nyroca","Bucephala clangula","Melanitta nigra","Mergus serrator","Clangula hyemalis","Alopochen aegyptiaca", "Aix galericulata","Cygnus atratus","Tadorna ferruginea","Branta leucopsis","Anser fabalis","Anser albifrons","Cygnus cygnus","Mergus merganser","Anser brachyrhynchus") #There are species that are not abundant enough to be used. 

for(normalize in normalize_seq){
	print(paste("normalize",normalize))

	if(normalize){
		end_nor="scaled"
	}else{
		end_nor="NOTscaled"
	}

print("Anas")
print(Sys.time())
#Load data
db_warm=read.csv(paste("IN/warmseason_anas_detailed_v2_",end_bio,"_wtoutrarespecies.txt",sep=""))
db_cold=read.csv(paste("IN/coldseason_anas_detailed_v2_",end_bio,"_wtoutrarespecies.txt",sep=""))
db_warm=db_warm[,c(2,3,4)] #Choose the columns we need, that is dates, sp_data_frame and abundance
db_cold=db_cold[,c(2,3,4)]
#Right format for synchrony scripts
names(db_warm)=c("dates","sp_data_frame","abundance")
names(db_cold)=c("dates","sp_data_frame","abundance")

#Subset of the dataframe so we don't have the species we want to ignore and some faulty points in 2016. Also, divide the database between before and after 2006
db_warm_all=subset(db_warm,!(as.character(sp_data_frame) %in%sp_to_ignore) &dates<2016)
db_warm_pre_2006=subset(db_warm,!(as.character(sp_data_frame) %in%sp_to_ignore)  &dates<=2006)
db_warm_post_2006=subset(db_warm,!(as.character(sp_data_frame) %in%sp_to_ignore)  &dates>2006&dates<2016)
db_cold_all=subset(db_cold,!(as.character(sp_data_frame) %in%sp_to_ignore)  &dates<2016)
db_cold_pre_2006=subset(db_cold,!(as.character(sp_data_frame) %in%sp_to_ignore) &dates<=2006)
db_cold_post_2006=subset(db_cold,!(as.character(sp_data_frame) %in%sp_to_ignore) &dates>2006&dates<2016)

#Normalize time series if necessary
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
#In this case, we just take the value from an already existing file.
        mat_save=read.table(paste("OUT/tab_data_frame_Gross_anas_",end_bio,"_",end_nor,"_with",anrands,".csv",sep=""),sep=";",dec=".",header=T)
        list_anas=list()
        for(v in 1:nrow(mat_save)){
                list_anas[[v]]=list(obs=as.numeric(mat_save[v,"obs"]),pval=as.numeric(mat_save[v,"pval"]),alternative=as.character(mat_save[v,"alternative"]),rands=as.numeric(c(mat_save[v,grep("rands",colnames(mat_save))])))
        }

}else{

#Compute synchrony values
synch_warm_all=community_sync_Gross(db_warm_all,nrands=anrands,method=amethod_w)
synch_warm_pre_2006=community_sync_Gross(db_warm_pre_2006,nrands=anrands,method=amethod_w)
synch_warm_post_2006=community_sync_Gross(db_warm_post_2006,nrands=anrands,method=amethod_w)

synch_cold_all=community_sync_Gross(db_cold_all,nrands=anrands,method=amethod_w)
synch_cold_pre_2006=community_sync_Gross(db_cold_pre_2006,nrands=anrands,method=amethod_w)
synch_cold_post_2006=community_sync_Gross(db_cold_post_2006,nrands=anrands,method=amethod_w)

list_anas=list(synch_cold_all,synch_cold_pre_2006,synch_cold_post_2006,synch_warm_all,synch_warm_pre_2006,synch_warm_post_2006)

#Save everything in a separate text file
mat_save=matrix(NA,nrow=length(list_anas),ncol=3+anrands+1) #3 for obs, pval, alternative ; nrands for all the values... of rands. We add 1 to nrands because we also use the observed value in the computation of the pvalues
colnames(mat_save)=c("obs",paste("rands",1:(anrands+1),sep=""),"pval","alternative")
for(v in 1:length(list_anas)){
        mat_save[v,"obs"]=list_anas[[v]]$obs
        mat_save[v,"pval"]=list_anas[[v]]$pval
        mat_save[v,"alternative"]=list_anas[[v]]$alternative
        for(r in 1:anrands){
                mat_save[v,paste("rands",r,sep="")]=list_anas[[v]]$rands[r]
        }
        mat_save[v,paste("rands",anrands+1,sep="")]=list_anas[[v]]$rands[anrands+1]
}
write.table(mat_save,paste("OUT/tab_data_frame_Gross_anas_",end_bio,"_",end_nor,"_with",anrands,".csv",sep=""),sep=";",col.names=TRUE,row.names=F,dec=".")
}

#With the values in mat_save, correct the p-value for FDR
mat=rep(NA,length(list_anas))
for(v in 1:length(list_anas)){
	mat[v]=list_anas[[v]]$pval
}
mat_adj=p.adjust(mat,method=type_correct)
for(v in 1:length(list_anas)){
	list_anas[[v]]$pval=mat_adj[v]
}

print("Calidris")
print(Sys.time())
#We look at only 6 species for the Calidris, we could also look at only 4
#Load data
db_warm=read.csv(paste("IN/warmseason_calidris_detailed_v2_",end_bio,"_wtoutrarespecies.txt",sep=""))
db_cold=read.csv(paste("IN/coldseason_calidris_detailed_v2_",end_bio,"_wtoutrarespecies.txt",sep=""))
db_warm=db_warm[,c(2,3,4)]
db_cold=db_cold[,c(2,3,4)]
#Right format for synchrony scripts
names(db_warm)=c("dates","sp_data_frame","abundance")
names(db_cold)=c("dates","sp_data_frame","abundance")

#Subset of the database to divide between before and after 2006, remove the species we need to ignore and remove faulty dates in 2016.
db_warm_all=subset(db_warm,!(as.character(sp_data_frame) %in%sp_to_ignore) &dates<2016)
db_warm_pre_2006=subset(db_warm,!(as.character(sp_data_frame) %in%sp_to_ignore)  &dates<=2006)
db_warm_post_2006=subset(db_warm,!(as.character(sp_data_frame) %in%sp_to_ignore)  &dates>2006&dates<2016)
db_cold_all=subset(db_cold,!(as.character(sp_data_frame) %in%sp_to_ignore)  &dates<2016)
db_cold_pre_2006=subset(db_cold,!(as.character(sp_data_frame) %in%sp_to_ignore) &dates<=2006)
db_cold_post_2006=subset(db_cold,!(as.character(sp_data_frame) %in%sp_to_ignore) &dates>2006&dates<2016)

#Normalize if necessary
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

#Load data from a previous text file if the analyses have already been done
if(doyouload){

        mat_save=read.table(paste("OUT/tab_data_frame_Gross_calidris_",end_bio,"_",end_nor,"_with",anrands,".csv",sep=""),sep=";",dec=".",header=T)
        list_calidris=list()
        for(v in 1:nrow(mat_save)){
                list_calidris[[v]]=list(obs=as.numeric(mat_save[v,"obs"]),pval=as.numeric(mat_save[v,"pval"]),alternative=as.character(mat_save[v,"alternative"]),rands=as.numeric(c(mat_save[v,grep("rands",colnames(mat_save))])))
        }

}else{
#Launch the Gross analyses
synch_warm_all=community_sync_Gross(db_warm_all,nrands=anrands,method=amethod_w)
synch_warm_pre_2006=community_sync_Gross(db_warm_pre_2006,nrands=anrands,method=amethod_w)
synch_warm_post_2006=community_sync_Gross(db_warm_post_2006,nrands=anrands,method=amethod_w)

synch_cold_all=community_sync_Gross(db_cold_all,nrands=anrands,method=amethod_w)
synch_cold_pre_2006=community_sync_Gross(db_cold_pre_2006,nrands=anrands,method=amethod_w)
synch_cold_post_2006=community_sync_Gross(db_cold_post_2006,nrands=anrands,method=amethod_w)

list_calidris=list(synch_cold_all,synch_cold_pre_2006,synch_cold_post_2006,synch_warm_all,synch_warm_pre_2006,synch_warm_post_2006)

mat_save=matrix(NA,nrow=length(list_calidris),ncol=3+anrands+1) #3 for obs, pval, alternative ; nrands for all the values... of rands. We add 1 to nrands because we also use the observed value in the computation of the pvalues
colnames(mat_save)=c("obs",paste("rands",1:(anrands+1),sep=""),"pval","alternative")
for(v in 1:length(list_calidris)){
        mat_save[v,"obs"]=list_calidris[[v]]$obs
        mat_save[v,"pval"]=list_calidris[[v]]$pval
        mat_save[v,"alternative"]=list_calidris[[v]]$alternative
        for(r in 1:anrands){
                mat_save[v,paste("rands",r,sep="")]=list_calidris[[v]]$rands[r]
        }
        mat_save[v,paste("rands",anrands+1,sep="")]=list_calidris[[v]]$rands[anrands+1]
}
write.table(mat_save,paste("OUT/tab_data_frame_Gross_calidris_",end_bio,"_",end_nor,"_with",anrands,".csv",sep=""),sep=";",col.names=TRUE,row.names=F,dec=".")
}

#Adjust p-values	
mat=rep(NA,length(list_calidris))
for(v in 1:length(list_calidris)){
        mat[v]=list_calidris[[v]]$pval
}
mat_adj=p.adjust(mat,method=type_correct)
for(v in 1:length(list_calidris)){
        list_calidris[[v]]$pval=mat_adj[v]
}


print("Waders")
print(Sys.time())
#Load data
db_warm=read.csv(paste("IN/warmseason_waders_detailed_v2_",end_bio,"_wtoutrarespecies.txt",sep=""))
db_cold=read.csv(paste("IN/coldseason_waders_detailed_v2_",end_bio,"_wtoutrarespecies.txt",sep=""))
db_warm=db_warm[,c(2,3,4)]
db_cold=db_cold[,c(2,3,4)]
#Right format for synchrony scripts
names(db_warm)=c("dates","sp_data_frame","abundance")
names(db_cold)=c("dates","sp_data_frame","abundance")
limicoles=unique(db_warm$sp_data_frame)

#Divide the data set between before and after 2006 and remove species and dates we do not want
db_warm_all=subset(db_warm,!(as.character(sp_data_frame) %in%sp_to_ignore) &dates<2016)
db_warm_pre_2006=subset(db_warm,!(as.character(sp_data_frame) %in%sp_to_ignore)  &dates<=2006)
db_warm_post_2006=subset(db_warm,!(as.character(sp_data_frame) %in%sp_to_ignore)  &dates>2006&dates<2016)
db_cold_all=subset(db_cold,!(as.character(sp_data_frame) %in%sp_to_ignore)  &dates<2016)
db_cold_pre_2006=subset(db_cold,!(as.character(sp_data_frame) %in%sp_to_ignore) &dates<=2006)
db_cold_post_2006=subset(db_cold,!(as.character(sp_data_frame) %in%sp_to_ignore) &dates>2006&dates<2016)

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

#load resuls from previous analyses
if(doyouload){

        mat_save=read.table(paste("OUT/tab_data_frame_Gross_waders_",end_bio,"_",end_nor,"_with",anrands,".csv",sep=""),sep=";",dec=".",header=T)
        list_waders=list()
        for(v in 1:nrow(mat_save)){
                list_waders[[v]]=list(obs=as.numeric(mat_save[v,"obs"]),pval=as.numeric(mat_save[v,"pval"]),alternative=as.character(mat_save[v,"alternative"]),rands=as.numeric(c(mat_save[v,grep("rands",colnames(mat_save))])))
        }

}else{

#Compute synchrony values
synch_warm_all=community_sync_Gross(db_warm_all,nrands=anrands,method=amethod_w)
synch_warm_pre_2006=community_sync_Gross(db_warm_pre_2006,nrands=anrands,method=amethod_w)
synch_warm_post_2006=community_sync_Gross(db_warm_post_2006,nrands=anrands,method=amethod_w)

synch_cold_all=community_sync_Gross(db_cold_all,nrands=anrands,method=amethod_w)
synch_cold_pre_2006=community_sync_Gross(db_cold_pre_2006,nrands=anrands,method=amethod_w)
synch_cold_post_2006=community_sync_Gross(db_cold_post_2006,nrands=anrands,method=amethod_w)

list_waders=list(synch_cold_all,synch_cold_pre_2006,synch_cold_post_2006,synch_warm_all,synch_warm_pre_2006,synch_warm_post_2006)

mat_save=matrix(NA,nrow=length(list_waders),ncol=3+anrands+1) #3 for obs, pval, alternative ; nrands for all the values... of rands. We add 1 to nrands because we also use the observed value in the computation of the pvalues
colnames(mat_save)=c("obs",paste("rands",1:(anrands+1),sep=""),"pval","alternative")
for(v in 1:length(list_waders)){
        mat_save[v,"obs"]=list_waders[[v]]$obs
        mat_save[v,"pval"]=list_waders[[v]]$pval
        mat_save[v,"alternative"]=list_waders[[v]]$alternative
        for(r in 1:anrands){
                mat_save[v,paste("rands",r,sep="")]=list_waders[[v]]$rands[r]
        }
        mat_save[v,paste("rands",anrands+1,sep="")]=list_waders[[v]]$rands[anrands+1]
}
write.table(mat_save,paste("OUT/tab_data_frame_Gross_waders_",end_bio,"_",end_nor,"_with",anrands,".csv",sep=""),sep=";",col.names=TRUE,row.names=F,dec=".")
}

mat=rep(NA,length(list_waders))
for(v in 1:length(list_waders)){
        mat[v]=list_waders[[v]]$pval
}
#Correct p-values
mat_adj=p.adjust(mat,method=type_correct)
for(v in 1:length(list_waders)){
        list_waders[[v]]$pval=mat_adj[v]
}

print("Ducks")
print(Sys.time())
#Load data
db_warm=read.csv(paste("IN/warmseason_ducks_detailed_v2_",end_bio,"_wtoutrarespecies.txt",sep=""))
db_cold=read.csv(paste("IN/coldseason_ducks_detailed_v2_",end_bio,"_wtoutrarespecies.txt",sep=""))
db_warm=db_warm[,c(2,3,4)]
db_cold=db_cold[,c(2,3,4)]
#Right format for synchrony scripts
names(db_warm)=c("dates","sp_data_frame","abundance")
names(db_cold)=c("dates","sp_data_frame","abundance")

#Divide the data set between before and after 2006 and remove species and dates we do not want
db_warm_all=subset(db_warm,!(as.character(sp_data_frame) %in%sp_to_ignore)  &dates<2016)
db_warm_pre_2006=subset(db_warm,!(as.character(sp_data_frame) %in%sp_to_ignore)   &dates<=2006)
db_warm_post_2006=subset(db_warm,!(as.character(sp_data_frame) %in%sp_to_ignore)  &dates>2006&dates<2016)
db_cold_all=subset(db_cold,!(as.character(sp_data_frame) %in%sp_to_ignore)  &dates<2016)
db_cold_pre_2006=subset(db_cold,!(as.character(sp_data_frame) %in%sp_to_ignore) &dates<=2006)
db_cold_post_2006=subset(db_cold,!(as.character(sp_data_frame) %in%sp_to_ignore)&dates>2006&dates<2016)

#Normalize if necessary
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

#load previous analyses
if(doyouload){

        mat_save=read.table(paste("OUT/tab_data_frame_Gross_ducks_",end_bio,"_",end_nor,"_with",anrands,".csv",sep=""),sep=";",dec=".",header=T)
        list_duck=list()
        for(v in 1:nrow(mat_save)){
                list_duck[[v]]=list(obs=as.numeric(mat_save[v,"obs"]),pval=as.numeric(mat_save[v,"pval"]),alternative=as.character(mat_save[v,"alternative"]),rands=as.numeric(c(mat_save[v,grep("rands",colnames(mat_save))])))
        }

}else{

#Compute synchrony values
synch_warm_all=community_sync_Gross(db_warm_all,nrands=anrands,method=amethod_w)
synch_warm_pre_2006=community_sync_Gross(db_warm_pre_2006,nrands=anrands,method=amethod_w)
synch_warm_post_2006=community_sync_Gross(db_warm_post_2006,nrands=anrands,method=amethod_w)

synch_cold_all=community_sync_Gross(db_cold_all,nrands=anrands,method=amethod_w)
synch_cold_pre_2006=community_sync_Gross(db_cold_pre_2006,nrands=anrands,method=amethod_w)
synch_cold_post_2006=community_sync_Gross(db_cold_post_2006,nrands=anrands,method=amethod_w)

list_duck=list(synch_cold_all,synch_cold_pre_2006,synch_cold_post_2006,synch_warm_all,synch_warm_pre_2006,synch_warm_post_2006)

mat_save=matrix(NA,nrow=length(list_duck),ncol=3+anrands+1) #3 for obs, pval, alternative ; nrands for all the values... of rands. We add 1 to nrands because we also use the observed value in the computation of the pvalues
colnames(mat_save)=c("obs",paste("rands",1:(anrands+1),sep=""),"pval","alternative")
for(v in 1:length(list_waders)){
        mat_save[v,"obs"]=list_duck[[v]]$obs
        mat_save[v,"pval"]=list_duck[[v]]$pval
        mat_save[v,"alternative"]=list_duck[[v]]$alternative
        for(r in 1:anrands){
                mat_save[v,paste("rands",r,sep="")]=list_duck[[v]]$rands[r]
        }
        mat_save[v,paste("rands",anrands+1,sep="")]=list_duck[[v]]$rands[anrands+1]
}
write.table(mat_save,paste("OUT/tab_data_frame_Gross_ducks_",end_bio,"_",end_nor,"_with",anrands,".csv",sep=""),sep=";",col.names=TRUE,row.names=F,dec=".")
}

#Adjust p-values
mat=rep(NA,length(list_duck))
for(v in 1:length(list_duck)){
        mat[v]=list_duck[[v]]$pval
}
mat_adj=p.adjust(mat,method=type_correct)
for(v in 1:length(list_duck)){
        list_duck[[v]]$pval=mat_adj[v]
}

print('Actual Figure 2, within groups')
#PLOT
upsi=0.05
color=rep(c("Black","Lightblue","Dodgerblue2"),2)
filename_pdf=paste("Fig2_new3_JAE_",end_bio,"_",end_nor,".pdf",sep="")
pdf(paste("Submission_JAE/Revisions_R2/",filename_pdf,sep=""),width=11,height=14)
#if(box_index){
#pdf(paste("Submission_JAE/Revisions/Fig2_new2_JAE",type_correct,amethod_b,"boxplot.pdf",sep="_"),width=11,height=14)
#}else{
#pdf(paste("Submission_JAE/Revisions/Fig2_new2_JAE",type_correct,amethod_b,"line.pdf",sep="_"),width=11,height=14)
#}

#Plot everything within groups
par(mfrow=c(2,2),mar=c(3,4.5,2,.25),oma=c(1,2,1,.25),mgp=c(3,1,0),xpd=NA)
plot(0,0,t="n",ylim=c(-1,1),xlim=c(0,7.5),xaxt="n",xlab="",ylab="Synchrony index",cex.axis=1.8,cex.lab=1.8,main="Taxonomic groups",cex.main=2,las=1)
mtext("Within",side=2,line=0.3,outer=T,cex=2,font=2,at=0.75)
mtext("a)",side=2,line=-4,at=0.99,cex=2,outer=T,las=1)
axis(1,at=c(2,6),lab=c("Cold","Warm"),cex.axis=1.8)
for (v in 1:6){
        essai_taxo=list_anas
        obs=essai_taxo[[v]]$obs
        p_s=essai_taxo[[v]]$pval
        x=v+(v>3)
        points(x,obs,pch=21,col=color[v],bg=color[v],cex=2)
	if(box_index){
        boxplot(essai_taxo[[v]]$rands[1:anrands],at=x,add=T,boxwex=0.25,range=0,xaxt="n",yaxt="n")
	}else{
        plou=quantile(essai_taxo[[v]]$rands[1:anrands],c(0.05,0.95))
        arrows(x,plou[1],x,plou[2],angle=90,length=0.1,code=3,col=color[v],lwd=2)
	}
        if (p_s<=thresh){
                points(x,as.numeric(obs)+upsi,pch='*',col="red",cex=2)
                }

        essai_taxo=list_calidris
        obs=essai_taxo[[v]]$obs
        p_s=essai_taxo[[v]]$pval
        x=v+(v>3)
        points(x+0.2,obs,pch=22,col=color[v],bg=color[v],cex=2)
	if(box_index){
        boxplot(essai_taxo[[v]]$rands[1:anrands],at=x+0.2,add=T,boxwex=0.25,range=0,xaxt="n",yaxt="n")
	}else{
        plou=quantile(essai_taxo[[v]]$rands[1:anrands],c(0.05,0.95))
        arrows(x+0.2,plou[1],x+0.2,plou[2],angle=90,length=0.1,code=3,col=color[v],lwd=2)
	}
        if (p_s<=thresh){
                points(x+0.2,as.numeric(obs)+upsi,pch='*',col="red",cex=2)
                }

        }

lines(c(0,7.5),c(0.,0.),lty=2,lwd=2)
legend("bottomleft",c("All","Pre-2006","Post-2006"),pch=NA,fill=c("black","Lightblue","Dodgerblue2"),pt.cex=2,bty="n",cex=2)
legend("bottomright",c("Anatini","Calidris"),pch=c(21,22),pt.bg=c("black"),pt.cex=2,bty="n",cex=2)

plot(0,0,t="n",ylim=c(-1,1),xlim=c(0,7.5),xaxt="n",xlab="",ylab="",yaxt="n",main='Functional groups',cex.main=2,las=1)
mtext("b)",side=2,line=-35,at=0.99,cex=2,outer=T,las=1)

axis(1,at=c(2,6),lab=c("Cold","Warm"),cex.axis=1.8)
for (v in 1:6){
        essai_taxo=list_waders
        obs=essai_taxo[[v]]$obs
        p_s=essai_taxo[[v]]$pval
        x=v+(v>3)
        points(x,obs,pch=23,col=color[v],bg=color[v],cex=2)
	if(box_index){
        boxplot(essai_taxo[[v]]$rands[1:anrands],at=x,add=T,boxwex=0.25,range=0,xaxt="n",yaxt="n")
        }else{
	plou=quantile(essai_taxo[[v]]$rands[1:anrands],c(0.05,0.95))
        arrows(x,plou[1],x,plou[2],angle=90,length=0.1,code=3,col=color[v],lwd=2)
	}
        if (p_s<=thresh){
                points(x,as.numeric(obs)+upsi,pch='*',col="red",cex=2)
                }

        essai_taxo=list_duck
        obs=essai_taxo[[v]]$obs
        p_s=essai_taxo[[v]]$pval
        x=v+(v>3)
        points(x+0.2,obs,pch=24,col=color[v],bg=color[v],cex=2)
	if(box_index){
        boxplot(essai_taxo[[v]]$rands[1:anrands],at=x+0.2,add=T,boxwex=0.25,range=0,yaxt="n",xaxt="n")
	}else{
        plou=quantile(essai_taxo[[v]]$rands[1:anrands],c(0.05,0.95))
        arrows(x+0.2,plou[1],x+0.2,plou[2],angle=90,length=0.1,code=3,col=color[v],lwd=2)
	}
        if (p_s<=thresh){
                points(x+0.2,as.numeric(obs)+upsi,pch='*',col="red",cex=2)
                }

        }

lines(c(0,7.5),c(0.,0.),lty=2,lwd=2)
legend("bottomright",c("Waders","Waterfowl"),pch=c(23,24),pt.bg="black",pt.cex=2,bty="n",cex=2)

print('Now computing between groups')
#Anas and Calidris
##Analyses are exactly the same, apart from the method used, which is now amethod_b
print("Anas Calidris")
print(Sys.time())
db_warm=read.csv(paste("IN/warmseason_",end_bio,"_asdataframe_summed_v2_wtoutrarespecies.csv",sep=""),sep=";",header=T)
db_cold=read.csv(paste("IN/coldseason_",end_bio,"_asdataframe_summed_v2_wtoutrarespecies.csv",sep=""),sep=";",header=T)

#Right format for synchrony scripts
names(db_warm)=c("dates","sp_data_frame","abundance")
names(db_cold)=c("dates","sp_data_frame","abundance")

db_warm_all=subset(db_warm,sp_data_frame %in% c("Anas","Calidris")&dates<2016)
db_warm_pre_2006=subset(db_warm,sp_data_frame %in% c("Anas","Calidris")&dates<=2006)
db_warm_post_2006=subset(db_warm,sp_data_frame %in% c("Anas","Calidris")&dates>2006&dates<2016)
db_cold_all=subset(db_cold,sp_data_frame %in% c("Anas","Calidris")&dates<2016)
db_cold_pre_2006=subset(db_cold,sp_data_frame %in% c("Anas","Calidris")&dates<=2006)
db_cold_post_2006=subset(db_cold,sp_data_frame %in% c("Anas","Calidris")&dates>2006&dates<2016)

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

        mat_save=read.table(paste("OUT/tab_data_frame_Gross_anascalidris_",end_bio,"_",end_nor,"_with",anrands,".csv",sep=""),sep=";",dec=".",header=T)
        essai_taxo=list()
        for(v in 1:nrow(mat_save)){
                essai_taxo[[v]]=list(obs=as.numeric(mat_save[v,"obs"]),pval=as.numeric(mat_save[v,"pval"]),alternative=as.character(mat_save[v,"alternative"]),rands=as.numeric(c(mat_save[v,grep("rands",colnames(mat_save))])))
        }

}else{


#Compute synchrony values
synch_warm_all=community_sync_Gross(db_warm_all,nrands=anrands,method=amethod_b)
synch_warm_pre_2006=community_sync_Gross(db_warm_pre_2006,nrands=anrands,method=amethod_b)
synch_warm_post_2006=community_sync_Gross(db_warm_post_2006,nrands=anrands,method=amethod_b)

synch_cold_all=community_sync_Gross(db_cold_all,nrands=anrands,method=amethod_b)
synch_cold_pre_2006=community_sync_Gross(db_cold_pre_2006,nrands=anrands,method=amethod_b)
synch_cold_post_2006=community_sync_Gross(db_cold_post_2006,nrands=anrands,method=amethod_b)

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
write.table(mat_save,paste("OUT/tab_data_frame_Gross_anascalidris_",end_bio,"_",end_nor,"_with",anrands,".csv",sep=""),sep=";",col.names=TRUE,row.names=F,dec=".")
}


#Plot everything
upsi=0.05

mat=rep(NA,length(essai_taxo))
for(v in 1:length(essai_taxo)){
        mat[v]=essai_taxo[[v]]$pval
}
mat_adj=p.adjust(mat,method=type_correct)
for(v in 1:length(essai_taxo)){
        essai_taxo[[v]]$pval=mat_adj[v]
}


color=rep(c("Black","Lightblue","Dodgerblue2"),2)
plot(0,0,t="n",ylim=c(-1,1),xlim=c(0,7.5),xaxt="n",xlab="",ylab="Synchrony index",cex.lab=1.8,cex.axis=1.8,las=1)
axis(1,at=c(2,6),lab=c("Cold","Warm"),cex.axis=1.8)
mtext("Between",side=2,line=0.3,outer=T,cex=2,font=2,at=0.25)
mtext("c)",side=2,line=-4,at=0.49,cex=2,outer=T,las=1)

for (v in 1:6){
        obs=essai_taxo[[v]]$obs
        p_s=essai_taxo[[v]]$pval
        x=v+(v>3)
        points(x,obs,pch=21,col=color[v],bg=color[v],cex=2)
	if(box_index){
        boxplot(essai_taxo[[v]]$rands[1:anrands],at=x,add=T,boxwex=0.25,range=0,yaxt="n",xaxt="n")
	}else{
        plou=quantile(essai_taxo[[v]]$rands[1:anrands],c(0.05,0.95))
        arrows(x,plou[1],x,plou[2],angle=90,length=0.1,code=3,col=color[v],lwd=2)
	}
        if (p_s<=thresh){
                points(x,as.numeric(obs)+upsi,pch='*',col="red",cex=2)
                }

        }
lines(c(0.0,7.5),c(0,0),lty=2,lwd=2)
ll=c(essai_taxo[[2]]$obs,essai_taxo[[3]]$obs)
legend("bottomleft",c("All","Pre-2006","Post-2006"),pch=NA,fill=c("black","Lightblue","Dodgerblue2"),pt.cex=2,bty="n",cex=2)
legend("topleft",c("Anatini/Calidris"),pch=21,pt.bg=c("black"),pt.cex=2,bty="n",cex=2)


print("Waders Ducks")
print(Sys.time())

db_warm_all=subset(db_warm,sp_data_frame %in% c("Waders","Ducks")&dates<2016)
db_warm_pre_2006=subset(db_warm,sp_data_frame %in% c("Waders","Ducks")&dates<=2006)
db_warm_post_2006=subset(db_warm,sp_data_frame %in% c("Waders","Ducks")&dates>2006&dates<2016)
db_cold_all=subset(db_cold,sp_data_frame %in% c("Waders","Ducks")&dates<2016)
db_cold_pre_2006=subset(db_cold,sp_data_frame %in% c("Waders","Ducks")&dates<=2006)
db_cold_post_2006=subset(db_cold,sp_data_frame %in% c("Waders","Ducks")&dates>2006&dates<2016) #Because there is a problem for 2016 (na values)

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

        mat_save=read.table(paste("OUT/tab_data_frame_Gross_wadersducks_",end_bio,"_",end_nor,"_with",anrands,".csv",sep=""),sep=";",dec=".",header=T)
        essai_taxo=list()
        for(v in 1:nrow(mat_save)){
                essai_taxo[[v]]=list(obs=as.numeric(mat_save[v,"obs"]),pval=as.numeric(mat_save[v,"pval"]),alternative=as.character(mat_save[v,"alternative"]),rands=as.numeric(c(mat_save[v,grep("rands",colnames(mat_save))])))
        }

}else{


#Compute synchrony values
synch_warm_all=community_sync_Gross(db_warm_all,nrands=anrands,method=amethod_b)
synch_warm_pre_2006=community_sync_Gross(db_warm_pre_2006,nrands=anrands,method=amethod_b)
synch_warm_post_2006=community_sync_Gross(db_warm_post_2006,nrands=anrands,method=amethod_b)

synch_cold_all=community_sync_Gross(db_cold_all,nrands=anrands,method=amethod_b)
synch_cold_pre_2006=community_sync_Gross(db_cold_pre_2006,nrands=anrands,method=amethod_b)
synch_cold_post_2006=community_sync_Gross(db_cold_post_2006,nrands=anrands,method=amethod_b)

#Plot everything
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
write.table(mat_save,paste("OUT/tab_data_frame_Gross_wadersducks_",end_bio,"_",end_nor,"_with",anrands,".csv",sep=""),sep=";",col.names=TRUE,row.names=F,dec=".")
}


mat=rep(NA,length(essai_taxo))
for(v in 1:length(essai_taxo)){
        mat[v]=essai_taxo[[v]]$pval
}
mat_adj=p.adjust(mat,method=type_correct)
for(v in 1:length(essai_taxo)){
        essai_taxo[[v]]$pval=mat_adj[v]
}


color=rep(c("Black","Lightblue","Dodgerblue2"),2)
plot(0,0,t="n",ylim=c(-1,1),xlim=c(0,7.5),xaxt="n",xlab="",ylab="",yaxt="n",las=1)
axis(1,at=c(2,6),lab=c("Cold","Warm"),cex.axis=2)
mtext("d)",side=2,line=-35,at=0.49,cex=2,outer=T,las=1)
for (v in 1:6){
        obs=essai_taxo[[v]]$obs
        p_s=essai_taxo[[v]]$pval
        x=v+(v>3)
        points(x,obs,pch=22,col=color[v],bg=color[v],cex=2)
	if(box_index){
	boxplot(essai_taxo[[v]]$rands[1:anrands],at=x,add=T,boxwex=0.25,range=0,yaxt="n",xaxt="n")
	}else{
	plou=quantile(essai_taxo[[v]]$rands[1:anrands],c(0.05,0.95)) 
        arrows(x,plou[1],x,plou[2],angle=90,length=0.1,code=3,col=color[v],lwd=2)
	}
        if (p_s<=thresh){
                points(x,as.numeric(obs)+upsi,pch='*',col="red",cex=2)
                }

        }
lines(c(0.0,7.5),c(0,0),lty=2,lwd=2)

legend("topleft",c("Waterfowl/Waders"),pch=c(22),pt.bg=c("black"),pt.cex=2,bty="n",cex=2)
dev.off()
}
