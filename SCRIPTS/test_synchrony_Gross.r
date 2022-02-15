##CP 16/05/2018
#Using the synchrony package, and more particularly the community.sync function from Gouhier and Guichard (2014) to compute confidence intervals for Gross index
#CP 04/07/2019: p-value is now relying on quantile. But we wan't compute FDR-correction if we only have the corresponding boolean. Back to the very first formulation we had: https://stats.stackexchange.com/questions/25927/doubling-the-tails-in-two-sample-permutation-test (answer from caracal) 

#source("../../../SCRIPTS/iaaft.R")
#source("SCRIPTS/iaaft.R")

community_sync_Gross=function (data, nrands = 0, alternative = c("two-tailed","greater","less"),method=c("shift","iaaft","ebisuzaki"), 
    quiet = FALSE, ...) 
{
	alternative="two-tailed" #Because I'm lazy
    alternatives = c("greater", "less","two-tailed")
    alternative = match.arg(tolower(alternative), alternatives)
#    data = as.matrix(data)
    results = list()
    results$obs = community_sync_Gross_aux(data)
    if (nrands > 0) {
	data_matrix=convert_to_matrix(data)
	dates=data$dates
	sp_data_frame=data$sp_data_frame
	nr = NROW(data_matrix)
        nc = NCOL(data_matrix)
        if (!quiet) 
            prog.bar = txtProgressBar(min = 0, max = nrands, 
                style = 3)
        results$rands = numeric(length = nrands) * NA
	for (i in 1:nrands) {
		if(method=="shift"){
                	lags = sample(1:nr, size = nc, replace = TRUE)
                	rand.mat = mlag(data_matrix, lags)
			rand.mat=data.frame(sp_data_frame,dates,c(rand.mat))
		}else{
			rand.mat=matrix(NA,nr,nc)
			iter=0
			for(sp in unique(as.character(sp_data_frame))){
				iter=iter+1
				sp_tmp=data_matrix[,sp]
				if(method=="iaaft"){
					rand.mat[,iter]=iaaft_surrogate(sp_tmp)
				}else if(method=="ebisuzaki"){
					rand.mat[,iter]=ebisuzaki_surrogate(sp_tmp)
				 }
				
			}
			rand.mat=data.frame(sp_data_frame,dates,c(rand.mat))
		}
		names(rand.mat)=c("sp_data_frame","dates","abundance")
            results$rands[i] = community_sync_Gross_aux(rand.mat)
            if (!quiet) 
                setTxtProgressBar(prog.bar, i)
        }
        results$rands[nrands + 1] = results$obs #CP: See North et al. (2002) and answers (2003). CP & FB agreed to keep North et al.'s formula (and Gouhier et al. too BTW). 
        if (alternative == "two-tailed") {
		results$pval= 2*min(sum(results$rands >= results$obs),sum(results$rands <= results$obs))/(nrands+1) #p2 in Submission_JAE/Revisions/pvalue 
        }else if(alternative=="greater"){
		results$pval = sum(results$rands >= results$obs)/(nrands+1)
	}else if(alternative=="less"){
		results$pval = sum(results$rands <= results$obs)/(nrands+1)
	}else{
		stop("Change your alternative")
	}
        results$alternative = alternative
    }
    return(results)
}

community_sync_Gross_aux=function (data){
	require('codyn')
	res=synchrony(data,time.var="dates",species.var="sp_data_frame",abundance.var="abundance",metric="Gross")
	return(res)
}

convert_to_matrix=function(data){
	sp=sort(unique(data$sp_data_frame))
	adates=sort(unique(data$dates))
	mm=matrix(0,ncol=length(sp),nrow=length(adates))
	colnames(mm)=sp
	rownames(mm)=adates
	for (s in sp){
		for (d in adates){
			if(length(which(data$sp_data_frame==s&data$dates==d))>0){
			mm[as.character(d),s]=data[which(data$sp_data_frame==s&data$dates==d),"abundance"]
			}
		}
	}
	return(mm)
}

mlag <- function (x, lag) {
  len=NROW(x)
  ncol=NCOL(x)
  index=matrix(rep(0:(len-1), ncol), nrow=len) 
  lag=lag %% len
  index=(index-matrix(rep(lag, each=len), nrow=len)) %% len + 1
  result=matrix(nrow=len, ncol=ncol, NA)
  for (i in 1:ncol)
    result[, i]=x[index[, i], i] 
  colnames(result)= 
  return(result)
}
