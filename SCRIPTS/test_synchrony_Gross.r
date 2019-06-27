##CP 16/05/2018
#Using the synchrony package, and more particularly the community.sync function from Gouhier and Guichard (2014) to compute confidence intervals for Gross index
#CP 26/06/2019: Maybe there was an error in the way we computed pval=rands<obs/(nrands+1). Is the +1 ok?

source("SCRIPTS/iaaft.R")

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
	#sp_data_frame=sapply(colnames(data_matrix),function(x) rep(x,dim(data_matrix)[1]))
	#adates=rep(c(data[2]),dim(data_matrix)[2])
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
        if (alternative == "two-tailed") {
            results$pval = sum(abs(results$rands) >= abs(results$obs))/(nrands) 
               # +1)
        }else if (alternative == "greater"){
            results$pval = sum(results$rands >= results$obs)/(nrands)
               # 1)
        }else{results$pval = sum(results$rands <= results$obs)/(nrands) }
           # 1)}
#        results$rands[nrands + 1] = results$obs #I guess this came from a previously written function (not by me), but this does not make sense à rands and obs may not represent the same hypothesis.
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
