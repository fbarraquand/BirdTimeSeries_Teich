#### First used by CP (2018) from Keitt's package mvcwt
#### CP 2019/07/03 Corrected BY to BH correction for homogeneity reason
#### CP 2019/09/05 Used las=1 to have horizontal y-axis tick labels
#### CP 2020/07/03 Trying out different colormaps
#### CP 2020/07 Was used previously to plot the values of the index with iaaft surrogates. Not used any more as plots use image_mvcwt_for_colormaps

image_mvcwt_with_iaaft=function (x, z.fun = "Re", bound = 1, reset.par = TRUE,col.palette="Spectral",inv=F,amain="",adj="BH",...) 
{
    require("viridis")
	test_col=9 #Should be 11
    z.fun = match.fun(z.fun)
	if(col.palette=="Viridis"){
    		pal = viridis(1024) #Test
	}else if(col.palette=="Spectral"){
	    pal = colorRampPalette(rev(brewer.pal(test_col, "Spectral")))(1024) #classical
	}else if(col.palette=="Inferno"){
	    pal = inferno(1024)
	}else if(col.palette=="Magma"){
	    pal = magma(1024)
	}else if(col.palette=="Plasma"){
	    pal = plasma(1024)
	}else if(col.palette=="Oranges"){
	    pal = colorRampPalette(rev(brewer.pal(test_col, "Oranges")))(1024)
	}else if(col.palette=="PuBu"){
	    pal = colorRampPalette(rev(brewer.pal(test_col, "PuBu")))(1024)
	}else if(col.palette=="YlGnBu"){
	    pal = colorRampPalette(rev(brewer.pal(test_col, "YlGnBu")))(1024)
	}else if(col.palette=="BrBG"){
	    pal = colorRampPalette(rev(brewer.pal(test_col, "BrBG")))(1024)
	}
	if(inv){
		pal=rev(pal)
	}
    with(x, {
        nvar = ifelse(length(dim(z)) == 3, dim(z)[3], 1)
        for (i in 1:nvar) { ####MAIN IMAGE, with actual observed values for Keitt's index
            image(x, y, z.fun(z[, , i]), log = "y", col = pal,
                axes = FALSE, xlab="",ylab="",main=amain,cex.main=1.5,...)
            if (i%%2) 
                axis(2,cex.axis=1.5,las=1)
            else axis(4) #### NOW, CONTOURS FOR PVALUES
            if (exists("z.boot") && !is.null(z.boot)) {
		#tab_= 2*min(sum(tab_values_iaaft[i,j,] >= ref_val[i,j]),sum(tab_values_iaaft[i,j,] < ref_val[i,j]))/(nrep+1)
                z.boot_tmp = 1 - abs(1 - 2 * z.boot) #First compute all p-values
		if(adj=="BH"){
			zb = p.adjust(as.vector(z.boot_tmp), method = "BH") #Adjust all p-values
		}else if(adj=="BY"){
			zb = p.adjust(as.vector(z.boot_tmp), method = "BY") #Adjust all p-values
		}else if (adj=="None"){
			zb=z.boot_tmp
                }
		dim(zb) = dim(z.boot)
	
		pval=z.boot #Now, keep Pr(X<=x)
		oneminuspval=1-z.boot #Pr(X>x)
		zb_small=array(NA,dim(z.boot_tmp))
		zb_big=array(NA,dim(z.boot_tmp))
		zb_small[which(pval<oneminuspval,arr.ind=T)]=zb[which(pval<oneminuspval,arr.ind=T)] #zb_small contains the p-value corresponding to small values of x, i.e. Pr(X<=x) < Pr(X>x)
		zb_big[which(pval>oneminuspval,arr.ind=T)]=zb[which(pval>oneminuspval,arr.ind=T)] #zb_big contains the p-value corresponding to higher values of x, i.e., Pr(X<=x) >Pr(X>x)
                contour(x, y, zb_small[, , i], levels = 0.1, lwd = 1.5,lty=1, 
                  add = TRUE, drawlabels = FALSE,col="blue")
                contour(x, y, zb_big[, , i], levels = 0.1, lwd = 1.5,lty=1, 
                  add = TRUE, drawlabels = FALSE,col="red")

            }
            if (is.finite(bound)) {
                lines(min(x) + bound * y, y, lty = 2, lwd = 4, 
                  col = "darkgrey")
                lines(max(x) - bound * y, y, lty = 2, lwd = 4, 
                  col = "darkgrey")
            }
            box()
        }
        axis(1,cex.axis=1.5)
        mtext("Years", 1, 1, outer = TRUE,cex=2)
        mtext("Scale (years)", 2, 1, outer = TRUE,cex=2)
    })
	par(mar=c(2,.5,2,4))
	seqi=seq(min(z.fun(x$z[, , 1]),na.rm=T),max(z.fun(x$z[, , 1]),na.rm=T),length.out=1024)
	image(x=1,y=seqi,matrix(seqi,nrow=1,ncol=1024),col=pal,axes=FALSE,xlab="",ylab="")
#WWW	image(x=1,y=seqi,matrix(seqi,nrow=1,ncol=1024),col=pal,axes=FALSE,xlab="",ylab="")
 axis(4,cex.axis=1.5,las=1)
#WWW box()

    return(invisible(x))
}
