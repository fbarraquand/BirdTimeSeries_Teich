#### First used by CP (2018) from Keitt's package mvcwt
#### CP 2019/07/03 Corrected BY to BH correction for homogeneity reason
#### CP 2019/09/05 Used las=1 to have horizontal y-axis tick labels
#### CP 2020/07/03 Trying out different colormaps

image_mvcwt_for_colormaps=function (x, z.fun = "Re", bound = 1, reset.par = TRUE,col.palette="Spectral",inv=F,...) 
{
    require("viridis")
	test_col=9 #Should be 11
    z.fun = match.fun(z.fun)
#    layout(matrix(c(1,2),nrow=1,ncol=2),widths=c(10,1),heights=1)
#    opar = par(no.readonly = TRUE)
#    if (reset.par) 
#        on.exit(par(opar))
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
#	par(mar=c(2, 2, 2, 2),oma=c(2,3,0,0))
#        par(mfrow = c(nvar, 1), mar = rep(0.2, 4), oma = rep(5, 
#            4))
        for (i in 1:nvar) {
            image(x, y, z.fun(z[, , i]), log = "y", col = pal,#main=col.palette, 
                axes = FALSE, xlab="",ylab="",...)
		print("it's ok here")
#		mtext("b)",side=2,line=-2,at=0.975,cex=2,outer=T,las=1)
#		lines(c(2006,2006),c(0.2,20),lwd=2,col="black") 
            if (i%%2) 
                axis(2,cex.axis=1.8,las=1)
            else axis(4)
            if (exists("z.boot") && !is.null(z.boot)) {


#                z.boot = 1 - abs(1 - 2 * z.boot)
		z.boot_tmp=1-2*z.boot
#		print(z.boot_tmp[171:173,200:212,1])
#		image(z.boot_tmp[,,1])
		z.boot_pos=array(NA,dim(z.boot))
		z.boot_neg=array(NA,dim(z.boot))
		z.boot_pos[which(z.boot_tmp>0,arr.ind=T)]=1-abs(1-2*z.boot[which(z.boot_tmp>0,arr.ind=T)])
		z.boot_neg[which(z.boot_tmp<0,arr.ind=T)]=1-abs(1-2*z.boot[which(z.boot_tmp<0,arr.ind=T)])

#		print("Pos")
#		print(z.boot_pos[171:173,200:212,1])
#		print("Neg")
#		print(z.boot_neg[171:173,200:212,1])

#                contour(x, y, z.boot[, , i], levels = 0.05, lty = 3,  #Removing the uncorrected p-values
#                  add = TRUE, drawlabels = FALSE)
                #zb = p.adjust(as.vector(z.boot), method = "BH")
                #dim(zb) = dim(z.boot)
		#print(which(zb[,,i]>0.95))
		#print(which(zb[,,i]<0.05))


                #contour(x, y, zb[, , i], levels = 0.1, lwd = 0.5,lty=1, 
                #  add = TRUE, drawlabels = FALSE)

                zb_pos = p.adjust(as.vector(z.boot_pos), method = "BH")
                dim(zb_pos) = dim(z.boot_pos)
                contour(x, y, zb_pos[, , i], levels = 0.05, lwd = 1.5,lty=1, 
                  add = TRUE, drawlabels = FALSE,col="darkblue")
                zb_neg = p.adjust(as.vector(z.boot_neg), method = "BH")
                dim(zb_neg) = dim(z.boot_neg)
                contour(x, y, zb_neg[, , i], levels = 0.05, lwd = 1.5,lty=1, 
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
        axis(1,cex.axis=1.8)
        mtext("Years", 1, 1, outer = TRUE,cex=2)
        mtext("Scale (years)", 2, 1, outer = TRUE,cex=2)
    })
#WWW	par(mar=c(2,.5,2,3))
	seqi=seq(min(z.fun(x$z[, , 1]),na.rm=T),max(z.fun(x$z[, , 1]),na.rm=T),length.out=1024)
	image(x=1,y=seqi,matrix(seqi,nrow=1,ncol=1024),col=pal,axes=FALSE,xlab="",ylab="")
#WWW	image(x=1,y=seqi,matrix(seqi,nrow=1,ncol=1024),col=pal,axes=FALSE,xlab="",ylab="")
#WWW axis(4,cex.axis=1.8,las=1)
#WWW box()

    return(invisible(x))
}
