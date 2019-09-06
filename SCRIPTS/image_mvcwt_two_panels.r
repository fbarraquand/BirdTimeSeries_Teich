#### First used by CP (2018) from Keitt's package mvcwt
#### 2019/07/03 Corrected BY to BH correction for homogeneity reeaon
#### CP 2019/09/05 Used las=1 to have horizontal y-axis tick labels

image_mvcwt_two_panels=function (x, z.fun = "Re", bound = 1, reset.par = TRUE,...) 
{
    z.fun = match.fun(z.fun)
#    layout(matrix(c(1,2),nrow=1,ncol=2),widths=c(10,1),heights=1)
#    opar = par(no.readonly = TRUE)
#    if (reset.par) 
#        on.exit(par(opar))
    pal = colorRampPalette(rev(brewer.pal(11, "Spectral")))(1024)
    with(x, {
        nvar = ifelse(length(dim(z)) == 3, dim(z)[3], 1)
#	par(mar=c(2, 2, 2, 2))#,oma=c(2,3,0,0))
#        par(mfrow = c(nvar, 1), mar = rep(0.2, 4), oma = rep(5, 
#            4))
        for (i in 1:nvar) {
            image(x, y, z.fun(z[, , i]), log = "y", col = pal, 
                axes = FALSE, xlab="",ylab="",...)
#		mtext("b)",side=2,line=-2,at=0.975,cex=2,outer=T,las=1)
#		lines(c(2006,2006),c(0.2,20),lwd=2,col="black") 
            if (i%%2) 
                axis(2,cex.axis=1.5,las=1)
            else axis(4)
            if (exists("z.boot") && !is.null(z.boot)) {
                z.boot = 1 - abs(1 - 2 * z.boot)
#                contour(x, y, z.boot[, , i], levels = 0.05, lty = 3,  #Removing the uncorrected p-values
#                  add = TRUE, drawlabels = FALSE)
                zb = p.adjust(as.vector(z.boot), method = "BH")
                dim(zb) = dim(z.boot)
                contour(x, y, zb[, , i], levels = 0.1, lwd = 0.5,lty=1, 
                  add = TRUE, drawlabels = FALSE)
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
#        mtext("Years", 1, 1, outer = TRUE,cex=2)
#        mtext("Scale (years)", 2, 1, outer = TRUE,cex=2)

		mtext("b)",side=2,line=-2,at=0.45,cex=1.5,outer=T,las=1)
        mtext("Years", 1, 3,cex=1.5)
        mtext("Scale (years)", 2, 3,cex=1.5)
#        mtext("Scale (years)", 2, 2.5,cex=1.5)


    })
	par(mar=c(2,.5,2,3))
	seqi=seq(min(z.fun(x$z[, , 1]),na.rm=T),max(z.fun(x$z[, , 1]),na.rm=T),length.out=1024)
	image(x=1,y=seqi,matrix(seqi,nrow=1,ncol=1024),col=pal,axes=FALSE,xlab="",ylab="")
axis(4,cex.axis=1.5,las=1)
box()

    return(invisible(x))
}
