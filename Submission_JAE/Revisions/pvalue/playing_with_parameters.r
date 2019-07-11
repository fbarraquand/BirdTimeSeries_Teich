###CP 11/07/2019 Illustration for the pvalue discussion on Monte-Carlo distributions

rm(list=ls())
graphics.off()

set.seed(666)

pval1=function(obs,plou){
	num=sum(abs(plou)>=abs(obs))+1
	den=length(plou)+1
	pval=num/den
	return(pval)
}

pval2=function(obs,plou){
	num1=(sum(plou>=obs)+1)/(length(plou)+1) #Removing the +1 does not change the pval, at least graphically
	num2=(sum(plou<=obs)+1)/(length(plou)+1)
	pval=2*min(num1,num2)
	return(pval)
}

obs=rbeta(1000,shape1=2,shape2=5)-0.6
obs[1]=0.25 #Should be rejected
obs[2]=-0.25 #Should be kept under H0

plou=rbeta(10000,shape1=2,shape2=5)-0.6
pdf("distrib.pdf")
dens=density(plou)
plot(dens,main="",xlab=expression(eta))
a_min=quantile(plou,0.05)
a_max=quantile(plou,0.95)
abline(v=c(a_min,a_max),lty=2)
abline(v=0,lwd=1.5)
abline(v=obs[1],col="red")
abline(v=obs[2],col="blue")
text(paste("pval1",format(pval1(obs[1],plou),digits=3)),x=0.1,y=2,col="red")
text(paste("pval2",format(pval2(obs[1],plou),digits=3)),x=0.1,y=1.9,col="red")
text(paste("pval1",format(pval1(obs[2],plou),digits=3)),x=0.1,y=1.7,col="blue")
text(paste("pval2",format(pval2(obs[2],plou),digits=3)),x=0.1,y=1.6,col="blue")
dev.off()

pval_1=rep(NA,length(obs))
pval_2=rep(NA,length(obs))
for(i in 1:length(obs)){
	pval_1[i]=pval1(obs[i],plou)
	pval_2[i]=pval2(obs[i],plou)
}

pdf("results_about_pvalues.pdf",width=10,height=4.5)
par(mfrow=c(1,3))
plot(0,0,t="n",xlab="",main="Prop pval<0.05",xlim=c(0.5,2.5),ylim=c(0,100),ylab="")
boxplot(100*sum(pval_1<0.05)/length(obs),at=1,add=T,names="pval1")
boxplot(100*sum(pval_2<0.05)/length(obs),at=2,add=T,names="pval2")
abline(h=5)
idx=sort(obs,index.return=T)$ix
plot(obs[idx],pval_1[idx],t="l",main="pval1",ylab="",xlab=expression(eta))
points(density(plou)$x,density(plou)$y/max(density(plou)$y))
plot(obs[idx],pval_2[idx],t="l",main="pval2",ylab="",xlab=expression(eta))
points(density(plou)$x,density(plou)$y/max(density(plou)$y))
legend("topright",c("pval","density distrib"),pch=c(NA,1),lty=c(1,NA),bty="n")
dev.off()
