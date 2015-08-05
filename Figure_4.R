rm(list = ls())


library(fields)
library(e1071) 

load("Xu_MASTER.RData")


recordings=c("WTP4_1","WTP4_2","WTP4_3","WTP4_4","WTP4_5","WTP4_6","WTP4_7","WTP4_8","WTP4_9","WTP4_10","WTP4_11","WTP4_12","WTP4_13","B2TGP4_1","B2TGP4_2","B2TGP4_3","B2TGP4_4",
"B2TGP4_5","B2TGP4_6","B2TGP4_7","B2TGP4_8","B2TGP4_9","B2TGP4_10","B2TGP4_11","B2TGP4_12","B2TGP4_13","B2TGP4_14","B2TGP4_15","B2TGP4_16","B2TGP4_17")

WT=subset(a,subset=(Phenotype=="WT"),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))
distance_WT=as.numeric(names(WT_split))

B2=subset(a,subset=(Phenotype=="B2TG"),select=c("distance","TC"))
B2_split=split(B2,as.factor(B2$distance))
distance_B2=as.numeric(names(B2_split))
range_wt=range(c(WT$TC,B2$TC))

Means=mat.or.vec(30,32)
Sd=mat.or.vec(30,32)
Skew=mat.or.vec(30,32)
num_points=mat.or.vec(30,32)

D=subset(a,subset=(Phenotype=="WT"),select=c("distance","TC"))
means_w=aggregate(TC~distance,D,mean)
l_w=aggregate(TC~distance,D,function(x) mean(x)-sd(x))
u_w=aggregate(TC~distance,D,function(x) mean(x)+sd(x))
D=subset(a,subset=(Phenotype=="B2TG"),select=c("distance","TC"))
means_b=aggregate(TC~distance,D,mean)
l_b=aggregate(TC~distance,D,function(x) mean(x)-sd(x))
u_b=aggregate(TC~distance,D,function(x) mean(x)+sd(x))

for(i in 1:30){

	D=subset(a,subset=(Animal==recordings[i]),select=c("distance","TC"))

	D_split=split(D,as.factor(D$distance))

	for(j in 1:32){

	Means[i,j]=mean(D_split[[j]][,2])
	Sd[i,j]=sd(D_split[[j]][,2])
	Skew[i,j]=skewness(D_split[[j]][,2])
	num_points[i,j]=length(D_split[[j]][,2])
	
	}
}	


distance=as.numeric(names(D_split))




cat("ok to here")
pdf("Figure_4.pdf")
attach(mtcars)
par(oma=c(1,1,2,1))
layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = TRUE))
par(mar=c(4,4,1,1.5))

max_x=max(c(means_w[,1]))
plot(as.vector(means_w[,1]),as.vector(means_w[,2]),type="p",pch=20,col="black",ylim=c(0,0.75),main="",xlab="",ylab="STTC",xlim=c(-5,max_x+10),cex=1,yaxt="n",cex.lab=1,cex.axis=1,xaxt="n")
axis(2,las=2,at=c(0,0.25,0.5,0.75),cex.axis=1)
axis(1,at=c(0,200,400,600,800),labels=FALSE)
#mtext(side = 2, "STTC", line = 2.75,cex=0.75)
arrows(as.numeric(as.vector(means_w[,1])),as.numeric(as.vector(l_w[,2])), as.numeric(as.vector(means_w[,1])),as.numeric(as.vector(u_w[,2])), angle=90, code=3, length=0.03,lwd=0.5)
points(as.vector(means_b[,1]),as.vector(means_b[,2]),col="blue",pch=20,cex=1)
arrows(as.numeric(as.vector(means_b[,1])),as.numeric(as.vector(l_b[,2])), as.numeric(as.vector(means_b[,1])),as.numeric(as.vector(u_b[,2])), angle=90, code=3, length=0.03,lwd=0.5,col="blue")

legend("topright",col=c("black","blue"),pch=20,cex=0.75,legend=c("WT       (0.17 Hz, n=13)",expression(paste(beta,"2(TG) (0.21 Hz, n=17)",sep=""))))



WT_split=split(WT,as.factor(WT$distance))
WT_box=rbind(WT_split[[1]],WT_split[[2]],WT_split[[4]],WT_split[[7]],WT_split[[10]],WT_split[[14]],WT_split[[19]],WT_split[[24]])
B2_split=split(B2,as.factor(B2$distance))
B2_box=rbind(B2_split[[1]],B2_split[[2]],B2_split[[4]],B2_split[[7]],B2_split[[10]],B2_split[[14]],B2_split[[19]],B2_split[[24]])


bplot(WT_box$TC,by=WT_box$distance,pos=c(0,100,200,300,400,500,600,700),xaxt="n",pch=17, outcol="black", ylim=range_wt,ylab="STTC",boxwex=0.2,xlim=c(0,800),cex=0.5)
axis(1,at=c(0,200,400,600,800),labels=FALSE)
bplot(B2_box$TC,by=B2_box$distance,pos=c(25,125,225,325,425,525,625,725),add=TRUE,col=c("blue"),pch=17,outcol="blue",cex=0.5,xaxt="n",boxwex=0.2)


plot(distance,num_points[1,],pch=18,ylim=range(num_points),xaxt="n",xlab="",ylab="number of data points")
axis(1,at=c(0,200,400,600,800),labels=FALSE)
for(i in 2:13){
points(distance,num_points[i,],pch=18)
}
for(i in 14:30){
points(distance,num_points[i,],pch=18,col="blue")
}


plot(distance,Means[1,],pch=18,ylim=range(Means),type="l",xaxt="n",xlab="", ylab="mean",xlim=c(0,800))
axis(1,at=c(0,200,400,600,800),labels=FALSE)
for(i in 2:13){
points(distance,Means[i,],pch=18,type="l")
}
for(i in 14:30){
points(distance,Means[i,],pch=18,col="blue",type="l")
}
legend(x=650,y=0.6,pch=17,lty=1,col=c("black","blue"),legend=c("WT","B2"))

plot(distance,Sd[1,],pch=18,ylim=range(Sd,na.rm=TRUE),type="l",ylab="standard deviation",xlab=expression(paste("Intercell distance (",mu,"m)",sep="")))
for(i in 2:13){
points(distance,Sd[i,],pch=18,type="l")
}
for(i in 14:30){
points(distance,Sd[i,],pch=18,col="blue",type="l")
}


plot(distance,Skew[1,],pch=18,ylim=range(Skew,na.rm=TRUE),type="l",xlab=expression(paste("Intercell distance (",mu,"m)",sep="")),ylab="skewness",xlim=c(0,800))
for(i in 2:13){
points(distance,Skew[i,],pch=18,type="l")
}
for(i in 14:30){
points(distance,Skew[i,],pch=18,col="blue",type="l")
}



text(grconvertX(c(0.04, 0.51,0.04,0.51,0.04,0.51), from='ndc'),
   grconvertY(c(0.96, 0.96,0.65,0.65,0.33,0.33), from='ndc'),  c('A',  'B',  'C','D','E','F'), xpd=NA, cex=2.2, font=2)



dev.off()






