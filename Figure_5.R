rm(list = ls())


library(fields)
library(e1071) 

load("Bshp11_MASTER.RData")


recordings=c("WT_1","WT_2","WT_3","WT_4","WT_5","cx45_1","cx45_2","cx45_3","cx45_4","cx3645_1","cx3645_2","cx3645_3","cx3645_4","cx3645_5","cx3645_6")

WT=subset(a,subset=(Phenotype=="WT"),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))
distance_WT=as.numeric(names(WT_split))

KO=subset(a,subset=(Phenotype=="cx45KO"),select=c("distance","TC"))
KO_split=split(KO,as.factor(KO$distance))
distance_KO=as.numeric(names(KO_split))

DKO=subset(a,subset=(Phenotype=="cx3645dKO"),select=c("distance","TC"))
DKO_split=split(DKO,as.factor(DKO$distance))
distance_DKO=as.numeric(names(DKO_split))

range_wt=range(c(WT$TC,KO$TC,DKO$TC))

Means=mat.or.vec(15,32)
Sd=mat.or.vec(15,32)
Skew=mat.or.vec(15,32)
num_points=mat.or.vec(15,32)
for(i in 1:15){

	D=subset(a,subset=(Animal==recordings[i]),select=c("distance","TC"))

	D_split=split(D,as.factor(D$distance))

	for(j in 1:length(names(D_split))){

	Means[i,j]=mean(D_split[[j]][,2])
	Sd[i,j]=sd(D_split[[j]][,2])
	Skew[i,j]=skewness(D_split[[j]][,2])
	num_points[i,j]=length(D_split[[j]][,2])
	}
}	


distance=as.numeric(names(D_split))




cat("ok to here")
pdf("Figure_5.pdf")
attach(mtcars)
par(oma=c(1,1,2,1))
layout(matrix(c(1,1,2,2,1,1,2,2,3,3,4,4,3,3,4,4,5,5,6,7,5,5,8,9), 6, 4, byrow = TRUE))
par(mar=c(4,4,1,1.5))


WT_split=split(WT,as.factor(WT$distance))
WT_box=rbind(WT_split[[1]],WT_split[[2]],WT_split[[4]],WT_split[[7]],WT_split[[10]],WT_split[[14]],WT_split[[19]],WT_split[[24]])
KO_split=split(KO,as.factor(KO$distance))
KO_box=rbind(KO_split[[1]],KO_split[[2]],KO_split[[4]],KO_split[[7]],KO_split[[10]],KO_split[[14]],KO_split[[19]],KO_split[[24]])
DKO_split=split(DKO,as.factor(DKO$distance))
DKO_box=rbind(DKO_split[[1]],DKO_split[[2]],DKO_split[[4]],DKO_split[[7]],DKO_split[[10]],DKO_split[[14]],DKO_split[[19]],DKO_split[[24]])



bplot(WT_box$TC,by=WT_box$distance,pos=c(0,100,200,300,400,500,600,700),xaxt="n",main="",pch=17, outcol="black",col="white", ylim=range_wt,boxwex=0.2,xlim=c(0,800),ylab="STTC")
axis(1,at=c(0,200,400,600,800),labels=FALSE)
bplot(KO_box$TC,by=KO_box$distance,pos=c(20,120,220,320,420,520,620,720),add=TRUE,col=c("red"),pch=17,outcol="red",cex=0.4,xaxt="n",boxwex=0.2,xlim=c(0,800))
bplot(DKO_box$TC,by=DKO_box$distance,pos=c(40,140,240,340,440,540,640,740),add=TRUE,col=c("blue"),pch=17,outcol="blue",cex=0.4,xaxt="n",boxwex=0.2,xlim=c(0,800))


plot(distance,num_points[1,],pch=18,ylim=range(num_points),xaxt="n",xlab="",ylab="number of data points")
axis(1,at=c(0,200,400,600,800),labels=FALSE)
for(i in 2:5){
points(distance,num_points[i,],pch=18)
}
for(i in 6:9){
points(distance,num_points[i,],pch=18,col="red")
}
for(i in 10:15){
points(distance,num_points[i,],pch=18,col="blue")
}

plot(distance,Means[1,],pch=18,ylim=range(Means),type="l",xaxt="n",xlab="",ylab="mean",xlim=c(0,800))
axis(1,at=c(0,200,400,600,800),labels=FALSE)
for(i in 2:5){
points(distance,Means[i,],pch=18,type="l")
}
for(i in 6:9){
points(distance,Means[i,],pch=18,col="red",type="l")
}
for(i in 10:15){
points(distance,Means[i,],pch=18,col="blue",type="l")
}
legend("topright",lty=1,col=c("black","red","blue"),legend=c("wild type","Cx45ko","Cx36/45dko"))

plot(distance,Sd[1,],pch=18,ylim=range(Sd,na.rm=TRUE),type="l",ylab="standard deviation",xlab=expression(paste("Intercell distance (",mu,"m)",sep="")))
for(i in 2:5){
points(distance,Sd[i,],pch=18,type="l")
}
for(i in 6:9){
points(distance,Sd[i,],pch=18,col="red",type="l")
}
for(i in 10:15){
points(distance,Sd[i,],pch=18,col="blue",type="l")
}

plot(distance,Skew[1,],pch=18,ylim=range(Skew,na.rm=TRUE),type="l",xlab=expression(paste("Intercell distance (",mu,"m)",sep="")),ylab="skewness",xlim=c(0,800))
for(i in 2:5){
points(distance,Skew[i,],pch=18,type="l")
}
for(i in 6:9){
points(distance,Skew[i,],pch=18,col="red",type="l")
}
for(i in 10:15){
points(distance,Skew[i,],pch=18,col="blue",type="l")
}

par(mar=c(2.5,4.3,2,0))
D=subset(a,subset=(Animal==recordings[5]),select=c("distance","TC"))
D_split=split(D,as.factor(D$distance))
plot(density(D_split[[11]][,2]),type="l",main=expression(paste("WT 5, 412",mu,"m",sep="")),xlab="",ylab="",xlim=c(-0.2,0.4),xaxt="n",cex.main=0.85)
axis(1,at=c(-0.2,0,0.2,0.4),labels=FALSE)

D=subset(a,subset=(Animal==recordings[7]),select=c("distance","TC"))
D_split=split(D,as.factor(D$distance))
plot(density(D_split[[15]][,2]),type="l",main=expression(paste("Cx45ko 2, 510",mu,"m",sep="")),col="red",xlab="",ylab="",xlim=c(-0.05,0.4),xaxt="n",cex.main=0.85)
axis(1,at=c(0,0.2,0.4),labels=FALSE)

D=subset(a,subset=(Animal==recordings[15]),select=c("distance","TC"))
D_split=split(D,as.factor(D$distance))
plot(density(D_split[[13]][,2]),type="l",main=expression(paste("Cx36/45dko 6, 447",mu,"m",sep="")),col="blue",xlab="",ylab="density",xlim=c(-0.2,0.4),xaxt="n",cex.main=0.85)
axis(1,at=c(-0.2,0,0.2,0.4),labels=TRUE)

D=subset(a,subset=(Animal==recordings[9]),select=c("distance","TC"))
D_split=split(D,as.factor(D$distance))
plot(density(D_split[[30]][,2]),type="l",main=expression(paste("Cx45ko 4, 806",mu,"m",sep="")),col="red",xlab="",ylab="density",xlim=c(-0.05,0.4),xaxt="n",cex.main=0.85)
axis(1,at=c(0,0.2,0.4),labels=TRUE)

text(grconvertX(c(0.04, 0.51,0.04,0.51,0.04,0.51,0.77,0.51,0.77), from='ndc'),
   grconvertY(c(0.96, 0.96,0.65,0.65,0.31,0.31,0.31,0.17,0.17), from='ndc'),  c('A',  'B',  'C','D','E','F','G','H','I'), xpd=NA, cex=2.2, font=2)

text(grconvertX(c(0.66, 0.9), from='ndc'),
   grconvertY(c(0.02, 0.02), from='ndc'),  c('STTC', 'STTC'), xpd=NA, cex=0.8, font=1)

dev.off()






