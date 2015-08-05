rec=28
ndist=32

library(MASS)
library(fGarch)
library(geoR)
library(QRM)
library(emg)
library(pracma)
library(e1071)
library(graphics)
set.seed(25)
library(nlme)
eulergamma=0.57721566490153
#start with distributions over the real line

recordings= c("WTP4_1","WTP4_2","WTP4_3","WTP4_4","WTP4_5","WTP4_6","WTP4_7","WTP4_8","WTP4_9","WTP4_10","WTP4_11","WTP4_12","B2TGP4_1","B2TGP4_3",
"B2TGP4_4","B2TGP4_5","B2TGP4_6","B2TGP4_7","B2TGP4_8","B2TGP4_9","B2TGP4_10","B2TGP4_11","B2TGP4_12","B2TGP4_13","B2TGP4_14","B2TGP4_15","B2TGP4_16","B2TGP4_17")

load("Xu_MASTER.RData")
#extract the recordings we're using (a couple are excluded as outliers)
r=subset(a,subset=(Animal %in% c("WTP4_1","WTP4_2","WTP4_3","WTP4_4","WTP4_5","WTP4_6","WTP4_7","WTP4_8","WTP4_9","WTP4_10","WTP4_11","WTP4_12","B2TGP4_1","B2TGP4_3",
"B2TGP4_4","B2TGP4_5","B2TGP4_6","B2TGP4_7","B2TGP4_8","B2TGP4_9","B2TGP4_10","B2TGP4_11","B2TGP4_12","B2TGP4_13","B2TGP4_14","B2TGP4_15","B2TGP4_16","B2TGP4_17")),select=c("distance","TC","Phenotype","Animal"))



gb_loc=mat.or.vec(rec,ndist)
gb_scale=mat.or.vec(rec,ndist)


#gumbel pdf 
fun_gb=function(x,scale,location){
z=(x-location)/scale
a=(1/scale)*exp(-(z+exp(-z)))
return(a)
}

wild_scale=rep(0,29)
wild_loc=rep(0,29)

b2_scale=rep(0,29)
b2_loc=rep(0,29)


f1=function(x,A,B,C){
y=A+B*exp(-C*x)}

f3=function(x,A,B,C){
y=A+B*exp(-C*x^2)}


fit_f1=function(WT_locs,dist,begin){
p=gnls(WT_locs~A+B*exp(-C*dist),start=c(A=begin[1],B=begin[2],C=begin[3]),na.action=na.omit)
return(p)
}


fit_f3=function(WT_locs,dist,begin){
p=gnls(WT_locs~A+B*exp(-C*(dist^2)),start=c(A=begin[1],B=begin[2],C=begin[3]),na.action=na.omit)
return(p)
}



#####################################################
# wild type
#####################################################
i=1

WT=subset(r,subset=(Animal==recordings[i]),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))

for(j in 1:29){

if(length(WT_split[[j]][,2])>10){

	beta=sqrt(6*var(WT_split[[j]][,2])/(pi*pi))
	mu=mean(WT_split[[j]][,2])-beta*eulergamma
	a2=fitdistr(WT_split[[j]][,2],fun_gb,start=list(scale=beta,location=mu))
	wild_scale[j]=a2$estimate[1]
	wild_loc[j]=a2$estimate[2]

}		
}


#####################################################
#beta 2 
#####################################################
i=16

WT=subset(r,subset=(Animal==recordings[i]),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))

for(j in 1:29){

if(length(WT_split[[j]][,2])>10){

	beta=sqrt(6*var(WT_split[[j]][,2])/(pi*pi))
	mu=mean(WT_split[[j]][,2])-beta*eulergamma
	a2=fitdistr(WT_split[[j]][,2],fun_gb,start=list(scale=beta,location=mu))
	b2_scale[j]=a2$estimate[1]
	b2_loc[j]=a2$estimate[2]

}		


}

################################################
#get the fits

distances=as.numeric(names(WT_split))
a1_w=fit_f1(wild_loc,distances[1:29],c(0.1,0.1,0.001))

a1_b=fit_f1(b2_loc,distances[1:29],c(0.1,0.1,0.001))



a3_w=fit_f3(wild_scale,distances[1:29],c(0.04,0.08,0.0001))

a3_b=fit_f3(b2_scale,distances[1:29],c(0.04,0.08,0.0001))



##############################################
#Plot

#############################################
pdf(file="Figure_6.pdf",width=8,height=10)
#recordings
recs=c(1,2,4,7,10,14,19,24)
par(mfrow=c(5,4))
par(oma=c(1,1,4,1))
par(mar=c(4,4,2,1))



# plot wild type 0 
#
#
i=1
j=1
WT=subset(r,subset=(Animal==recordings[i]),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))


	
	
x3=rGumbel(10000,mu=wild_loc[j],sigma=wild_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)

plot(d,main=expression(paste("Electrode separation: 0",mu,"m",sep="")),xlim=c(min(b),1.1),ylim=c,xlab="",ylab="",yaxt="n",cex.main=1,lwd=1)
	rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3)
	points(g$x,g$y,type="l",lty=2,lwd=1.5) #gumbel
	legend("topright",c("u=0.354","v=0.175"),bty="n")
		

# plot wild type 100 
#
#	
#
j=2
	
x3=rGumbel(10000,mu=wild_loc[j],sigma=wild_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)

plot(d,main=expression(paste("100",mu,"m",sep="")),xlim=c(min(b),1),ylim=c,xlab="",ylab="",yaxt="n",cex.main=1,lwd=1)			
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3)			
points(g$x,g$y,type="l",lty=2,lwd=1.5) #gumbel
legend("topright",c("u=0.214","v=0.094"),bty="n")
	

		
# plot beta 2 0 
#
#
i=16
j=1
WT=subset(r,subset=(Animal==recordings[i]),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))


	
	
x3=rGumbel(10000,mu=wild_loc[j],sigma=wild_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)

plot(d,main=expression(paste("0",mu,"m",sep="")),xlim=c(min(b),1),ylim=c,xlab="",ylab="",yaxt="n",cex.main=1,lwd=1,col="blue")
	rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3,col="blue")
	points(g$x,g$y,type="l",col="blue",lty=2,lwd=1.5) #gumbel
	legend("topright",c("u=0.355","v=0.158"),bty="n")
		

# plot beta 2 100 
#
#	
#
j=2
	
x3=rGumbel(10000,mu=wild_loc[j],sigma=wild_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)

plot(d,main=expression(paste("100",mu,"m",sep="")),xlim=c(min(b),1),ylim=c,xlab="",ylab="",yaxt="n",cex.main=1,lwd=1,col="blue")			
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3,col="blue")			
points(g$x,g$y,type="l",col="blue",lty=2,lwd=1.5) #gumbel
legend("topright",c("u=0.098","v=0.102"),bty="n")
			

# plot wild type 200 
#
#
i=1
j=4

WT=subset(r,subset=(Animal==recordings[i]),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))
	
	
x3=rGumbel(10000,mu=wild_loc[j],sigma=wild_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)




plot(d,main=expression(paste("200",mu,"m",sep="")),xlim=b,ylim=c,xlab="",ylab="",yaxt="n",cex.main=1,lwd=1)
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3)
points(g$x,g$y,type="l",lty=2,lwd=1.5) #gumbel
legend("topright",c("u=0.051","v=0.052"),bty="n")

#plot wild type 300
#
j=7			

WT=subset(r,subset=(Animal==recordings[i]),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))

	
x3=rGumbel(10000,mu=wild_loc[j],sigma=wild_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)



plot(d,main=expression(paste("300",mu,"m",sep="")),xlim=b,ylim=c,xlab="",ylab="",yaxt="n",cex.main=1,lwd=1)
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3)
points(g$x,g$y,type="l",lty=2,lwd=1.5) #gumbel
legend("topright",c("u=0.009","v=0.034"),bty="n")

# plot beta 2 200 
#
#
i=16
j=4

WT=subset(r,subset=(Animal==recordings[i]),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))


	
	
x3=rGumbel(10000,mu=wild_loc[j],sigma=wild_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)




plot(d,main=expression(paste("200",mu,"m",sep="")),xlim=b,ylim=c,xlab="",ylab="",yaxt="n",cex.main=1,lwd=1,col="blue")
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3,col="blue")
points(g$x,g$y,type="l",col="blue",lty=2,lwd=1.5) #gumbel
legend("topright",c("u=0.017","v=0.036"),bty="n")

#plot beta 2 300
#
j=7			

WT=subset(r,subset=(Animal==recordings[i]),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))

	
	
x3=rGumbel(10000,mu=wild_loc[j],sigma=wild_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)

plot(d,main=expression(paste("300",mu,"m",sep="")),xlim=b,ylim=c,xlab="",ylab="",yaxt="n",cex.main=1,lwd=1,col="blue")
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3,col="blue")
points(g$x,g$y,type="l",col="blue",lty=2,lwd=1.5) #gumbel
legend("topright",c("u=-0.010","v=  0.024"),bty="n")
			



#plot wild type 400 
#
#
#

i=1
j=10

WT=subset(r,subset=(Animal==recordings[i]),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))


	
x3=rGumbel(10000,mu=wild_loc[j],sigma=wild_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)


plot(d,main=expression(paste("400",mu,"m",sep="")),xlim=b,ylim=c,xlab="",ylab="",yaxt="n",cex.main=1,lwd=1)
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3)			
points(g$x,g$y,type="l",lty=2,lwd=1.5) #gumbel
legend("topright",c("u=-0.010","v=  0.026"),bty="n")

#plot wild type 500
#
#			
j=14


WT=subset(r,subset=(Animal==recordings[i]),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))


	
x3=rGumbel(10000,mu=wild_loc[j],sigma=wild_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)
			
plot(d,main=expression(paste("500",mu,"m",sep="")),xlim=b,ylim=c,xlab="",ylab="",yaxt="n",cex.main=1,lwd=1)	
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3)			
points(g$x,g$y,type="l",lty=2,lwd=1.5) #gumbel
legend("topright",c("u=-0.012","v=  0.022"),bty="n")
			
#plot beta 2 400 
#
#
#

i=16
j=10

WT=subset(r,subset=(Animal==recordings[i]),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))


	
x3=rGumbel(10000,mu=wild_loc[j],sigma=wild_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)


plot(d,main=expression(paste("400",mu,"m",sep="")),xlim=b,ylim=c,xlab="",ylab="",yaxt="n",cex.main=1,lwd=1,col="blue")
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3,col="blue")			
points(g$x,g$y,type="l",col="blue",lty=2,lwd=1.5) #gumbel
legend("topright",c("u=-0.005","v=  0.017"),bty="n")

#plot beta 2 500
#
#			
j=14



	
x3=rGumbel(10000,mu=wild_loc[j],sigma=wild_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)
			
plot(d,main=expression(paste("500",mu,"m",sep="")),xlim=b,ylim=c,xlab="",ylab="",yaxt="n",cex.main=1,lwd=1,col="blue")	
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3,col="blue")			
points(g$x,g$y,type="l",col="blue",lty=2,lwd=1.5) #gumbel
legend("topright",c("u=-0.009","v=  0.018"),bty="n")
		


#plot wild type 600
#
#
i=1
j=19


WT=subset(r,subset=(Animal==recordings[i]),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))


	
x3=rGumbel(10000,mu=wild_loc[j],sigma=wild_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)

plot(d,main=expression(paste("600",mu,"m",sep="")),xlim=b,ylim=c,ylab="",xlab="",yaxt="n",cex.main=1,lwd=1)
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3)
points(g$x,g$y,type="l",lty=2,lwd=1.5) #gumbel
legend("topright",c("u=-0.018","v=  0.024"),bty="n")
			

# plot wild type 700
#
#
j=24

	
x3=rGumbel(10000,mu=wild_loc[j],sigma=wild_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)


plot(d,main=expression(paste("700",mu,"m",sep="")),xlim=b,ylim=c,xlab="STTC",ylab="Density",yaxt="n",cex.main=1,lwd=1)				
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3)
points(g$x,g$y,type="l",lty=2,lwd=1.5) #gumbel
#legend("topright",lwd=c(3,1.5,1.5,1.5,1.5),col=c("black","red",259,"brown","blue"),legend=c("Data","Normal","Skew-normal","EMG","Gumbel"))
legend("topright",c("u=-0.022","v=  0.025"),bty="n")
		

#plot beta 2  600
#
#
i=16
j=19


WT=subset(r,subset=(Animal==recordings[i]),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))


	
x3=rGumbel(10000,mu=wild_loc[j],sigma=wild_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)

plot(d,main=expression(paste("600",mu,"m",sep="")),xlim=b,ylim=c,ylab="",xlab="",yaxt="n",cex.main=1,lwd=1,col="blue")
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3,col="blue")
points(g$x,g$y,type="l",lty=2,col="blue",lwd=1.5) #gumbel
legend("topright",c("u=-0.008","v=  0.016"),bty="n")
			

# plot beta 2 700
#
#
j=24

	
x3=rGumbel(10000,mu=wild_loc[j],sigma=wild_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)

plot(d,main=expression(paste("700",mu,"m",sep="")),xlim=b,ylim=c,xlab="",ylab="",yaxt="n",cex.main=1,lwd=1,col="blue")				
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3,col="blue")
points(g$x,g$y,type="l",col="blue",lty=2,lwd=1.5) #gumbel
legend("topright",lty=c(1,2),col=c("blue"),legend=c("Data","Gumbel"))
legend("right",c("u=-0.007","v=  0.021"),bty="n")
		



par(mar=c(4,4,3,1))
dist=c(0,100,200,300,400,500,600,700)
u=c(0.355,0.098,0.017,-0.010,-0.005,-0.009,-0.008,-0.007)
v=c(0.158,0.102,0.036,0.024,0.017,0.018,0.016,0.021)

distances=as.numeric(names(WT_split))


plot(distances[1:29],wild_loc,xlab=expression(paste("Intercell distance (",mu,"m)",sep="")),ylim=c(-0.03,0.36),ylab="u",main="u",pch=19)
points(distances[1:29],f1(distances[1:29],a1_w$coefficients[1],a1_w$coefficients[2],a1_w$coefficients[3]),type="l")
plot(distances[1:29],wild_scale,xlab="",ylab="v",main="v",ylim=c(0.0,0.18),pch=19)
points(distances[1:29],f3(distances[1:29],a3_w$coefficients[1],a3_w$coefficients[2],a3_w$coefficients[3]),type="l")

plot(distances[1:29],b2_loc,xlab="",ylim=c(-0.03,0.36),ylab="u",main="u",col="blue",pch=19)
points(distances[1:29],f1(distances[1:29],a1_b$coefficients[1],a1_b$coefficients[2],a1_b$coefficients[3]),type="l",col="blue")
plot(distances[1:29],b2_scale,xlab="",ylab="v",main="v",ylim=c(0,0.18),col="blue",pch=19)
points(distances[1:29],f3(distances[1:29],a3_b$coefficients[1],a3_b$coefficients[2],a3_b$coefficients[3]),type="l",col="blue")
legend("topright",pch=c(19,NA),lty=c(NA,1),legend=c("MLE","L-S"),col="blue")

text(grconvertX(c(0.04, 0.04,0.54,0.54), from='ndc'),
   grconvertY(c(0.92, 0.17,0.92,0.17), from='ndc'),  c('A','C','B','D'), xpd=NA, cex=2, font=2)

text(grconvertX(c(0.29, 0.76), from='ndc'),
   grconvertY(c(0.98, 0.98), from='ndc'),  c('Wild type',expression(paste(beta,"2(TG)"))), xpd=NA, cex=2, font=2)
	
dev.off()


