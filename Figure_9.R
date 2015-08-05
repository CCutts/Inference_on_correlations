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

recordings= c("WT_1","WT_2","WT_3","WT_4","WT_5","cx45_1","cx45_2","cx45_3","cx45_4","cx3645_1","cx3645_2","cx3645_3","cx3645_5","cx3645_6")
load("Bshp11_MASTER.RData")

#extract the recordings were using (a couple are excluded as outliers)

r=subset(a,subset=(Animal %in% c("WT_1","WT_2","WT_3","WT_4","WT_5","cx45_1","cx45_2","cx45_3","cx45_4","cx3645_1","cx3645_2","cx3645_3","cx3645_5","cx3645_6")),select=c("distance","TC","Phenotype","Animal"))



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
i=2

WT=subset(r,subset=(Animal==recordings[i]),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))
distance_recording=as.numeric(names(WT_split))

wild_scale=rep(0,24)
wild_loc=rep(0,24)

for(j in 1:24){

if(length(WT_split[[j]][,2])>2){

	beta=sqrt(6*var(WT_split[[j]][,2])/(pi*pi))
	mu=mean(WT_split[[j]][,2])-beta*eulergamma
	a2=fitdistr(WT_split[[j]][,2],fun_gb,start=list(scale=beta,location=mu))
	wild_scale[j]=a2$estimate[1]
	wild_loc[j]=a2$estimate[2]

}		
}


#####################################################
#Cx45ko
#####################################################
i=8

WT=subset(r,subset=(Animal==recordings[i]),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))
distance_recording=as.numeric(names(WT_split))

cx45_scale=rep(0,24)
cx45_loc=rep(0,24)

for(j in 1:24){

if(length(WT_split[[j]][,2])>2){

	beta=sqrt(6*var(WT_split[[j]][,2])/(pi*pi))
	mu=mean(WT_split[[j]][,2])-beta*eulergamma
	a2=fitdistr(WT_split[[j]][,2],fun_gb,start=list(scale=beta,location=mu))
	cx45_scale[j]=a2$estimate[1]
	cx45_loc[j]=a2$estimate[2]

}		


}




#####################################################
#Cx3645ko
#####################################################
i=13

WT=subset(r,subset=(Animal==recordings[i]),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))
distance_recording=as.numeric(names(WT_split))

cx36_scale=rep(0,24)
cx36_loc=rep(0,24)

for(j in 1:24){
	
	if(length(WT_split[[j]][,2])>2){
		
		beta=sqrt(6*var(WT_split[[j]][,2])/(pi*pi))
		mu=mean(WT_split[[j]][,2])-beta*eulergamma
		a2=fitdistr(WT_split[[j]][,2],fun_gb,start=list(scale=beta,location=mu))
		cx36_scale[j]=a2$estimate[1]
		cx36_loc[j]=a2$estimate[2]
		
	}		
	
	
}



################################################
#get the fits

distances=as.numeric(names(WT_split))
a1_w=fit_f1(wild_loc,distances[1:24],c(0.1,0.1,0.001))
a1_45=fit_f1(cx45_loc,distances[1:24],c(0.1,0.1,0.001))
a1_36=fit_f1(cx36_loc,distances[1:24],c(0.1,0.1,0.001))




a3_w=fit_f3(wild_scale,distances[1:24],c(0.04,0.08,0.0001))
a3_45=fit_f3(cx45_scale,distances[1:24],c(0.04,0.08,0.0001))
a3_36=fit_f3(cx36_scale,distances[1:24],c(0.04,0.08,0.0001))


##############################################
#Plot

#############################################
pdf(file="Figure_9.pdf",width=16,height=16)
#recordings
recs=c(1,2,4,7,10,14,19,24)
layout(matrix(c(rep(1,3),rep(2,3),3,rep(4,3),rep(5,3),6,rep(7,3),rep(8,3),rep(1,3),rep(2,3),3,rep(4,3),rep(5,3),6,rep(7,3),rep(8,3),rep(9,3),rep(10,3),11,rep(12,3),rep(13,3),14,rep(15,3),rep(16,3),rep(9,3),rep(10,3),11,rep(12,3),rep(13,3),14,rep(15,3),rep(16,3),
rep(17,3),rep(18,3),19,rep(20,3),rep(21,3),22,rep(23,3),rep(24,3),rep(17,3),rep(18,3),19,rep(20,3),rep(21,3),22,rep(23,3),rep(24,3),rep(25,3),rep(26,3),27,rep(28,3),rep(29,3),30,rep(31,3),rep(32,3),rep(25,3),rep(26,3),27,rep(28,3),rep(29,3),30,rep(31,3),rep(32,3),
rep(33,3),rep(34,3),35,rep(36,3),rep(37,3),38,rep(39,3),rep(40,3),rep(41,3),rep(42,3),43,rep(44,3),rep(45,3),46,rep(47,3),rep(48,3),rep(41,3),rep(42,3),43,rep(44,3),rep(45,3),46,rep(47,3),rep(48,3)), 11, 20, byrow = TRUE))

par(oma=c(1,1,4,1))
par(mar=c(4,4,2,1))



# plot wild type 0 
#
#
i=2
j=1
WT=subset(r,subset=(Animal==recordings[i]),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))


	
	
x3=rGumbel(10000,mu=wild_loc[j],sigma=wild_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)

plot(d,main=expression(paste("Electrode separation: 0",mu,"m",sep="")),xlim=c(min(b),1.1),ylim=c,xlab="",ylab="",yaxt="n",cex.main=1.55,lwd=1.5,cex.axis=1.5)
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3)
	points(g$x,g$y,type="l",lty=2,lwd=1.5) #gumbel
	legend("topright",c("u=0.440","v=0.160"),bty="n",cex=1.5)
		

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

plot(d,main=expression(paste("100",mu,"m",sep="")),xlim=c(min(b),1),ylim=c,xlab="",ylab="",yaxt="n",cex.main=1.6,lwd=1.5,cex.axis=1.5)			
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3)			
points(g$x,g$y,type="l",lty=2,lwd=1.5) #gumbel
legend("topright",c("u=0.236","v=0.138"),bty="n",cex=1.5)



plot(1, type="n", axes=F, xlab="", ylab="")

######################################################		
# plot Cx45ko 0 
#
#
i=8
j=1
WT=subset(r,subset=(Animal==recordings[i]),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))


	
	
x3=rGumbel(10000,mu=cx45_loc[j],sigma=cx45_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)
par(mar=c(4,4,2,1))
plot(d,main=expression(paste("0",mu,"m",sep="")),xlim=c(min(b),1),ylim=c,xlab="",ylab="",yaxt="n",cex.main=1.6,lwd=1.5,col="red",cex.axis=1.5)
	rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3,col="red")
	points(g$x,g$y,type="l",col="red",lty=2,lwd=1.5) #gumbel
	legend("topright",c("u=0.481","v=0.118"),bty="n",cex=1.5)

# Cx45 100 
#
#	
#
j=2
	
x3=rGumbel(10000,mu=cx45_loc[j],sigma=cx45_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)

plot(d,main=expression(paste("100",mu,"m",sep="")),xlim=c(min(b),1),ylim=c,xlab="",ylab="",yaxt="n",cex.main=1.6,lwd=1.5,col="red",cex.axis=1.5)			
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3,col="red")			
points(g$x,g$y,type="l",col="red",lty=2,lwd=1.5) #gumbel
legend("topright",c("u=0.372","v=0.145"),bty="n",cex=1.5)
			


plot(1, type="n", axes=F, xlab="", ylab="")

######################################################		
# plot Cx3645ko 0 
#
#
i=13
j=1
WT=subset(r,subset=(Animal==recordings[i]),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))




x3=rGumbel(10000,mu=cx36_loc[j],sigma=cx36_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)

plot(d,main=expression(paste("0",mu,"m",sep="")),xlim=c(min(b),1),ylim=c,xlab="",ylab="",yaxt="n",cex.main=1.6,lwd=1.5,cex.axis=1.5,col="blue")
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3,col="blue")
points(g$x,g$y,type="l",col="blue",lty=2,lwd=1.5) #gumbel
legend("topright",c("u=0.168","v=0.169"),bty="n",cex=1.5)

# Cx3645 100 
#
#	
#
j=2

x3=rGumbel(10000,mu=cx36_loc[j],sigma=cx36_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)

plot(d,main=expression(paste("100",mu,"m",sep="")),xlim=c(min(b),1),ylim=c,xlab="",ylab="",yaxt="n",cex.main=1.6,cex.axis=1.5,lwd=1.5,col="blue")			
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3,col="blue")			
points(g$x,g$y,type="l",col="blue",lty=2,lwd=1.5) #gumbel
legend("topright",c("u=0.100","v=0.118"),bty="n",cex=1.5)


###########################################################################################
# plot wild type 200 
#
#
i=2
j=4

WT=subset(r,subset=(Animal==recordings[i]),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))


	
	
x3=rGumbel(10000,mu=wild_loc[j],sigma=wild_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)




plot(d,main=expression(paste("200",mu,"m",sep="")),xlim=b,ylim=c,xlab="",ylab="",yaxt="n",cex.main=1.6,lwd=1.5,cex.axis=1.5)
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3)
points(g$x,g$y,type="l",lty=2,lwd=1.5) #gumbel
legend("topright",c("u=0.175","v=0.121"),bty="n",cex=1.5)

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



plot(d,main=expression(paste("300",mu,"m",sep="")),xlim=b,ylim=c,xlab="",ylab="",yaxt="n",cex.main=1.6,lwd=1.5,cex.axis=1.5)
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3)
points(g$x,g$y,type="l",lty=2,lwd=1.5) #gumbel
legend("topright",c("u=0.113","v=0.092"),bty="n",cex=1.5)



plot(1, type="n", axes=F, xlab="", ylab="")

##################################################
# plot Cx45 200 
#
#
i=8
j=4

WT=subset(r,subset=(Animal==recordings[i]),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))


	
	
x3=rGumbel(10000,mu=cx45_loc[j],sigma=cx45_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)




plot(d,main=expression(paste("200",mu,"m",sep="")),xlim=b,ylim=c,xlab="",ylab="",yaxt="n",cex.main=1.6,lwd=1.5,cex.axis=1.5,col="red")
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3,col="red")
points(g$x,g$y,type="l",col="red",lty=2,lwd=1.5) #gumbel
legend("topright",c("u=0.187","v=0.095"),bty="n",cex=1.5)

#plot Cx45 300
#
j=7			

WT=subset(r,subset=(Animal==recordings[i]),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))

	
	
x3=rGumbel(10000,mu=wild_loc[j],sigma=wild_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)

plot(d,main=expression(paste("300",mu,"m",sep="")),xlim=b,ylim=c,xlab="",ylab="",yaxt="n",cex.main=1.6,lwd=1.5,cex.axis=1.5,col="red")
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3,col="red")
points(g$x,g$y,type="l",col="red",lty=2,lwd=1.5) #gumbel
legend("topright",c("u=0.090","v=0.058"),bty="n",cex=1.5)
			


plot(1, type="n", axes=F, xlab="", ylab="")

##################################################
# plot Cx3645 200 
#
#
i=13
j=4

WT=subset(r,subset=(Animal==recordings[i]),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))




x3=rGumbel(10000,mu=cx36_loc[j],sigma=cx36_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)




plot(d,main=expression(paste("200",mu,"m",sep="")),xlim=b,ylim=c,xlab="",ylab="",yaxt="n",cex.main=1.6,lwd=1.5,cex.axis=1.5,col="blue")
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3,col="blue")
points(g$x,g$y,type="l",col="blue",lty=2,lwd=1.5) #gumbel
legend("topright",c("u=0.096","v=0.123"),bty="n",cex=1.5)

#plot Cx45 300
#
j=7			

WT=subset(r,subset=(Animal==recordings[i]),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))


x3=rGumbel(10000,mu=cx36_loc[j],sigma=cx36_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)

plot(d,main=expression(paste("300",mu,"m",sep="")),xlim=b,ylim=c,xlab="",ylab="",yaxt="n",cex.main=1.6,lwd=1.5,cex.axis=1.5,col="blue")
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3,col="blue")
points(g$x,g$y,type="l",col="blue",lty=2,lwd=1.5) #gumbel
legend("topright",c("u=0.069","v=0.100"),bty="n",cex=1.5)



#plot wild type 400 
#
#
#

i=2
j=10

WT=subset(r,subset=(Animal==recordings[i]),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))


	
x3=rGumbel(10000,mu=wild_loc[j],sigma=wild_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)


plot(d,main=expression(paste("400",mu,"m",sep="")),xlim=b,ylim=c,xlab="",ylab="",yaxt="n",cex.main=1.6,lwd=1.5,cex.axis=1.5)
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3)			
points(g$x,g$y,type="l",lty=2,lwd=1.5) #gumbel
legend("topright",c("u=0.062","v=0.047"),bty="n",cex=1.5)

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
			
plot(d,main=expression(paste("500",mu,"m",sep="")),xlim=b,ylim=c,xlab="",ylab="",yaxt="n",cex.main=1.6,lwd=1.5,cex.axis=1.5)	
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3)			
points(g$x,g$y,type="l",lty=2,lwd=1.5) #gumbel
legend("topright",c("u=0.041","v=0.032"),bty="n",cex=1.5)
	

plot(1, type="n", axes=F, xlab="", ylab="")

###########################################
#plot cx45  400 
#
#
#

i=8
j=10

WT=subset(r,subset=(Animal==recordings[i]),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))


	
x3=rGumbel(10000,mu=cx45_loc[j],sigma=cx45_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)


plot(d,main=expression(paste("400",mu,"m",sep="")),xlim=b,ylim=c,xlab="",ylab="",yaxt="n",cex.main=1.6,lwd=1.5,cex.axis=1.5,col="red")
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3,col="red")			
points(g$x,g$y,type="l",col="red",lty=2,lwd=1.5) #gumbel
legend("topright",c("u=0.043","v=0.043"),bty="n",cex=1.5)

#plot cx45 500
#
#			
j=14



	
x3=rGumbel(10000,mu=cx45_loc[j],sigma=cx45_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)
			
plot(d,main=expression(paste("500",mu,"m",sep="")),xlim=b,ylim=c,xlab="",ylab="",yaxt="n",cex.main=1.6,cex.axis=1.5,lwd=1.5,col="red")	
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3,col="red")			
points(g$x,g$y,type="l",col="red",lty=2,lwd=1.5) #gumbel
legend("topright",c("u=0.024","v=0.040"),bty="n",cex=1.5)

plot(1, type="n", axes=F, xlab="", ylab="")
		
###########################################
#plot cx3645  400 
#
#
#

i=13
j=10

WT=subset(r,subset=(Animal==recordings[i]),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))



x3=rGumbel(10000,mu=cx36_loc[j],sigma=cx36_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)


plot(d,main=expression(paste("400",mu,"m",sep="")),xlim=b,ylim=c,xlab="",ylab="",yaxt="n",cex.main=1.6,lwd=1.5,col="blue",cex.axis=1.5)
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3,col="blue")			
points(g$x,g$y,type="l",col="blue",lty=2,lwd=1.5) #gumbel
legend("topright",c("u=0.065","v=0.102"),bty="n",cex=1.5)

#plot cx3645 500
#
#			
j=14




x3=rGumbel(10000,mu=cx36_loc[j],sigma=cx36_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)

plot(d,main=expression(paste("500",mu,"m",sep="")),xlim=b,ylim=c,xlab="",ylab="",yaxt="n",cex.main=1.6,lwd=1.5,col="blue",cex.axis=1.5)	
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3,col="blue")			
points(g$x,g$y,type="l",col="blue",lty=2,lwd=1.5) #gumbel
legend("topright",c("u=0.063","v=0.101"),bty="n",cex=1.5)


#plot wild type 600
#
#
i=2
j=19


WT=subset(r,subset=(Animal==recordings[i]),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))


	
x3=rGumbel(10000,mu=wild_loc[j],sigma=wild_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)

plot(d,main=expression(paste("600",mu,"m",sep="")),xlim=b,ylim=c,ylab="",xlab="",yaxt="n",cex.main=1.6,lwd=1.5,cex.axis=1.5)
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3)
points(g$x,g$y,type="l",lty=2,lwd=1.5) #gumbel
legend("topright",c("u=0.038","v=0.021"),bty="n",cex=1.5)
		

# plot wild type 700
#
#
j=24

	
x3=rGumbel(10000,mu=wild_loc[j],sigma=wild_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)


plot(d,main=expression(paste("700",mu,"m",sep="")),xlim=b,ylim=c,xlab="STTC",ylab="Density",yaxt="n",cex.main=1.6,lwd=1.5,cex.axis=1.6)				
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3)
points(g$x,g$y,type="l",lty=2,lwd=1.5) #gumbel
#legend("topright",lwd=c(3,1.5,1.5,1.5,1.5),col=c("black","red",259,"brown","blue"),legend=c("Data","Normal","Skew-normal","EMG","Gumbel"))
legend("topright",c("u=0.076","v=0.067"),bty="n",cex=1.5)
		
plot(1, type="n", axes=F, xlab="", ylab="")
#plot cx45  600
#
#
i=8
j=19


WT=subset(r,subset=(Animal==recordings[i]),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))


	
x3=rGumbel(10000,mu=cx45_loc[j],sigma=cx45_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)

plot(d,main=expression(paste("600",mu,"m",sep="")),xlim=b,ylim=c,ylab="",xlab="",yaxt="n",cex.main=1.6,lwd=1.5,cex.axis=1.5,col="red")
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3,col="red")
points(g$x,g$y,type="l",lty=2,col="red",lwd=1.5) #gumbel
legend("topright",c("u=0.009","v=0.038"),bty="n",cex=1.5)
			

# plot cx45 700
#
#
j=24

	
x3=rGumbel(10000,mu=cx45_loc[j],sigma=cx45_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)

plot(d,main=expression(paste("700",mu,"m",sep="")),xlim=b,ylim=c,xlab="",ylab="",yaxt="n",cex.main=1.6,lwd=1.5,cex.axis=1.5,col="red")				
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3,col="red")
points(g$x,g$y,type="l",col="red",lty=2,lwd=1.5) #gumbel
legend("right",c("u=0.023","v=0.023"),bty="n",cex=1.5)
legend("topright",lty=c(1,2),col=c("red"),legend=c("Data","Gumbel"),cex=1.4)
plot(1, type="n", axes=F, xlab="", ylab="")
#plot cx3645  600
#
#
i=13
j=19


WT=subset(r,subset=(Animal==recordings[i]),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))



x3=rGumbel(10000,mu=cx36_loc[j],sigma=cx36_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)

plot(d,main=expression(paste("600",mu,"m",sep="")),xlim=b,ylim=c,ylab="",xlab="",yaxt="n",cex.main=1.6,lwd=1.5,cex.axis=1.5,col="blue")
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3,col="blue")
points(g$x,g$y,type="l",lty=2,col="blue",lwd=1.5) #gumbel
legend("topright",c("u=0.036","v=0.100"),bty="n",cex=1.5)


# plot beta 2 700
#
#
j=24


x3=rGumbel(10000,mu=cx36_loc[j],sigma=cx36_scale[j])
d=density(WT_split[[j]][,2])			
g=density(x3)
b=range(x3,d$x,na.rm=TRUE)
c=range(d$y,g$y)

plot(d,main=expression(paste("700",mu,"m",sep="")),xlim=b,ylim=c,xlab="",ylab="",yaxt="n",cex.main=1.6,lwd=1.5,cex.axis=1.5,col="blue")				
rug(WT_split[[j]][,2], ticksize = 0.1, side = 1,lwd=0.3,col="blue")
points(g$x,g$y,type="l",col="blue",lty=2,lwd=1.5) #gumbel

legend("topright",c("u=0.089","v=0.087"),bty="n",cex=1.5)





cat("ok to here")



plot(1, type="n", axes=F, xlab="", ylab="")
plot(1, type="n", axes=F, xlab="", ylab="")
plot(1, type="n", axes=F, xlab="", ylab="")
plot(1, type="n", axes=F, xlab="", ylab="")
plot(1, type="n", axes=F, xlab="", ylab="")
plot(1, type="n", axes=F, xlab="", ylab="")
plot(1, type="n", axes=F, xlab="", ylab="")
plot(1, type="n", axes=F, xlab="", ylab="")

dist=c(0,100,200,300,400,500,600,700)
u=c(0.355,0.098,0.017,-0.010,-0.005,-0.009,-0.008,-0.007)
v=c(0.158,0.102,0.036,0.024,0.017,0.018,0.016,0.021)

distances=as.numeric(names(WT_split))
par(las=1)

plot(distances[1:24],wild_loc,xlab=expression(paste("Intercell distance (",mu,"m)",sep="")),ylim=c(0,0.52),ylab="u",main="u",pch=19,cex.axis=1.5,cex.lab=1.5,lwd=1.5,cex.main=2)
points(distances[1:24],f1(distances[1:24],a1_w$coefficients[1],a1_w$coefficients[2],a1_w$coefficients[3]),type="l")
plot(distances[1:24],wild_scale,xlab="",ylab="v",main="v",ylim=c(0.02,0.18),pch=19,cex.axis=1.5,cex.lab=1.5,lwd=1.5,cex.main=2)
points(distances[1:24],f3(distances[1:24],a3_w$coefficients[1],a3_w$coefficients[2],a3_w$coefficients[3]),type="l")

plot(1, type="n", axes=F, xlab="", ylab="")

plot(distances[1:24],cx45_loc,xlab="",ylim=c(0,0.52),ylab="u",main="u",col="red",pch=19,cex.axis=1.5,cex.lab=1.5,lwd=1.5,cex.main=2)
points(distances[1:24],f1(distances[1:24],a1_45$coefficients[1],a1_45$coefficients[2],a1_45$coefficients[3]),type="l",col="red")
plot(distances[1:24],cx45_scale,xlab="",ylab="v",main="v",ylim=c(0.02,0.18),col="red",pch=19,cex.axis=1.5,cex.lab=1.5,lwd=1.5,cex.main=2)
points(distances[1:24],f3(distances[1:24],a3_45$coefficients[1],a3_45$coefficients[2],a3_45$coefficients[3]),type="l",col="red")

plot(1, type="n", axes=F, xlab="", ylab="")

plot(distances[1:24],cx36_loc,xlab="",ylim=c(0,0.52),ylab="u",main="u",col="blue",pch=19,cex.axis=1.5,cex.lab=1.5,lwd=1.5,cex.main=2)
points(distances[1:24],f1(distances[1:24],a1_36$coefficients[1],a1_36$coefficients[2],a1_36$coefficients[3]),type="l",col="blue")
plot(distances[1:24],cx36_scale,xlab="",ylab="v",main="v",ylim=c(0.02,0.18),col="blue",pch=19,cex.axis=1.5,cex.lab=1.5,lwd=1.5,cex.main=2)
points(distances[1:24],f3(distances[1:24],a3_36$coefficients[1],a3_36$coefficients[2],a3_36$coefficients[3]),type="l",col="blue")

legend("topright",pch=c(19,NA),lty=c(NA,1),legend=c("MLE","L-S"),col="blue",cex=1.75)

text(grconvertX(c(0.02, 0.02,0.36,0.36,0.7,0.7), from='ndc'),
   grconvertY(c(0.95, 0.18,0.95,0.18,0.95,0.18), from='ndc'),  c('A','D','B','E','C','F'), xpd=NA, cex=3.5, font=2)

text(grconvertX(c(0.15, 0.51 , 0.85), from='ndc'),
   grconvertY(c(0.99), from='ndc'),  c('Wild type','Cx45ko','Cx36/45dko'), xpd=NA, cex=2.4, font=2)
	
dev.off()


