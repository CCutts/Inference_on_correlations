rm(list = ls())
##########################
#Load all the data
library(QRM)



library(Scale)
library(scales)
library(fields)


#get the data 

load("Bshp11_MASTER.RData")
#extract the recordings we are using (a couple are excluded as outliers)
r=subset(a,subset=(Animal %in% c("WT_1","WT_2","WT_3","WT_4","WT_5","cx45_1","cx45_2","cx45_3","cx45_4","cx3645_1","cx3645_2","cx3645_3","cx3645_5","cx3645_6")),select=c("distance","TC","Phenotype","Animal"))

e=which(r[,4]=="WT_1")
r[e,4]=1
e=which(r[,4]=="WT_2")
r[e,4]=2
e=which(r[,4]=="WT_3")
r[e,4]=3
e=which(r[,4]=="WT_4")
r[e,4]=4
e=which(r[,4]=="WT_5")
r[e,4]=5
e=which(r[,4]=="cx45_1")
r[e,4]=6
e=which(r[,4]=="cx45_2")
r[e,4]=7
e=which(r[,4]=="cx45_3")
r[e,4]=8
e=which(r[,4]=="cx45_4")
r[e,4]=9
e=which(r[,4]=="cx3645_1")
r[e,4]=10
e=which(r[,4]=="cx3645_2")
r[e,4]=11
e=which(r[,4]=="cx3645_3")
r[e,4]=12

e=which(r[,4]=="cx3645_5")
r[e,4]=13
e=which(r[,4]=="cx3645_6")
r[e,4]=14

n=dim(r)[1]
rec=14
nWT=5
nKO=4
t=3
y=r[,2]
x=r[,1]
z=x^2


Rec=as.numeric(r[,4])



#Get our vector of distances
WT=subset(r,subset=(Phenotype=="WT"),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))
distances=as.numeric(names(WT_split))


library(rstan)

################################################################################################################################
#Load the model_fits

#Full
fit=read_stan_csv(c("Bshp_F_gb_v4_1.csv","Bshp_F_gb_v4_2.csv","Bshp_F_gb_v4_3.csv"))

Model_F=extract(fit,permuted=T)




#######################################################################################################################################
#We have 11 different models/parameter sets and we estimate y_sim for each of them 


y_sim=mat.or.vec(4,length(y))

#Some need to be done accounting for recording/phenotype
#n.b can be sped up, but don't see point right now. 
for(i in 1:length(y)){
j=Rec[i]

#FULL model
alpha=sample(Model_F$alpha[,j],1,replace=TRUE) 
beta=sample(Model_F$beta[,j],1,replace=TRUE) 
c=sample(Model_F$c[,j],1,replace=TRUE) 
A=sample(Model_F$A[,j],1,replace=TRUE)
B=sample(Model_F$B[,j],1, replace=TRUE)
C=sample(Model_F$C[,j],1 , replace=TRUE)
mean=alpha+beta*exp(-c*x[i])
sd=A+B*exp(-C*z[i])
if(sd<0){
sd=0
}
y_sim[1,i]=rGumbel(1, mean, sd)
}

###########################################################################################################################################
r_sim=r
r_sim[,2]=y_sim[1,]

WT=subset(r,subset=(Phenotype=="WT"),select=c("distance","TC"))
WT_split=split(WT,as.factor(WT$distance))
distance_WT=as.numeric(names(WT_split))

B2=subset(r,subset=(Phenotype=="cx45KO"),select=c("distance","TC"))
B2_split=split(B2,as.factor(B2$distance))
distance_B2=as.numeric(names(B2_split))

B3=subset(r,subset=(Phenotype=="cx3645dKO"),select=c("distance","TC"))
B3_split=split(B3,as.factor(B3$distance))
distance_B3=as.numeric(names(B3_split))

WT_sim=subset(r_sim,subset=(Phenotype=="WT"),select=c("distance","TC"))
WT_sim_split=split(WT_sim,as.factor(WT_sim$distance))
distance_WT_sim=as.numeric(names(WT_sim_split))

B2_sim=subset(r_sim,subset=(Phenotype=="B2TG"),select=c("distance","TC"))
B2_sim_split=split(B2_sim,as.factor(B2_sim$distance))
distance_sim_B2=as.numeric(names(B2_sim_split))

WT_split=split(WT,as.factor(WT$distance))
WT_box=rbind(WT_split[[1]],WT_split[[2]],WT_split[[4]],WT_split[[7]],WT_split[[10]],WT_split[[14]],WT_split[[19]],WT_split[[24]],WT_split[[32]])
WT_sim_split=split(WT_sim,as.factor(WT_sim$distance))
WT_sim_box=rbind(WT_sim_split[[1]],WT_sim_split[[2]],WT_sim_split[[4]],WT_sim_split[[7]],WT_sim_split[[10]],WT_sim_split[[14]],WT_sim_split[[19]],WT_sim_split[[24]],WT_sim_split[[32]])

B2_split=split(B2,as.factor(B2$distance))
B2_box=rbind(B2_split[[1]],B2_split[[2]],B2_split[[4]],B2_split[[7]],B2_split[[10]],B2_split[[14]],B2_split[[19]],B2_split[[24]],B2_split[[32]])
B2_sim_split=split(B2_sim,as.factor(B2_sim$distance))
B2_sim_box=rbind(B2_sim_split[[1]],B2_sim_split[[2]],B2_sim_split[[4]],B2_sim_split[[7]],B2_sim_split[[10]],B2_sim_split[[14]],B2_sim_split[[19]],B2_sim_split[[24]],B2_sim_split[[32]])

B3_split=split(B3,as.factor(B3$distance))
B3_box=rbind(B3_split[[1]],B3_split[[2]],B3_split[[4]],B3_split[[7]],B3_split[[10]],B3_split[[14]],B3_split[[19]],B3_split[[24]],B3_split[[32]])
B3_sim_split=split(B3_sim,as.factor(B3_sim$distance))
B3_sim_box=rbind(B3_sim_split[[1]],B3_sim_split[[2]],B3_sim_split[[4]],B3_sim_split[[7]],B3_sim_split[[10]],B3_sim_split[[14]],B3_sim_split[[19]],B3_sim_split[[24]],B3_sim_split[[32]])



#Create our matrices which we can split up
F=r
F[,2]=y_sim[1,]


#Get the medians by recording by distance

# data
m_rd=mat.or.vec(14,length(distances))
m_rd_F=mat.or.vec(14,length(distances))


for(i in 1:14){

recording=subset(r,subset=(Animal==i),select=c("distance","TC"))
recording_split=split(recording,as.factor(recording$distance))
distance_recording=as.numeric(names(recording_split))

for(j in 1:length(distance_recording)){
m_rd[i,j]=median(recording_split[[j]][,2])
}

recording=subset(F,subset=(Animal==i),select=c("distance","TC"))
recording_split=split(recording,as.factor(recording$distance))
distance_recording=as.numeric(names(recording_split))

for(j in 1:length(distance_recording)){
m_rd_F[i,j]=median(recording_split[[j]][,2])
}

}

##################################################################################################################################################################################################
#now q1


# data
q1_rd=mat.or.vec(14,length(distances))
q1_rd_F=mat.or.vec(14,length(distances))


for(i in 1:14){

recording=subset(r,subset=(Animal==i),select=c("distance","TC"))
recording_split=split(recording,as.factor(recording$distance))
distance_recording=as.numeric(names(recording_split))

for(j in 1:length(distance_recording)){
q1_rd[i,j]=quantile(recording_split[[j]][,2],0.25)
}

recording=subset(F,subset=(Animal==i),select=c("distance","TC"))
recording_split=split(recording,as.factor(recording$distance))
distance_recording=as.numeric(names(recording_split))

for(j in 1:length(distance_recording)){
q1_rd_F[i,j]=quantile(recording_split[[j]][,2],0.25)
}

}


##################################################################################################################################################################################################
#now q3


# data
q3_rd=mat.or.vec(14,length(distances))
q3_rd_F=mat.or.vec(14,length(distances))


for(i in 1:14){

recording=subset(r,subset=(Animal==i),select=c("distance","TC"))
recording_split=split(recording,as.factor(recording$distance))
distance_recording=as.numeric(names(recording_split))

for(j in 1:length(distance_recording)){
q3_rd[i,j]=quantile(recording_split[[j]][,2],0.75)
}

recording=subset(F,subset=(Animal==i),select=c("distance","TC"))
recording_split=split(recording,as.factor(recording$distance))
distance_recording=as.numeric(names(recording_split))

for(j in 1:length(distance_recording)){
q3_rd_F[i,j]=quantile(recording_split[[j]][,2],0.75)
}

}

######################################################################################################################################################################################
#######PDM fits
N=2000


sim=extract(fit,permuted=T)


chai_stat_lv2_alpha=rep(0,N)
chai_stat_lv2_beta=rep(0,N)
chai_stat_lv2_c=rep(0,N)
chai_stat_lv2_A=rep(0,N)
chai_stat_lv2_B=rep(0,N)
chai_stat_lv2_C=rep(0,N)
#note that this can very easily be sped up, but it will not be v slow as it is 

d_alpha=rep(0,rec)
d_beta=rep(0,rec)
d_c=rep(0,rec)
d_A=rep(0,rec)
d_B=rep(0,rec)
d_C=rep(0,rec)

for(k in 1:N){


#we sample from each recording in turn according to order to create a vector of alpha, beta etc
#means our loop is much smaller 

alpha=rep(0,rec)
theta_a=rep(0,rec)
sd_a=rep(0,rec)

theta_a[1:nWT]=sample(sim$theta_a[,1],nWT,replace=TRUE)
sd_a[1:nWT]=sample(sim$sd_a[,1],nWT,replace=TRUE)
theta_a[(nWT+1):(nWT+nKO)]=sample(sim$theta_a[,2],(nKO),replace=TRUE)
sd_a[(nWT+1):(nWT+nKO)]=sample(sim$sd_a[,2],(nKO),replace=TRUE)
theta_a[(nWT+nKO+1):rec]=sample(sim$theta_a[,3],(rec-nWT-nKO),replace=TRUE)
sd_a[(nWT+nKO+1):rec]=sample(sim$sd_a[,3],(rec-nWT-nKO),replace=TRUE)

beta=rep(0,rec)
theta_b=rep(0,rec)
sd_b=rep(0,rec)

theta_b[1:nWT]=sample(sim$theta_b[,1],nWT,replace=TRUE)
sd_b[1:nWT]=sample(sim$sd_b[,1],nWT,replace=TRUE)
theta_b[(nWT+1):(nWT+nKO)]=sample(sim$theta_b[,2],(nKO),replace=TRUE)
sd_b[(nWT+1):(nWT+nKO)]=sample(sim$sd_b[,2],(nKO),replace=TRUE)
theta_b[(nWT+nKO+1):rec]=sample(sim$theta_b[,3],(rec-nWT-nKO),replace=TRUE)
sd_b[(nWT+nKO+1):rec]=sample(sim$sd_b[,3],(rec-nWT-nKO),replace=TRUE)


c=rep(0,rec)
theta_c=rep(0,rec)
sd_c=rep(0,rec)

theta_c[1:nWT]=sample(sim$theta_c[,1],nWT,replace=TRUE)
sd_c[1:nWT]=sample(sim$sd_c[,1],nWT,replace=TRUE)
theta_c[(nWT+1):(nWT+nKO)]=sample(sim$theta_c[,2],(nKO),replace=TRUE)
sd_c[(nWT+1):(nWT+nKO)]=sample(sim$sd_c[,2],(nKO),replace=TRUE)
theta_c[(nWT+nKO+1):rec]=sample(sim$theta_c[,3],(rec-nWT-nKO),replace=TRUE)
sd_c[(nWT+nKO+1):rec]=sample(sim$sd_c[,3],(rec-nWT-nKO),replace=TRUE)


A=rep(0,rec)
theta_A=rep(0,rec)
sd_A=rep(0,rec)

theta_A[1:nWT]=sample(sim$theta_A[,1],nWT,replace=TRUE)
sd_A[1:nWT]=sample(sim$sd_A[,1],nWT,replace=TRUE)
theta_A[(nWT+1):(nWT+nKO)]=sample(sim$theta_A[,2],(nKO),replace=TRUE)
sd_A[(nWT+1):(nWT+nKO)]=sample(sim$sd_A[,2],(nKO),replace=TRUE)
theta_A[(nWT+nKO+1):rec]=sample(sim$theta_A[,3],(rec-nWT-nKO),replace=TRUE)
sd_A[(nWT+nKO+1):rec]=sample(sim$sd_A[,3],(rec-nWT-nKO),replace=TRUE)


B=rep(0,rec)
theta_B=rep(0,rec)
sd_B=rep(0,rec)

theta_B[1:nWT]=sample(sim$theta_B[,1],nWT,replace=TRUE)
sd_B[1:nWT]=sample(sim$sd_B[,1],nWT,replace=TRUE)
theta_B[(nWT+1):(nWT+nKO)]=sample(sim$theta_B[,2],(nKO),replace=TRUE)
sd_B[(nWT+1):(nWT+nKO)]=sample(sim$sd_B[,2],(nKO),replace=TRUE)
theta_B[(nWT+nKO+1):rec]=sample(sim$theta_B[,3],(rec-nWT-nKO),replace=TRUE)
sd_B[(nWT+nKO+1):rec]=sample(sim$sd_B[,3],(rec-nWT-nKO),replace=TRUE)


C=rep(0,rec)
theta_C=rep(0,rec)
sd_C=rep(0,rec)

theta_C[1:nWT]=sample(sim$theta_C[,1],nWT,replace=TRUE)
sd_C[1:nWT]=sample(sim$sd_C[,1],nWT,replace=TRUE)
theta_C[(nWT+1):(nWT+nKO)]=sample(sim$theta_C[,2],(nKO),replace=TRUE)
sd_C[(nWT+1):(nWT+nKO)]=sample(sim$sd_C[,2],(nKO),replace=TRUE)
theta_C[(nWT+nKO+1):rec]=sample(sim$theta_C[,3],(rec-nWT-nKO),replace=TRUE)
sd_C[(nWT+nKO+1):rec]=sample(sim$sd_C[,3],(rec-nWT-nKO),replace=TRUE)


for(i in 1:rec){
alpha[i]=sample(sim$alpha[,i],1,replace=TRUE)
beta[i]=sample(sim$beta[,i],1,replace=TRUE)
c[i]=sample(sim$c[,i],1,replace=TRUE)
A[i]=sample(sim$A[,i],1,replace=TRUE)
B[i]=sample(sim$B[,i],1,replace=TRUE)
C[i]=sample(sim$C[,i],1,replace=TRUE)
d_alpha[i]=(alpha[i]-theta_a[i])/sd_a[i]
d_beta[i]=(beta[i]-theta_b[i])/sd_b[i]
d_c[i]=(c[i]-theta_c[i])/sd_c[i]
d_A[i]=(A[i]-theta_A[i])/sd_A[i]
d_B[i]=(B[i]-theta_B[i])/sd_B[i]
d_C[i]=(C[i]-theta_C[i])/sd_C[i]
}

#remeber that stan parameterises a normal distribution in terms of s.d. not variance 

chai_stat_lv2_alpha[k]=sum(d_alpha^2)
chai_stat_lv2_beta[k]=sum(d_beta^2)
chai_stat_lv2_c[k]=sum(d_c^2)
chai_stat_lv2_A[k]=sum(d_A^2)
chai_stat_lv2_B[k]=sum(d_B^2)
chai_stat_lv2_C[k]=sum(d_C^2)

}
###############################################################################################################################################
#Plot perturbed priors:


library(pscl)
a=read_stan_csv(c("Bshp_F_gb_v4_1.csv","Bshp_F_gb_v4_2.csv","Bshp_F_gb_v4_3.csv"))
b=read_stan_csv(c("Bshp_F_gb_norm_conj_v2_1.csv","Bshp_F_gb_norm_conj_v2_2.csv","Bshp_F_gb_norm_conj_v2_3.csv"))
c=read_stan_csv(c("Bshp_F_v4_1.csv","Bshp_F_v4_2.csv","Bshp_F_v4_3.csv"))


a1=extract(a,permuted=TRUE)
b1=extract(b,permuted=TRUE)
c1=extract(c,permuted=TRUE)




#################################################################################################################################################################################################
#Big graph for paper 


pdf(file="Figure_11.pdf",width=10,height=10)
attach(mtcars)
a=layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),6,3,byrow=TRUE))
par(mar=c(4,4,1,1))
par(oma=c(1,3,6,1))




#########################################################################################################################################################################
#BOX PLOT - wild type
#########################################################################################################################################################################
range_wt=range(c(WT$TC,WT_sim$TC))
range_b2=range(c(B2$TC,B2_sim$TC))

d_points=c(0,100,200,300,400,500,600,700,860)


bplot(WT_box$TC,by=WT_box$distance,pos=d_points,xaxt="n",col="grey",main="",pch=17, outcol="grey", ylim=range_wt,cex=0.4,ylab="STTC",xlim=c(0,900),boxwex=0.2)
axis(1,at=c(0,100,200,300,400,500,600,700,800),labels=FALSE)
bplot(WT_sim_box$TC,by=WT_sim_box$distance,pos=d_points+25,add=TRUE,col=c("green"),pch=17,outcol="green",cex=0.4,xaxt="n",boxwex=0.2)
abline(h=1,lty=2)
legend("topright",pch=17,col=c("grey","green"),legend=c("data","model"),bg="white")


bplot(B2_box$TC,by=B2_box$distance,pos=d_points,xaxt="n",col="grey",main="",pch=17, outcol="grey", ylim=range_b2,cex=0.4,ylab="",xlab=expression(paste("Electrode separation (",mu,"m)")),boxwex=0.2)
axis(1,at=c(0,100,200,300,400,500,600,700,800),labels=TRUE)
bplot(B2_sim_box$TC,by=B2_sim_box$distance,pos=d_points+25,add=TRUE,col=c("green"),pch=17,outcol="green",cex=0.2,xaxt="n",boxwex=0.2)
abline(h=1,lty=2)

bplot(B3_box$TC,by=B3_box$distance,pos=d_points,xaxt="n",col="grey",main="",pch=17, outcol="grey", ylim=range_b3,cex=0.4,ylab="",xlab=expression(paste("Electrode separation (",mu,"m)")),boxwex=0.2)
axis(1,at=c(0,100,200,300,400,500,600,700,800),labels=TRUE)
bplot(B3_sim_box$TC,by=B3_sim_box$distance,pos=d_points+25,add=TRUE,col=c("green"),pch=17,outcol="green",cex=0.2,xaxt="n",boxwex=0.2)
abline(h=1,lty=2)



####################################################################

plot(distance_recording,m_rd[1,],type="l",ylim=range(c(q1_rd[1,],q3_rd_F[1,])),ylab="STTC",xlab="",lwd=1.5,xaxt="n",lty=5)
axis(1,at=c(0,100,200,300,400,500,600,700,800),labels=FALSE)
points(distance_recording,m_rd_F[1,],type="l",lty=1,lwd=1.5)
polygon(c(distance_recording,rev(distance_recording)),c(q1_rd[1,],rev(q3_rd[1,])),col=rgb(0.192,0.192,0.192,0.2),border=NA)
arrows(as.numeric(as.vector(distance_recording)),as.numeric(as.vector(q1_rd_F[1,])), as.numeric(as.vector(distance_recording)),as.numeric(as.vector(q3_rd_F[1,])), angle=90, code=3, length=0.03,lwd=1)

plot(distance_recording,m_rd[8,],type="l",ylim=range(c(q1_rd[8,],q3_rd_F[8,])),ylab="",xlab="",col="red",lwd=1.5,xaxt="n",lty=5)
axis(1,at=c(0,100,200,300,400,500,600,700,800),labels=FALSE)
points(distance_recording,m_rd_F[8,],type="l",lty=1,col="red",lwd=1.5)
polygon(c(distance_recording,rev(distance_recording)),c(q1_rd[8,],rev(q3_rd[8,])),col=rgb(1,0,0,0.2),border=NA)
arrows(as.numeric(as.vector(distance_recording)),as.numeric(as.vector(q1_rd_F[8,])), as.numeric(as.vector(distance_recording)),as.numeric(as.vector(q3_rd_F[8,])), angle=90, code=3, length=0.03,lwd=1,col="red")


plot(distance_recording,m_rd[11,],type="l",ylim=range(c(q1_rd[11,],q3_rd_F[11,])),ylab="",xlab="",col="blue",lwd=1.5,xaxt="n",lty=5)
axis(1,at=c(0,100,200,300,400,500,600,700,800),labels=FALSE)
points(distance_recording,m_rd_F[27,],type="l",lty=1,col="blue",lwd=1.5)
polygon(c(distance_recording,rev(distance_recording)),c(q1_rd[11,],rev(q3_rd[11,])),col=rgb(0,0,1,0.2),border=NA)
arrows(as.numeric(as.vector(distance_recording)),as.numeric(as.vector(q1_rd_F[11,])), as.numeric(as.vector(distance_recording)),as.numeric(as.vector(q3_rd_F[11,])), angle=90, code=3, length=0.03,lwd=1,col="blue")
legend("topright",lty=c(5,1),lwd=c(1.5,1.5),legend=c("data","model"),col="blue",cex=0.9)

#text(grconvertX(c(0.3, 0.74), from='ndc'),
  # grconvertY(c(0.98, 0.98), from='ndc'),  c("Recording 27", "Recording 10"), xpd=NA, cex=1.5, font=2)
#text(grconvertX(c(0.02, 0.02, 0.02), from='ndc'),
  # grconvertY(c(0.14, 0.33, 0.52,0.72,0.91), from='ndc'),  c("Max","Min","q3","q1","Median"), xpd=NA, cex=1, font=2,srt=90)

#################################################################################

x=seq(from=0,to=100,by=0.1)
hist(chai_stat_lv2_alpha,freq=FALSE,xlim=c(0,80),xlab=expression(mu[alpha]),cex.lab=1.5,xaxt="n",ylab="")
axis(1,at=c(0,20,40,60,80),labels=TRUE)
points(x,dchisq(x,df=28),col="red",type="l")
legend("right",lty=1,col="red",legend=c("Chi-squared df=28"),cex=1,bty="n")
hist(chai_stat_lv2_beta,freq=FALSE,xlim=c(0,80),xlab=expression(mu[beta]),cex.lab=1.5,xaxt="n",ylab="",ylim=c(0,0.06))
axis(1,at=c(0,20,40,60,80),labels=FALSE)
points(x,dchisq(x,df=28),col="red",type="l")
hist(chai_stat_lv2_c,freq=FALSE,xlim=c(0,80),xlab=expression(mu[gamma]),xaxt="n",ylab="",cex.lab=1.5)
axis(1,at=c(0,20,40,60,80),labels=FALSE)
points(x,dchisq(x,df=28),col="red",type="l")
hist(chai_stat_lv2_A,freq=FALSE,xlim=c(0,80),xlab=expression(mu[A]),xaxt="n",ylab="",cex.lab=1.5)
axis(1,at=c(0,20,40,60,80),labels=FALSE)
points(x,dchisq(x,df=28),col="red",type="l")
hist(chai_stat_lv2_B,freq=FALSE,xlim=c(0,80),xlab=expression(mu[B]),xaxt="n",ylab="",ylim=c(0,0.06),cex.lab=1.5)
axis(1,at=c(0,20,40,60,80),labels=FALSE)
points(x,dchisq(x,df=28),col="red",type="l")
hist(chai_stat_lv2_C,freq=FALSE,xlim=c(0,80),xlab=expression(mu[C]),xaxt="n",ylab="",cex.lab=1.5,ylim=c(0,0.06))
axis(1,at=c(0,20,40,60,80),labels=FALSE)
points(x,dchisq(x,df=28),col="red",type="l")

###############################################################################

plot(density(a1$theta_a[,1]),xlim=c(-0.05,0.15),ylim=range(c(density(a1$theta_a[,1])$y,density(b1$theta_a[,1])$y,density(a1$theta_a[,2])$y,density(b1$theta_a[,2])$y,density(a1$theta_a[,3])$y,density(b1$theta_a[,3])$y)),xlab=expression(mu[alpha]),ylab="density",main="",cex.lab=1.5)
points(density(b1$theta_a[,1]),type="l",lty=2)
points(density(c1$theta_a[,1]),type="l",lty=3)
points(density(a1$theta_a[,2]),type="l",col="red")
points(density(b1$theta_a[,2]),type="l",col="red",lty=2)
points(density(c1$theta_a[,2]),type="l",col="red",lty=3)
points(density(a1$theta_a[,3]),type="l",col="blue")
points(density(b1$theta_a[,3]),type="l",col="blue",lty=2)
points(density(c1$theta_a[,3]),type="l",col="blue",lty=3)

plot(density(a1$theta_b[,1]),xlim=c(0,0.7),ylim=range(c(density(a1$theta_b[,1])$y,density(b1$theta_b[,1])$y,density(a1$theta_b[,2])$y,density(b1$theta_b[,2])$y,density(a1$theta_b[,3])$y,density(b1$theta_b[,3])$y)),xlab=expression(mu[beta]),ylab="",main="",cex.lab=1.5)
points(density(b1$theta_b[,1]),type="l",lty=2)
points(density(c1$theta_b[,1]),type="l",lty=3)
points(density(a1$theta_b[,2]),type="l",col="red")
points(density(b1$theta_b[,2]),type="l",col="red",lty=2)
points(density(c1$theta_b[,2]),type="l",col="red",lty=3)
points(density(a1$theta_b[,3]),type="l",col="blue")
points(density(b1$theta_b[,3]),type="l",col="blue",lty=2)
points(density(c1$theta_b[,3]),type="l",col="blue",lty=3)

plot(density(a1$theta_c[,1]),xlim=c(0,0.015),ylim=range(c(density(a1$theta_c[,1])$y,density(b1$theta_c[,1])$y,density(c1$theta_c[,1])$y,density(a1$theta_c[,2])$y,density(b1$theta_c[,2])$y,density(c1$theta_c[,2])$y,density(a1$theta_c[,3])$y,density(b1$theta_c[,3])$y,density(c1$theta_c[,3])$y)),xlab=expression(mu[gamma]),ylab="",main="",cex.lab=1.5)
points(density(b1$theta_c[,1]),type="l",lty=2)
points(density(c1$theta_c[,1]),type="l",lty=3)
points(density(a1$theta_c[,2]),type="l",col="red")
points(density(b1$theta_c[,2]),type="l",col="red",lty=2)
points(density(c1$theta_c[,2]),type="l",col="red",lty=3)
points(density(a1$theta_c[,3]),type="l",col="blue")
points(density(b1$theta_c[,3]),type="l",col="blue",lty=2)
points(density(c1$theta_c[,3]),type="l",col="blue",lty=3)


plot(density(a1$theta_A[,1]),xlim=c(0,0.15),ylim=range(c(density(a1$theta_A[,1])$y,density(b1$theta_A[,1])$y,density(c1$theta_A[,1])$y,density(a1$theta_A[,2])$y,density(b1$theta_A[,2])$y,density(c1$theta_A[,2])$y,density(a1$theta_A[,3])$y,density(b1$theta_A[,3])$y,density(c1$theta_A[,3])$y)),xlab=expression(mu[A]),ylab="",main="",cex.lab=1.5)
points(density(b1$theta_A[,1]),type="l",lty=2)
points(density(c1$theta_A[,1]),type="l",lty=3)
points(density(a1$theta_A[,2]),type="l",col="red")
points(density(b1$theta_A[,2]),type="l",col="red",lty=2)
points(density(c1$theta_A[,2]),type="l",col="red",lty=3)
points(density(a1$theta_A[,3]),type="l",col="blue")
points(density(b1$theta_A[,3]),type="l",col="blue",lty=2)
points(density(c1$theta_A[,3]),type="l",col="blue",lty=3)
legend("topright",lty=c(1,3,2),legend=c("Gumbel","Normal","Perturbed priors"))


plot(density(a1$theta_B[,1]),xlim=c(0,0.2),ylim=range(c(density(a1$theta_B[,1])$y,density(b1$theta_B[,1])$y,density(c1$theta_B[,1])$y,density(a1$theta_B[,2])$y,density(b1$theta_B[,2])$y,density(c1$theta_B[,2])$y,density(a1$theta_B[,3])$y,density(b1$theta_B[,3])$y,density(c1$theta_B[,3])$y)),xlab=expression(mu[B]),ylab="",main="",cex.lab=1.5)
points(density(b1$theta_B[,1]),type="l",lty=2)
points(density(c1$theta_B[,1]),type="l",lty=3)
points(density(a1$theta_B[,2]),type="l",col="red")
points(density(b1$theta_B[,2]),type="l",col="red",lty=2)
points(density(c1$theta_B[,2]),type="l",col="red",lty=3)
points(density(a1$theta_B[,3]),type="l",col="blue")
points(density(b1$theta_B[,3]),type="l",col="blue",lty=2)
points(density(c1$theta_B[,3]),type="l",col="blue",lty=3)
legend("topleft",lty=c(1),col=c("black","blue"),legend=c("Wild type",expression(paste(beta,"2TG"))))


plot(density(a1$theta_C[,1]),xlim=c(0,0.00006),ylim=range(c(density(a1$theta_C[,1])$y,density(b1$theta_C[,1])$y,density(a1$theta_C[,2])$y,density(b1$theta_C[,2])$y,density(a1$theta_C[,3])$y,density(b1$theta_C[,3])$y,density(c1$theta_C[,3])$y)),xlab=expression(mu[C]),ylab="",main="",cex.lab=1.5)
points(density(b1$theta_C[,1]),type="l",lty=2)
points(density(c1$theta_C[,1]),type="l",lty=3)
points(density(a1$theta_C[,2]),type="l",col="red")
points(density(b1$theta_C[,2]),type="l",col="red",lty=2)
points(density(c1$theta_C[,2]),type="l",col="red",lty=3)
points(density(a1$theta_C[,3]),type="l",col="blue")
points(density(b1$theta_C[,3]),type="l",col="blue",lty=2)
points(density(c1$theta_C[,3]),type="l",col="blue",lty=3)

text(grconvertX(c(0.03), from='ndc'),
grconvertY(c(0.92, 0.77, 0.61,0.32), from='ndc'),  c("A","B" ,"C" ,"D"), xpd=NA, cex=2, font=2)


#### THINK HOW TO LABEL PHENOTYPES IN TOP TWO ROWS: THEY DO NOT APPLY TO THE BOTTOM FOUR ROWS
text(grconvertX(c(0.27, 0.75), from='ndc'),
  grconvertY(c(0.95,0.95), from='ndc'),  c("Wild type",expression(paste(beta,"2(TG)"))), xpd=NA, cex=1.5, font=2)
dev.off()

