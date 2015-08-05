library(rstan)
library(coda)
posterior=read_stan_csv(c("Bshp_F_gb_v4_1.csv","Bshp_F_gb_v4_2.csv","Bshp_F_gb_v4_3.csv"))

num_points=15000*3

a1=extract(posterior,pars=c("theta_a"),permuted=FALSE)

pdf(file="Figure_10.pdf",width=6,height=8)
attach(mtcars)
a=layout(matrix(c(1,2,3,4,5,6,0,0,7,8,9,10,11,12),7,2,byrow=TRUE), widths=c(1,1), heights=c(4,4,4,1,4,4,4))
par(oma=c(0.5,3.5,4.5,0.5))
par(mar=c(3.5,3,1,1))
par(mgp=c(1.75,0.5,0))


plot(density(a1[,,1]),xlim=c(-0.11,0.15),ylim=range(density(a1[,,1])$y,density(a1[,,2])$y,density(a1[,,3])$y),xlab=expression(mu[alpha]),ylab="density",main="",cex.lab=1.5)
points(density(a1[,,2]),type="l",col="red")
points(density(a1[,,3]),type="l",col="blue")

q=HPDinterval(as.mcmc(as.vector(a1[,,1])),prob=0.95)
abline(v=q[1],lty=2)
abline(v=q[2],lty=2)
q=HPDinterval(as.mcmc(as.vector(a1[,,2])),prob=0.95)
abline(v=q[1],lty=2,col="red")
abline(v=q[2],lty=2,col="red")
q=HPDinterval(as.mcmc(as.vector(a1[,,3])),prob=0.95)
abline(v=q[1],lty=2,col="blue")
abline(v=q[2],lty=2,col="blue")
#legend("topright",expression(mu[alpha]),cex=2,bty="n")




a1=extract(posterior,pars=c("theta_A"),permuted=FALSE)
plot(density(a1[,,1]),xlim=c(0,0.15),ylim=range(density(a1[,,1])$y,density(a1[,,2])$y,density(a1[,,3])$y),xlab=expression(mu[A]),ylab="",main="",cex.lab=1.5)
points(density(a1[,,2]),type="l",col="red")
points(density(a1[,,3]),type="l",col="blue")
q=HPDinterval(as.mcmc(as.vector(a1[,,1])),prob=0.95)
abline(v=q[1],lty=2)
abline(v=q[2],lty=2)
q=HPDinterval(as.mcmc(as.vector(a1[,,2])),prob=0.95)
abline(v=q[1],lty=2,col="red")
abline(v=q[2],lty=2,col="red")
q=HPDinterval(as.mcmc(as.vector(a1[,,3])),prob=0.95)
abline(v=q[1],lty=2,col="blue")
abline(v=q[2],lty=2,col="blue")

#legend("topright",expression(mu[A]),cex=2,bty="n")






a1=extract(posterior,pars=c("theta_b"),permuted=FALSE)
plot(density(a1[,,1]),xlim=c(0,1),ylim=range(density(a1[,,1])$y,density(a1[,,2])$y,density(a1[,,3])$y),xlab=expression(mu[beta]),ylab="",main="",cex.lab=1.5)
points(density(a1[,,2]),type="l",col="red")
points(density(a1[,,3]),type="l",col="blue")
q=HPDinterval(as.mcmc(as.vector(a1[,,1])),prob=0.95)
abline(v=q[1],lty=2)
abline(v=q[2],lty=2)
q=HPDinterval(as.mcmc(as.vector(a1[,,2])),prob=0.95)
abline(v=q[1],lty=2,col="red")
abline(v=q[2],lty=2,col="red")
q=HPDinterval(as.mcmc(as.vector(a1[,,3])),prob=0.95)
abline(v=q[1],lty=2,col="blue")
abline(v=q[2],lty=2,col="blue")


#legend("topright",expression(mu[beta]),cex=2,bty="n")

a1=extract(posterior,pars=c("theta_B"),permuted=FALSE)
plot(density(a1[,,1]),xlim=c(0,0.25),ylim=range(density(a1[,,1])$y,density(a1[,,2])$y,density(a1[,,3])$y),xlab=expression(mu[B]),ylab="",main="",cex.lab=1.5)
points(density(a1[,,2]),type="l",col="red")
points(density(a1[,,3]),type="l",col="blue")
q=HPDinterval(as.mcmc(as.vector(a1[,,1])),prob=0.95)
abline(v=q[1],lty=2)
abline(v=q[2],lty=2)
q=HPDinterval(as.mcmc(as.vector(a1[,,2])),prob=0.95)
abline(v=q[1],lty=2,col="red")
abline(v=q[2],lty=2,col="red")
q=HPDinterval(as.mcmc(as.vector(a1[,,3])),prob=0.95)
abline(v=q[1],lty=2,col="blue")
abline(v=q[2],lty=2,col="blue")
#legend("topright",expression(mu[B]),cex=2,bty="n")


a1=extract(posterior,pars=c("theta_c"),permuted=FALSE)
plot(density(a1[,,1]),xlim=c(0,0.012),ylim=range(density(a1[,,1])$y,density(a1[,,2])$y,density(a1[,,3])$y),xlab=expression(mu[gamma]),ylab="",main="",cex.lab=1.5)
points(density(a1[,,2]),type="l",col="red")
points(density(a1[,,3]),type="l",col="blue")
q=HPDinterval(as.mcmc(as.vector(a1[,,1])),prob=0.95)
abline(v=q[1],lty=2)
abline(v=q[2],lty=2)
q=HPDinterval(as.mcmc(as.vector(a1[,,2])),prob=0.95)
abline(v=q[1],lty=2,col="red")
abline(v=q[2],lty=2,col="red")
q=HPDinterval(as.mcmc(as.vector(a1[,,3])),prob=0.95)
abline(v=q[1],lty=2,col="blue")
abline(v=q[2],lty=2,col="blue")
#legend("topright",expression(mu[gamma]),cex=2,bty="n")




a1=extract(posterior,pars=c("theta_C"),permuted=FALSE)
plot(density(a1[,,1]),xlim=c(0,0.00005),ylim=range(density(a1[,,1])$y,density(a1[,,2])$y,density(a1[,,3])$y),xlab=expression(mu[C]),ylab="",main="",cex.lab=1.5)
points(density(a1[,,2]),type="l",col="red")
points(density(a1[,,3]),type="l",col="blue")
q=HPDinterval(as.mcmc(as.vector(a1[,,1])),prob=0.95)
abline(v=q[1],lty=2)
abline(v=q[2],lty=2)
q=HPDinterval(as.mcmc(as.vector(a1[,,2])),prob=0.95)
abline(v=q[1],lty=2,col="red")
abline(v=q[2],lty=2,col="red")
q=HPDinterval(as.mcmc(as.vector(a1[,,3])),prob=0.95)
abline(v=q[1],lty=2,col="blue")
abline(v=q[2],lty=2,col="blue")
#legend("topright",expression(mu[C]),cex=2,bty="n")


a1=extract(posterior,pars=c("sd_a"),permuted=FALSE)
plot(density(a1[,,1]),xlim=c(0,0.2),ylim=range(density(a1[,,1])$y,density(a1[,,2])$y,density(a1[,,3])$y),xlab=expression(sigma[alpha]),ylab="",main="",cex.lab=1.5)
points(density(a1[,,2]),type="l",col="red")
points(density(a1[,,3]),type="l",col="blue")
q=HPDinterval(as.mcmc(as.vector(a1[,,1])),prob=0.95)
abline(v=q[1],lty=2)
abline(v=q[2],lty=2)
q=HPDinterval(as.mcmc(as.vector(a1[,,2])),prob=0.95)
abline(v=q[1],lty=2,col="red")
abline(v=q[2],lty=2,col="red")
q=HPDinterval(as.mcmc(as.vector(a1[,,3])),prob=0.95)
abline(v=q[1],lty=2,col="blue")
abline(v=q[2],lty=2,col="blue")
#legend(x="bottomright",lty=1,col=c("black","blue"),legend=c("Wild type",expression(paste(beta,"2(TG)"))))
#legend(x="right",lty=2,col=c("black"),legend=c("95 % HDR"))
#legend("topright",expression(sigma[alpha]),cex=2,bty="n")


a1=extract(posterior,pars=c("sd_A"),permuted=FALSE)
plot(density(a1[,,1]),xlim=c(0,0.2),ylim=range(density(a1[,,1])$y,density(a1[,,2])$y,density(a1[,,3])$y),xlab=expression(sigma[A]),ylab="",main="",cex.lab=1.5)
points(density(a1[,,2]),type="l",col="red")
points(density(a1[,,3]),type="l",col="blue")
q=HPDinterval(as.mcmc(as.vector(a1[,,1])),prob=0.95)
abline(v=q[1],lty=2)
abline(v=q[2],lty=2)
q=HPDinterval(as.mcmc(as.vector(a1[,,2])),prob=0.95)
abline(v=q[1],lty=2,col="red")
abline(v=q[2],lty=2,col="red")
q=HPDinterval(as.mcmc(as.vector(a1[,,3])),prob=0.95)
abline(v=q[1],lty=2,col="blue")
abline(v=q[2],lty=2,col="blue")#legend("topright",expression(sigma[A]),cex=2,bty="n")




a1=extract(posterior,pars=c("sd_b"),permuted=FALSE)
plot(density(a1[,,1]),xlim=c(0,1),ylim=range(density(a1[,,1])$y,density(a1[,,2])$y,density(a1[,,3])$y),xlab=expression(sigma[beta]),ylab="",main="",cex.lab=1.5)
points(density(a1[,,2]),type="l",col="red")
points(density(a1[,,3]),type="l",col="blue")
q=HPDinterval(as.mcmc(as.vector(a1[,,1])),prob=0.95)
abline(v=q[1],lty=2)
abline(v=q[2],lty=2)
q=HPDinterval(as.mcmc(as.vector(a1[,,2])),prob=0.95)
abline(v=q[1],lty=2,col="red")
abline(v=q[2],lty=2,col="red")
q=HPDinterval(as.mcmc(as.vector(a1[,,3])),prob=0.95)
abline(v=q[1],lty=2,col="blue")
abline(v=q[2],lty=2,col="blue")
#legend("topright",expression(sigma[beta]),cex=2,bty="n")




a1=extract(posterior,pars=c("sd_B"),permuted=FALSE)
plot(density(a1[,,1]),xlim=c(0,0.25),ylim=range(density(a1[,,1])$y,density(a1[,,2])$y,density(a1[,,3])$y),xlab=expression(sigma[B]),ylab="",main="",cex.lab=1.5)
points(density(a1[,,2]),type="l",col="red")
points(density(a1[,,3]),type="l",col="blue")
q=HPDinterval(as.mcmc(as.vector(a1[,,1])),prob=0.95)
abline(v=q[1],lty=2)
abline(v=q[2],lty=2)
q=HPDinterval(as.mcmc(as.vector(a1[,,2])),prob=0.95)
abline(v=q[1],lty=2,col="red")
abline(v=q[2],lty=2,col="red")
q=HPDinterval(as.mcmc(as.vector(a1[,,3])),prob=0.95)
abline(v=q[1],lty=2,col="blue")
abline(v=q[2],lty=2,col="blue")
#legend("topright",expression(sigma[B]),cex=2,bty="n")
#legend(x="topright",lty=1,col=c("black","blue"),legend=c("Wild type",expression(paste(beta,"2(TG)"))))


a1=extract(posterior,pars=c("sd_c"),permuted=FALSE)
plot(density(a1[,,1]),xlim=c(0,0.012),ylim=range(density(a1[,,1])$y,density(a1[,,2])$y,density(a1[,,3])$y),xlab=expression(sigma[gamma]),ylab="",main="",cex.lab=1.5)
points(density(a1[,,2]),type="l",col="red")
points(density(a1[,,3]),type="l",col="blue")
q=HPDinterval(as.mcmc(as.vector(a1[,,1])),prob=0.95)
abline(v=q[1],lty=2)
abline(v=q[2],lty=2)
q=HPDinterval(as.mcmc(as.vector(a1[,,2])),prob=0.95)
abline(v=q[1],lty=2,col="red")
abline(v=q[2],lty=2,col="red")
q=HPDinterval(as.mcmc(as.vector(a1[,,3])),prob=0.95)
abline(v=q[1],lty=2,col="blue")
abline(v=q[2],lty=2,col="blue")
 #legend("topright",expression(sigma[gamma]),cex=2,bty="n")





a1=extract(posterior,pars=c("sd_C"),permuted=FALSE)
plot(density(a1[,,1]),xlim=c(0,0.00004),ylim=range(density(a1[,,1])$y,density(a1[,,2])$y,density(a1[,,3])$y),xlab=expression(sigma[C]),ylab="",main="",cex.lab=1.5)
points(density(a1[,,2]),type="l",col="red")
points(density(a1[,,3]),type="l",col="blue")
q=HPDinterval(as.mcmc(as.vector(a1[,,1])),prob=0.95)
abline(v=q[1],lty=2)
abline(v=q[2],lty=2)
q=HPDinterval(as.mcmc(as.vector(a1[,,2])),prob=0.95)
abline(v=q[1],lty=2,col="red")
abline(v=q[2],lty=2,col="red")
q=HPDinterval(as.mcmc(as.vector(a1[,,3])),prob=0.95)
abline(v=q[1],lty=2,col="blue")
abline(v=q[2],lty=2,col="blue")
#legend("topright",expression(sigma[C]),cex=2,bty="n")
legend(x="topright",lty=c(1,1,1,2),col=c("black","red","blue","black"),legend=c("Wild type","Cx45ko","Cx36/45dko","95% HDR"),cex=0.85,bty="n")




text(grconvertX(c(0.35,0.80), from='ndc'),
   grconvertY(c(0.97,0.97), from='ndc'),  c('Location parameters','Scale parameters'), xpd=NA, cex=1.5, font=2)

text(grconvertX(c(0.55), from='ndc'),
   grconvertY(c(0.93), from='ndc'),  c('Mean'), xpd=NA, cex=1.25, font=2)

text(grconvertX(c(0.55), from='ndc'),
   grconvertY(c(0.47), from='ndc'),  c('Standard deviation'), xpd=NA, cex=1.25, font=2)


text(grconvertX(c(0.03,0.03,0.03), from='ndc'),
   grconvertY(c(0.10,0.25,0.39,0.58,0.72,0.87), from='ndc'),  c('Gradient','Range','Baseline'), xpd=NA, cex=1.5, font=2,srt=90,ps=0.6)


dev.off()





