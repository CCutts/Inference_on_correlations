data{

int<lower=0> n; 
int<lower=0> t; 
vector[n] x; 
vector[n] y; 

vector[n] z;
int Rec[n];
int<lower=0> nWT;
int<lower=0> nKO;
int<lower=0> rec;
}


parameters{

real<lower=-1,upper=1> alpha[rec];
real<lower=0,upper=2> beta[rec];
real<lower=0,upper=0.5> c[rec];
real<lower=0,upper=1> A[rec];
real<lower=0,upper=1> B[rec];
real<lower=0,upper=0.5> C[rec]; 
real<lower=-1,upper=1> theta_a[3];
real<lower=0,upper=2> theta_b[3];
real<lower=0,upper=0.5> theta_c[3];
real<lower=0,upper=1> sd_a[3];
real<lower=0,upper=1> sd_b[3];
real<lower=0,upper=0.25> sd_c[3];
real<lower=0,upper=1> theta_A[3];
real<lower=0,upper=1> theta_B[3];
real<lower=0,upper=0.5> theta_C[3];
real<lower=0,upper=0.25 > sd_A[3];
real<lower=0,upper=0.25 > sd_B[3];
real<lower=0, upper=0.25 > sd_C[3];
}



model {
vector[n] Ymu;
vector[n] Ysd;
int j;	
		
		for(i in 1:nWT) {
		alpha[i]~ normal(theta_a[1],sd_a[1]) T[-1,1];
		beta[i]~normal(theta_b[1],sd_b[1]) T[0,2];
		c[i]~normal(theta_c[1],sd_c[1]) T[0,0.5];
		A[i]~normal(theta_A[1],sd_A[1])T[0,1];
		B[i]~normal(theta_B[1],sd_B[1])T[0,1];
		C[i]~normal(theta_C[1],sd_C[1])T[0,0.5];
		}

		for(i in (nWT+1):(nWT+nKO)){
			alpha[i]~normal(theta_a[2],sd_a[2]) T[-1,1];
			beta[i]~normal(theta_b[2],sd_b[2]) T[0,2];
			c[i]~normal(theta_c[2],sd_c[2]) T[0,0.5];
			A[i]~normal(theta_A[2],sd_A[2])T[0,1];
			B[i]~normal(theta_B[2],sd_B[2])T[0,1];
			C[i]~normal(theta_C[2],sd_C[2])T[0,0.5];		
		}
		for(i in (nWT+nKO+1):rec){
			alpha[i]~normal(theta_a[3],sd_a[3]) T[-1,1];
			beta[i]~normal(theta_b[3],sd_b[3]) T[0,2];
			c[i]~normal(theta_c[3],sd_c[3]) T[0,0.5];
			A[i]~normal(theta_A[3],sd_A[3])T[0,1];
			B[i]~normal(theta_B[3],sd_B[3])T[0,1];
			C[i]~normal(theta_C[3],sd_C[3])T[0,0.5];		
		}

		

		for(p in 1:n) {
		j<-Rec[p];
		Ymu[p]<-alpha[j]+beta[j]*exp(-c[j]*x[p]);
		Ysd[p]<-A[j]+B[j]*exp(-C[j]*z[p]);		
		}
		
		
		
		y~gumbel(Ymu,Ysd);
		
	
		
}

generated quantities{
vector[n] log_lik;
int<lower=0> q;

for(j in 1:n){
q<-Rec[j];
log_lik[j]<- gumbel_log(y[j],alpha[q]+beta[q]*exp(-c[q]*x[j]),A[q]+B[q]*exp(-C[q]*z[j]));
}
}



