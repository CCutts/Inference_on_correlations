ata{

int<lower=0> n; 
int<lower=0> t; 
vector[n] x; 
vector[n] y; 

vector[n] z;
int Rec[n];
int<lower=0> nWT;
int<lower=0> rec;
}


parameters{

real<lower=-1,upper=1> alpha[rec];
real<lower=0,upper=2> beta[rec];
real<lower=0,upper=0.5> c[rec];
real<lower=0,upper=1> A[rec];
real<lower=0,upper=1> B[rec];
real<lower=0,upper=0.5> C[rec]; 
real<lower=-1,upper=1> theta_a[2];
real<lower=0,upper=2> theta_b[2];
real<lower=0,upper=0.5> theta_c[2];
real<lower=0,upper=1> theta_A[2];
real<lower=0,upper=1> theta_B[2];
real<lower=0,upper=0.5> theta_C[2];
real<lower=0,upper=1> sd_a[2];
real<lower=0,upper=1> sd_b[2];
real<lower=0,upper=0.25> sd_c[2];
real<lower=0,upper=0.25 > sd_A[2];
real<lower=0,upper=0.25 > sd_B[2];
real<lower=0, upper=0.25 > sd_C[2];

}


model {
vector[n] Ymu;
vector[n] Ysd;
int j;	
		
		theta_a[1] ~ normal(0,0.5) T[-1,1];
		theta_a[2] ~ normal(0,0.5) T[-1,1];
		theta_b[1] ~ normal(1,0.5) T[0,2];
		theta_b[2] ~ normal(1,0.5) T[0,2];
		theta_c[1] ~ normal(0.25,0.25) T[0,0.5];
		theta_c[2] ~ normal(0.25,0.25) T[0,0.5];
		theta_A[1] ~ normal(0.5,0.5) T[0,1];
		theta_A[2] ~ normal(0.5,0.5) T[0,1];
		theta_B[1] ~ normal(0.5,0.5) T[0,1];
		theta_B[2] ~ normal(0.5,0.5) T[0,1];
		theta_C[1] ~ normal(0.25,0.25) T[0,0.5];
		theta_C[2] ~ normal(0.25,0.25) T[0,0.5];
		sd_a ~ inv_gamma(3,0.05);
		sd_b ~ inv_gamma(3,0.05);
		sd_c ~ inv_gamma(3,0.01);
		sd_A ~ inv_gamma(3,0.01);
		sd_B ~ inv_gamma(3,0.01);
		sd_C ~ inv_gamma(3,0.01);

		for(i in 1:nWT) {
		alpha[i]~ normal(theta_a[1],sd_a[1]) T[-1,1];
		beta[i]~normal(theta_b[1],sd_b[1]) T[0,2];
		c[i]~normal(theta_c[1],sd_c[1]) T[0,0.5];
		A[i]~normal(theta_A[1],sd_A[1])T[0,1];
		B[i]~normal(theta_B[1],sd_B[1])T[0,1];
		C[i]~normal(theta_C[1],sd_C[1])T[0,0.5];
		}

		for(i in (nWT+1):rec){
			alpha[i]~normal(theta_a[2],sd_a[2]) T[-1,1];
			beta[i]~normal(theta_b[2],sd_b[2]) T[0,2];
			c[i]~normal(theta_c[2],sd_c[2]) T[0,0.5];
			A[i]~normal(theta_A[2],sd_A[2])T[0,1];
			B[i]~normal(theta_B[2],sd_B[2])T[0,1];
			C[i]~normal(theta_C[2],sd_C[2])T[0,0.5];		
		}
		

		for(p in 1:n) {
		j<-Rec[p];
		Ymu[p]<-alpha[j]+beta[j]*exp(-c[j]*x[p]);
		Ysd[p]<-A[j]+B[j]*exp(-C[j]*z[p]);		
		}
		
		
		
		y~gumbel(Ymu,Ysd);
		
	
		
}



