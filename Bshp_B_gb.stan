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
real<lower=-1,upper=1> theta_a[1];
real<lower=0,upper=2> theta_b[1];
real<lower=0,upper=1> theta_c[1];
real<lower=0,upper=1> sd_a[1];
real<lower=0,upper=1> sd_b[1];
real<lower=0,upper=0.25> sd_c[1];
real<lower=0,upper=1> theta_A[1];
real<lower=0,upper=1> theta_B[1];
real<lower=0,upper=0.5> theta_C[1];
real<lower=0,upper=0.25 > sd_A[1];
real<lower=0,upper=0.25 > sd_B[1];
real<lower=0, upper=0.25 > sd_C[1];
}






model {
vector[n] Ymu;
vector[n] Ysd;
int j;	
		
for(p in 1:rec){	
		alpha[p]~normal(theta_a[1],sd_a[1]) T[-1,1];
		beta[p]~normal(theta_b[1],sd_b[1]) T[0,2];
		c[p]~normal(theta_c[1],sd_c[1]) T[0,0.5];
		A[p]~normal(theta_A[1],sd_A[1])T[0,1];
		B[p]~normal(theta_B[1],sd_B[1])T[0,1];
		C[p]~normal(theta_C[1],sd_C[1])T[0,0.5];
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
log_lik[j]<-gumbel_log(y[j],alpha[q]+beta[q]*exp(-c[q]*x[j]),A[q]+B[q]*exp(-C[q]*z[j]));
}
}

