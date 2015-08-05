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

real<lower=-1,upper=1> theta_a[3];
real<lower=0,upper=2> theta_b[3];
real<lower=0,upper=1> theta_c[3];
real<lower=0,upper=1> theta_A[3];
real<lower=0,upper=1> theta_B[3];
real<lower=0,upper=0.5> theta_C[3];

}


model {
vector[n] Ymu;
vector[n] Ysd;
int j;		
		for(p in 1:n) {
			j<-Rec[p];
			if(j < (nWT+1)){
				Ymu[p]<-theta_a[1]+theta_b[1]*exp(-theta_c[1]*x[p]);
				Ysd[p]<-theta_A[1]+theta_B[1]*exp(-theta_C[1]*z[p]);		
			}
			else if(j<(nWT+nKO+1)){
				Ymu[p]<-theta_a[2]+theta_b[2]*exp(-theta_c[2]*x[p]);
				Ysd[p]<-theta_A[2]+theta_B[2]*exp(-theta_C[2]*z[p]);		
			}
			else{
				Ymu[p]<-theta_a[3]+theta_b[3]*exp(-theta_c[3]*x[p]);
				Ysd[p]<-theta_A[3]+theta_B[3]*exp(-theta_C[3]*z[p]);
			}	
		}
				
		
		y~gumbel(Ymu,Ysd);	
		
}


generated quantities{
vector[n] log_lik;
int<lower=0> q;

for(j in 1:n){
q<-Rec[j];
if(q <(nWT+1)){
log_lik[j]<- gumbel_log(y[j],theta_a[1]+theta_b[1]*exp(-theta_c[1]*x[j]),theta_A[1]+theta_B[1]*exp(-theta_C[1]*z[j]));
}
else if (j<(nWT+nKO+1)){
log_lik[j]<- gumbel_log(y[j],theta_a[2]+theta_b[2]*exp(-theta_c[2]*x[j]),theta_A[2]+theta_B[2]*exp(-theta_C[2]*z[j]));
}
else{
log_lik[j]<- gumbel_log(y[j],theta_a[3]+theta_b[3]*exp(-theta_c[3]*x[j]),theta_A[3]+theta_B[3]*exp(-theta_C[3]*z[j]));

}

}
}


