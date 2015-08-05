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

real<lower=-1,upper=1> alpha;
real<lower=0,upper=2> beta;
real<lower=0,upper=0.5> c;
real<lower=0,upper=1> A;
real<lower=0,upper=1> B;
real<lower=0,upper=0.5> C; 

}

model {
vector[n] Ymu;
vector[n] Ysd;
			
		Ymu<-alpha+beta*exp(-c*x);
		Ysd<-A+B*exp(-C*z);							
		y~gumbel(Ymu,Ysd);
				
}


generated quantities{
vector[n] log_lik;


for(j in 1:n){
log_lik[j]<- gumbel_log(y[j],alpha+beta*exp(-c*x[j]),A+B*exp(-C*z[j]));
}
}

