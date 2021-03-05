data {
int<lower = 0> K;               // number of stress levels
int<lower = 0> numdelta;        // number of levels for delta
vector<lower = 0>[K] N;         // N = (N1, N2, ..., NK) - units at risk for each stress levels
vector<lower = 0>[K] tau;       // tau = (tau1, tau2, ..., tauK) - endtime for each stress levels
vector<lower = 0>[K] x;         // x = (x1, x2, ..., xK) - the value of the stress levels
vector<lower = 0>[K] n;         // n = (n1, n2, ..., nK) - observed number of subjects failures
vector<lower = 0>[K] ys;        // ys = (ys1, ys2, ..., ysK) - average survival time at each levels
vector<lower = 0, upper = 1>[numdelta] delta;     //vector delta: power parameter in power prior
real<lower = 0> alphapriorsd; 
real<lower = 0> betapriorsd; 
}

parameters {
vector[numdelta] alpha;              //vector alpha
vector<upper = 0>[numdelta] beta;    //vector beta
}

model{
alpha ~ normal(0, alphapriorsd); 
beta ~ normal(0, betapriorsd); // prior
{
vector[K] acc;
for (i in 1:numdelta) {
  for(j in 1:K){
    acc[j] = -( ((N[j]-n[j])*tau[j]+ys[j])/exp(alpha[i]+beta[i]*x[j]) + n[j]*(alpha[i]+beta[i]*x[j]) );
  }
  target += delta[i]*sum(acc);       // likelihood  
 }
 }
}
