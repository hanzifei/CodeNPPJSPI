functions{
  real interpolateC(real deltaCur, int numdelta, vector deltaKnot, vector logCKnot){
  real logCest;
  for(id in 1:(numdelta-1)){
    if(deltaCur >= deltaKnot[id] && deltaCur < deltaKnot[id+1]){
      logCest = logCKnot[id]+ (deltaCur-deltaKnot[id])*(logCKnot[id+1]-logCKnot[id])/(deltaKnot[id+1]-deltaKnot[id]);
    }  // Interpolation function, given a sequence of logCKnot
  }
  return logCest;
}
}


data {
int<lower = 0> K0;              // number of stress levels in D0
int<lower = 0> K;               // number of stress levels in D
vector<lower = 0>[K0] N0;       // N0 = (N10, N20, ..., NK0) - units at risk for each stress level in D0, constant stress
vector<lower = 0>[K] N;         // N = (N1, N2, ..., NK) - number of units at risk for each stress level in D, step stress
vector<lower = 0>[K0] tau0;     // tau0 = (tau10, tau20, ..., tauK0) - endtime for each stress levels, constant stress
vector<lower = 0>[K] tau;       // tau = (tau1, tau2, ..., tauK) - endtime for each stress levels, step stress
vector<lower = 0>[K] taum1;     // taum1 = (0, tau1, ..., tauK-1) - endtime for i-1 stress level, with first element 0
vector<lower = 0>[K0] x0;       // x0 = (x10, x20, ..., xK0) - the value of the stress levels in D0, constant stress
vector<lower = 0>[K] x;         // x = (x1, x2, ..., xK) - the value of the stress levels in D, step stress
vector<lower = 0>[K0] n0;       // n = (n10, n20, ..., nK0) - observed number of subjects failures each level, constant stress
vector<lower = 0>[K] n;         // n = (n1, n2, ..., nK) - observed number of subjects failures each level, step stress
vector<lower = 0>[K0] ys0;      // ys0 = (ys01, ys02, ..., ys0K) - sum of survival time at each levels constant stress
vector<lower = 0>[K] ys;        // ys = (ys1, ys2, ..., ysK) - sum of survival time at each levels step stress
real<lower = 0> a_delta_prior;  // Prior of the first parameter in a Beta distribution
real<lower = 0> b_delta_prior;  // Prior of the second parameter in a Beta distribution
real<lower = 0> alphapriorsd; 
real<lower = 0> betapriorsd;
int<lower = 0> numdelta;        // number of levels for delta
vector<lower = 0>[numdelta] deltaKnot;  
vector[numdelta] logCKnot;
}



parameters {
real<lower = 1e-10, upper = 1> delta; //delta
real alpha;              // alpha
real<upper = 0> beta;    // beta
}

model{
delta ~ beta(a_delta_prior, b_delta_prior);
alpha ~ normal(0, alphapriorsd); 
//lnegbeta ~ normal(0, 1000); // prior
beta ~normal(0, betapriorsd);
{
vector[K0] acc0;
vector[K] acc;
real C; 
  for(j in 1:K0){
    acc0[j] = -((N0[j]-n0[j])*tau0[j]+ys0[j])/exp(alpha+beta*x0[j])-n0[j]*(alpha+beta*x0[j]);
  }
  for(i in 1:K){
    acc[i] = -n[i]*(alpha+beta*x[i])-
    ((ys[i]-n[i]*taum1[i])+(N[i]-n[i])*(tau[i]-taum1[i]))/exp(alpha+beta*x[i]);
  }
  C =  interpolateC(delta, numdelta, deltaKnot, logCKnot);  //delta*llikMLE0-(log(delta)+log(detH))/2;
  target += delta*sum(acc0)+sum(acc)-C;       // likelihood  
 }
}
