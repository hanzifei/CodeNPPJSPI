library(rstan) 
setwd("~/Desktop/NppALT/ALTRegExp")
## Directory of the Stan Source files 
rstan_options(auto_write = TRUE)
options(timeout= 1e20)  # To avoid loading error

#### If run the program first time, need to run the following
# StanALTFix=rstan::stan_model(file = "StanALTFix.stan", verbose = FALSE)
# save('StanALTFix', file = "StanALTFix.RData")
# 
# StanPathLogC=rstan::stan_model(file = "StanPathLogC.stan", verbose = FALSE)
# save('StanPathLogC', file = "StanPathLogC.RData")
# 
# StanALTRandPath=rstan::stan_model(file = "StanALTRandPath.stan", verbose = FALSE)
# save('StanALTRandPath', file = "StanALTRandPath.RData")


#### If we save the Stan into RData, we don't need to re-compile (if no version update)
### So simply loading may works 

### If loading these files has error on your computer, please re-compile with stan_model
#load(file = "StanALTFix.RData")
#load(file = "StanPathLogC.RData")
#load(file = "StanALTRandPath.RData")


#### DATA0: K0 levels constant stress test
#### DATA1: K levels step-stress test
#### delta: The level of discounting 
#### iter: number of MCMC iteration
#### cores: Number of cores used in parallel
#### priorsd: Standard deviation of the prior for \beta, which follows N(0, priorsd^2)
#### D0id, D1id: The Data id used in simulation studies; 
####             It is 1 for this specific example since we only have 1 data
ExpRegFix_MCMC <- function(DATA0, DATA1, delta, iter = 5000, cores = 1, 
                           priorsd = c(50, 50), D0id = 1, D1id = 1){
  D1Sum = DATA1$DataSum[DATA1$DataSum$id == D1id, ]
  D0Sum = DATA0$DataSum[DATA0$DataSum$id == D0id, ]
  
  x = DATA1$x; n = as.numeric(D1Sum[1,2:(DATA1$K+1)]); tau = DATA1$tau
  x0 = DATA0$x; n0 = as.numeric(D0Sum[1,2:(DATA0$K+1)]); tau0 = DATA0$tau
  K = DATA1$K; K0 = DATA0$K
  ys = D1Sum$ys; ys0 = D0Sum$ys
  
  StanFix_dat = list(K0 = K0, K = K, N0 = DATA0$N, 
                     N = c(DATA1$N, DATA1$N-D1Sum$n1[1], DATA1$N-D1Sum$n1[1]-D1Sum$n2[1]), 
                     tau0 = tau0, tau = tau, taum1 = c(0, tau[1:(DATA1$K-1)]), 
                     x0 = x0, x =x, n0 = n0, n = n, ys0 = ys0, ys = ys, delta = delta, 
                     alphapriorsd = priorsd[1], betapriorsd = priorsd[2])
  
  FixStan = suppressWarnings(sampling(object = StanALTFix, data = StanFix_dat, iter = iter, chains = 1, 
                                      refresh = 0, cores = cores))
  
  MCMCalpha = rstan::extract(FixStan)$alpha; MCMCbeta = rstan::extract(FixStan)$beta
  
  llikFix <- function(par){
    alpha = par[1]; beta = par[2]
    muexp = exp(alpha+beta*x)                    # Vector of K elements 
    #  const = lgamma(N+1)-lgamma(N-sum(n)+1)    # Constant 
    Ni = DATA1$N-cumsum(c(n))                    # Vector of K elements; N(1), N(2), N(3),..., N(K)
    TAU = c(0, tau)                              # Vector of K+1 elements; 0, tau(1), ..., tau(K)
    DiffTAU = diff(TAU)                          # Vector of K elements; tau(1), tau(2)-tau(1), ..., tau(K)-tau(K-1)
    tauminus1 = TAU[1:K]                         # Vector of K elements; 0, tau(1), ..., tau(K-1)
    A = sum(n*(alpha+beta*x))                    # sum of vector of K elements; sum(n(i)*mu(i))
    #  B = sum((y-rep(tauminus1, times = n))/(rep(mu, times = n)))
    B = sum((ys-n*tauminus1)/muexp)
    # sum of vector of sum(n(i)) elements; order matched  
    C = sum(DiffTAU*Ni/muexp)                     # sum of vector of K elements; N(i)*(tau(i)-tau(i-1))/mu(i)
    ans = -A-B-C
    return(ans)
  }
  meanalpha = mean(MCMCalpha)
  meanbeta = mean(MCMCbeta)
  estmat = cbind(MCMCalpha, MCMCbeta)
  D = -2*apply(estmat, 1, llikFix)
  Destbar = -2*llikFix(c(meanalpha, meanbeta))
  DIC = 2*mean(D)-Destbar
  return(list(alpha = MCMCalpha, beta = MCMCbeta, DIC = DIC))
}

#### NPP Using Path Sampling to calculate C (The "path sampling" is the Algorithm in Appendix B2 since it is similar, we use an aliased name)
#### DATA0: K0 levels constant stress test
#### DATA1: K levels step-stress test
#### iter: number of MCMC iteration
#### cores: Number of cores used in parallel
#### delta_prior: hyperparameters for beta prior of \delta
#### priorsd: Standard deviation of the prior for \beta, which follows N(0, priorsd^2)
#### deltaKnot: the knots when we calculate the C(\delta)
#### numbdelta: number of knots we used; it is length(deltaKnot)
#### iterC: Number of MCMC iterations we use to calculate C(deltaKnot)
#### coresC: Number of cores used to calculate C(deltaKnot) in parallel
#### D0id, D1id: The Data id used in simulation studies; 
####             It is 1 for this specific example since we only have 1 data

ExpRegRandPath_MCMC <- function(DATA0, DATA1, iter = 5000, cores = 1, delta_prior = c(0.5, 0.5),
                                priorsd = c(10, 10), deltaKnot = (0:200/200)^2, numdelta = 201, 
                                iterC = 2000, coresC = 8, D0id = 1, D1id = 1){
  D0Sum = DATA0$DataSum[DATA0$DataSum$id == D0id, ]
  D1Sum = DATA1$DataSum[DATA1$DataSum$id == D1id, ]
  
  x = DATA1$x; n = as.numeric(D1Sum[1,2:(DATA1$K+1)]); tau = DATA1$tau
  x0 = DATA0$x; n0 = as.numeric(D0Sum[1,2:(DATA0$K+1)]); tau0 = DATA0$tau
  K = DATA1$K; K0 = DATA0$K
  ys = D1Sum$ys; ys0 = D0Sum$ys; N0 = DATA0$N
  
  CT_dat <- list(K = K0, numdelta = numdelta, N = N0, tau = tau0, x = x0, 
                 n = n0, ys = ys0, delta = deltaKnot, alphapriorsd = priorsd[1], betapriorsd = priorsd[2])
  
  STstan = suppressWarnings(sampling(object = StanPathLogC, data = CT_dat, iter = iterC, 
                                     chains = 1, refresh = 0, cores = coresC))
  
  alphaMat = rstan::extract(STstan)$alpha; betaMat = rstan::extract(STstan)$beta
  mu1 = alphaMat+betaMat*x0[1]  
  mu2 = alphaMat+betaMat*x0[2]                   
  A1 = ((N0[1]-n0[1])*tau0[1]+ys0[1])/exp(mu1)                     
  A2 = ((N0[2]-n0[2])*tau0[2]+ys0[2])/exp(mu2)                    
  A = A1+A2
  B = n0[1]*mu1+n0[2]*mu2     
  llikMat = -(A+B)
  loglikmean <- apply(llikMat, 2, mean)
  dvec <- c(deltaKnot[1], diff(deltaKnot))
  lgC_knot <- cumsum(dvec*loglikmean)
  
  StanPath_dat <- list(K0 = K0, K = K, N0 = N0, 
                       N = c(DATA1$N, DATA1$N-DATA1$DataSum$n1[1], DATA1$N-DATA1$DataSum$n1[1]-DATA1$DataSum$n2[1]), 
                       tau0 = tau0, tau = tau, taum1 = c(0, tau[1:(DATA1$K-1)]), x0 = x0, 
                       x = x, n0 = n0, n = n, ys0 = ys0, ys = ys, 
                       a_delta_prior = delta_prior[1], b_delta_prior = delta_prior[2], 
                       alphapriorsd = priorsd[1], betapriorsd = priorsd[2], 
                       numdelta = numdelta, deltaKnot = deltaKnot, logCKnot = lgC_knot)
  
  RandomStan = sampling(object = StanALTRandPath, data = StanPath_dat, 
                        iter = iter, chains = 1, refresh = 0, cores = cores)
  
  ExtractMCMC = rstan::extract(RandomStan)
  MCMCalpha = ExtractMCMC$alpha
  MCMCbeta = ExtractMCMC$beta
  MCMCdelta = ExtractMCMC$delta
  
  llST <- function(par){
    alpha = par[1]; beta = par[2]
    muexp = exp(alpha+beta*x)                    
    Ni = DATA1$N-cumsum(c(n))                    
    TAU = c(0, tau)                              
    DiffTAU = diff(TAU)                          
    tauminus1 = TAU[1:K]                         
    A = sum(n*(alpha+beta*x))               
    B = sum((ys-n*tauminus1)/muexp)
    C = sum(DiffTAU*Ni/muexp)                   
    ans = -A-B-C
    return(ans)
  }
  meanalpha = mean(MCMCalpha)
  meanbeta = mean(MCMCbeta)
  estmat = cbind(MCMCalpha, MCMCbeta)
  D = -2*apply(estmat, 1, llST)
  Destbar = -2*llST(c(meanalpha, meanbeta))
  DIC = 2*mean(D)-Destbar
  return(list(alpha = MCMCalpha, beta = MCMCbeta, delta = MCMCdelta, DIC = DIC))
}

### Empirical 95% CI based on posterior HPD
emp.hpd <- function(x, conf=0.95){
  conf <- min(conf, 1-conf)
  n <- length(x)
  nn <- round( n*conf )
  x <- sort(x)
  xx <- x[ (n-nn+1):n ] - x[1:nn]
  m <- min(xx)
  nnn <- which(xx==m)[1]
  return( c( x[ nnn ], x[ n-nn+nnn ] ) )
}

########################################################################################################
########################## ABOVE ARE FUNCTIONS MODIFIED FOR THE DATA ###################################
########################################################################################################
########################## Data Analysis Starts Below ##################################################
########################################################################################################

### Time to failures 
### Historical Data
Time01 = c(0.611, 2.337, 5.461, 5.726, 7.138, 9.832, 10.021, 12.136, 12.375, 13.996)
Time02 = c(0.081, 0.828, 1.116, 1.839, 5.043, 6.717, 7.258, 8.472)

## Current Data
Time11 = c(1.515, 2.225, 4.629, 4.654, 6.349, 8.003, 8.262, 10.416, 11.381, 12.433, 14.755)
Time12 = c(15.164, 15.355, 15.953, 16.654, 16.735, 18.796, 19.248, 19.295)
Time13 = c(20.110, 20.318, 21.228, 21.543, 22.227, 24.541)


#### Historical Dataset. In order to follow format to adopt functions used in simulation
#### 2-Level Constant Stress Test
D0solar = list()
D0solar$Data = data.frame(id = 1, n1 = 10, n2 = 8, c1 = 10, c2 = 2, 
                          level = c(rep(1, length(Time01)), rep(2, length(Time02))), 
                          y = c(Time01, Time02))
D0solar$DataSum = data.frame(id = 1, n1 = 10, n2 = 8, c1 = 10, c2 = 2, 
                             level = c(1, 2), ys = c(sum(Time01), sum(Time02)))
D0solar$N = c(20, 10); D0solar$K = 2; D0solar$tau = c(15, 10); D0solar$x = c(0.1, 0.6)

#### Current Dataset. In order to follow format to adopt functions used in simulation
D1solar = list()
D1solar$Data = data.frame(id = 1, n1 = 11, n2 = 8, n3 = 6, c = 5, 
                          level = c(rep(1, length(Time11)), rep(2, length(Time12)), rep(3, length(Time13))), 
                          y = c(Time11, Time12, Time13))
D1solar$DataSum = data.frame(id = 1, n1 = 11, n2 = 8, n3 = 6, c = 5, 
                             level = c(1, 2, 3), 
                             ys = c(sum(Time11), sum(Time12), sum(Time13)) )
D1solar$N = 30; D1solar$K = 3; D1solar$tau = c(15, 20, 25); D1solar$x = c(0.1, 0.5, 0.9)


#################################################################################################################
#######################################  DATA FITTING ###########################################################
#################################################################################################################
#### Full Borrowing 
set.seed(123)
SolarFB = ExpRegFix_MCMC(DATA0 = D0solar, DATA1 = D1solar, delta = 1, iter = 10000, 
                         cores = 1, priorsd = c(50, 50), D0id = 1, D1id = 1)
plot(SolarFB$alpha, type = 'l'); round(c(mean(SolarFB$alpha), emp.hpd(SolarFB$alpha)), 3)
#  3.530 3.040 4.002
plot(SolarFB$beta, type = 'l'); round(c(mean(SolarFB$beta), emp.hpd(SolarFB$beta)), 3)
# -2.228 -3.135 -1.135

#### No Borrowing 
set.seed(1234)
SolarNB = ExpRegFix_MCMC(DATA0 = D0solar, DATA1 = D1solar, delta = 0, iter = 10000, 
                         cores = 1, priorsd = c(50, 50), D0id = 1, D1id = 1)
plot(SolarNB$alpha, type = 'l'); round(c(mean(SolarNB$alpha), emp.hpd(SolarNB$alpha)), 3)
# 3.692 3.091 4.332
plot(SolarNB$beta, type = 'l'); round(c(mean(SolarNB$beta), emp.hpd(SolarNB$beta)), 3)
# -2.280 -3.423 -1.092


#### NPP with Path Sampling 
set.seed(232)
PathNPP1 = ExpRegRandPath_MCMC(DATA0 = D0solar, DATA1 = D1solar, iter = 10000, cores = 4, delta_prior = c(1, 1),
                               deltaKnot = ((0:200)/200)^2, numdelta = 201, 
                               iterC = 2000, coresC = 8, priorsd = c(50, 50), D0id = 1, D1id = 1)

round(c(mean(PathNPP1$alpha), emp.hpd(PathNPP1$alpha)),3) 
# 3.571 3.061 4.112
round(c(mean(PathNPP1$beta), emp.hpd(PathNPP1$beta)), 3)
# -2.230 -3.231 -1.176
round(c(mean(PathNPP1$delta), emp.hpd(PathNPP1$delta)), 3)
# 0.593 0.141 0.999

### Some rough plots
par(mgp=c(1.9,0.8,0), mai = c(0.4, 0.4, 0.1, 0.1), mfrow = c(1,3))
plot(PathNPP1$alpha, type = 'p', ylab = expression(alpha), cex = 0.1)
plot(PathNPP1$beta,  ylab = expression(beta), type = 'p', cex = 0.1)  
plot(PathNPP1$delta, ylab = expression(delta), type = 'p', cex = 0.1)


acf(PathNPP1$alpha)
acf(PathNPP1$beta)  
acf(PathNPP1$delta)
