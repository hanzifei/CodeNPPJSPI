#################################################################################################
#################################################################################################
#### Example 2: Normalized Power Prior in Borrowing Historical Control Data in Vaccine Trial ####
#### Binary Endpoints (response, non-response)                                               ####
#################################################################################################
#################################################################################################

## Calculate the (empirical) HPD from posterior sample
## Source: package TeachingDemos
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

###########################################################################################
#### Original Joint Power Prior Posterior Sampling Function for Bernoulli Data
#### Use likelihood as Product of Bernoulli
###########################################################################################
BerJPP1_MCMC <- function(Data.Cur = c(100, 50), Data.Hist = c(100, 50), 
                         prior = list(p.alpha = 1, p.beta = 1, 
                                      delta.alpha = 1, delta.beta = 1), 
                         MCMCmethod = 'IND', rw.logit.delta = 0.1, 
                         ind.delta.alpha= 1, ind.delta.beta= 1, nsample = 5000, 
                         control.mcmc = list(delta.ini = NULL, burnin = 0, thin = 1))
{
  
  y0 <- Data.Hist[2]; n0 <- Data.Hist[1]
  y1 <- Data.Cur[2]; n1 <- Data.Cur[1]
  
  
  LogPostBerDelta <- function(x){
    lbeta(x*y0 + y1 + prior$p.alpha, x*(n0-y0)+n1-y1+prior$p.beta)+
      (prior$delta.alpha-1)*log(x)+(prior$delta.beta-1)*log(1-x)
  }
  
  if(is.null(control.mcmc$delta.ini)) delta.ini = 0.5
  delta_cur <- delta.ini
  delta <- rep(delta.ini, nsample)
  counter <- 0
  niter <- nsample*control.mcmc$thin + control.mcmc$burnin
  
  for (i in 1:niter){
    ### Update delta with RW MH for Logit delta
    if(MCMCmethod == 'RW'){
      lgdelta_cur <- log(delta_cur/(1-delta_cur))
      lgdelta_prop <- rnorm(1, mean = lgdelta_cur, sd = sqrt(rw.logit.delta))
      delta_prop <- exp(lgdelta_prop)/(1+exp(lgdelta_prop))
      
      llik.prop <- LogPostBerDelta(delta_prop)
      llik.cur <- LogPostBerDelta(delta_cur)
      
      logr <- min(0, (llik.prop-llik.cur+log(delta_prop)+log(1-delta_prop)-log(delta_cur)-log(1-delta_cur)))
    }
    if(MCMCmethod == 'IND'){
      delta_prop <- rbeta(1, shape1 = ind.delta.alpha, shape2 = ind.delta.beta)
      llik.prop <- LogPostBerDelta(delta_prop)
      llik.cur <- LogPostBerDelta(delta_cur)
      
      logr <- min(0, (llik.prop-llik.cur+
                        dbeta(delta_cur, shape1 = ind.delta.alpha, shape2 = ind.delta.beta, log = TRUE) -
                        dbeta(delta_prop, shape1 = ind.delta.alpha, shape2 = ind.delta.beta, log = TRUE)))
    }
    
    if(runif(1) <= exp(logr)){
      delta_cur = delta_prop; counter = counter+1
    }
    
    if( i > control.mcmc$burnin & (i-control.mcmc$burnin)%%control.mcmc$thin==0) {
      delta[(i-control.mcmc$burnin)/control.mcmc$thin] <- delta_cur
    }
  }
  
  # Generate lambda conditional on delta 
  p <- rbeta(nsample, shape1 = delta*y0+y1+prior$p.alpha, 
             shape2 = delta*(n0-y0)+(n1-y1)+prior$p.beta)    
  meanp <- mean(p)
  D <- -2*(y1*log(p)+(n1-y1)*log(1-p))
  Dpbar <- -2*(y1*log(meanp)+(n1-y1)*log(1-meanp))
  DIC <- 2*mean(D)-Dpbar
  
  return(list(p = p, delta = delta, acceptrate = counter/niter, DIC = DIC))
}

###########################################################################################
#### Original Joint Power Prior Posterior Sampling Function for Bernoulli Data
#### Use binomial likelihood (including the (n choose x) constant)
###########################################################################################
BerJPP2_MCMC <- function(Data.Cur = c(100, 50), Data.Hist = c(100, 50), 
                         prior = list(p.alpha = 1, p.beta = 1, 
                                      delta.alpha = 1, delta.beta = 1), 
                         MCMCmethod = 'IND', rw.logit.delta = 0.1, 
                         ind.delta.alpha= 1, ind.delta.beta= 1, nsample = 5000, 
                         control.mcmc = list(delta.ini = NULL, burnin = 0, thin = 1))
{
  
  y0 <- Data.Hist[2]; n0 <- Data.Hist[1]
  y1 <- Data.Cur[2]; n1 <- Data.Cur[1]
  
  LogPostBerDelta <- function(x){
      x*lchoose(n0, y0)+  
      lbeta(x*y0 + y1 + prior$p.alpha, x*(n0-y0)+n1-y1+prior$p.beta)+
      (prior$delta.alpha-1)*log(x)+(prior$delta.beta-1)*log(1-x)
  }
  
  if(is.null(control.mcmc$delta.ini)) delta.ini = 0.5
  delta_cur <- delta.ini
  delta <- rep(delta.ini, nsample)
  counter <- 0
  niter <- nsample*control.mcmc$thin + control.mcmc$burnin
  
  for (i in 1:niter){
    ### Update delta with RW MH for Logit delta
    if(MCMCmethod == 'RW'){
      lgdelta_cur <- log(delta_cur/(1-delta_cur))
      lgdelta_prop <- rnorm(1, mean = lgdelta_cur, sd = sqrt(rw.logit.delta))
      delta_prop <- exp(lgdelta_prop)/(1+exp(lgdelta_prop))
      
      llik.prop <- LogPostBerDelta(delta_prop)
      llik.cur <- LogPostBerDelta(delta_cur)
      
      logr <- min(0, (llik.prop-llik.cur+log(delta_prop)+log(1-delta_prop)-log(delta_cur)-log(1-delta_cur)))
    }
    if(MCMCmethod == 'IND'){
      delta_prop <- rbeta(1, shape1 = ind.delta.alpha, shape2 = ind.delta.beta)
      llik.prop <- LogPostBerDelta(delta_prop)
      llik.cur <- LogPostBerDelta(delta_cur)
      
      logr <- min(0, (llik.prop-llik.cur+
                        dbeta(delta_cur, shape1 = ind.delta.alpha, shape2 = ind.delta.beta, log = TRUE) -
                        dbeta(delta_prop, shape1 = ind.delta.alpha, shape2 = ind.delta.beta, log = TRUE)))
    }
    
    if(runif(1) <= exp(logr)){
      delta_cur = delta_prop; counter = counter+1
    }
    
    if( i > control.mcmc$burnin & (i-control.mcmc$burnin)%%control.mcmc$thin==0) {
      delta[(i-control.mcmc$burnin)/control.mcmc$thin] <- delta_cur
    }
  }
  # Generate lambda conditional on delta 
  p <- rbeta(nsample, shape1 = delta*y0+y1+prior$p.alpha, 
             shape2 = delta*(n0-y0)+(n1-y1)+prior$p.beta)    
  meanp <- mean(p)
  D <- -2*(y1*log(p)+(n1-y1)*log(1-p))
  Dpbar <- -2*(y1*log(meanp)+(n1-y1)*log(1-meanp))
  DIC <- 2*mean(D)-Dpbar
  
  return(list(p = p, delta = delta, acceptrate = counter/niter, DIC = DIC))
}

###########################################################################################
#### Normalized Power Prior Posterior Sampling Function for Bernoulli Data
###########################################################################################
BerNPP_MCMC <- function(Data.Cur = c(100, 50), Data.Hist = c(100, 50), 
                        CompStat = list(y0 = NULL, n0 = NULL, y1 = NULL, n1 = NULL), 
                        prior = list(p.alpha = 1, p.beta = 1, 
                                     delta.alpha = 1, delta.beta = 1), 
                        MCMCmethod = 'IND', rw.logit.delta = 0.1, 
                        ind.delta.alpha = 1, ind.delta.beta = 1, nsample = 5000, 
                        control.mcmc = list(delta.ini = NULL, burnin = 0, thin = 1))
{
  if(missing(CompStat)){
    y0 <- Data.Hist[2]; n0 <- Data.Hist[1]
    y1 <- Data.Cur[2]; n1 <- Data.Cur[1]
  }else{
    y0 <- CompStat$y0
    n0 <- CompStat$n0
    y1 <- CompStat$y1
    n1 <- CompStat$n1
  }

  LogPostBerDelta <- function(x){
    lbeta(x*y0 + y1 + prior$p.alpha, x*(n0-y0)+n1-y1+prior$p.beta) -
      lbeta(x*y0 + prior$p.alpha, x*(n0-y0) + prior$p.beta) + 
      (prior$delta.alpha-1)*log(x)+(prior$delta.beta-1)*log(1-x)
  }
  
  if(is.null(control.mcmc$delta.ini)) delta.ini = 0.5
  delta_cur <- delta.ini
  delta <- rep(delta.ini, nsample)
  counter <- 0
  niter <- nsample*control.mcmc$thin + control.mcmc$burnin
  
  for (i in 1:niter){
    ### Update delta with RW MH for Logit delta
    if(MCMCmethod == 'RW'){
      lgdelta_cur <- log(delta_cur/(1-delta_cur))
      lgdelta_prop <- rnorm(1, mean = lgdelta_cur, sd = sqrt(rw.logit.delta))
      delta_prop <- exp(lgdelta_prop)/(1+exp(lgdelta_prop))
      
      llik.prop <- LogPostBerDelta(delta_prop)
      llik.cur <- LogPostBerDelta(delta_cur)
      
      logr <- min(0, (llik.prop-llik.cur+log(delta_prop)+log(1-delta_prop)-log(delta_cur)-log(1-delta_cur)))
    }
    if(MCMCmethod == 'IND'){
      delta_prop <- rbeta(1, shape1 = ind.delta.alpha, shape2 = ind.delta.beta)
      llik.prop <- LogPostBerDelta(delta_prop)
      llik.cur <- LogPostBerDelta(delta_cur)
      
      logr <- min(0, (llik.prop-llik.cur+
                        dbeta(delta_cur, shape1 = ind.delta.alpha, shape2 = ind.delta.beta, log = TRUE) -
                        dbeta(delta_prop, shape1 = ind.delta.alpha, shape2 = ind.delta.beta, log = TRUE)))}
    
    if(runif(1) <= exp(logr)){
      delta_cur = delta_prop; counter = counter+1
    }
    
    if( i > control.mcmc$burnin & (i-control.mcmc$burnin)%%control.mcmc$thin==0) {
      delta[(i-control.mcmc$burnin)/control.mcmc$thin] <- delta_cur
    }
  }
  # Generate lambda conditional on delta 
  p <- rbeta(nsample, shape1 = delta*y0+y1+prior$p.alpha, 
             shape2 = delta*(n0-y0)+(n1-y1)+prior$p.beta)    
  meanp <- mean(p)
  D <- -2*(y1*log(p)+(n1-y1)*log(1-p))
  Dpbar <- -2*(y1*log(meanp)+(n1-y1)*log(1-meanp))
  DIC <- 2*mean(D)-Dpbar
  return(list(p = p, delta = delta, acceptrate = counter/niter, DIC = DIC))
}


###########################################################################################
#### Original Joint Power Prior, Log-Posterior Densities on a Grid for Bernoulli Data
#### Use likelihood as Product of Bernoulli
###########################################################################################
LPDeltaBerJPP1 <- function(Data.Cur, Data.Hist, npoints = 1000, 
                           CompStat = list(y0 = NULL, n0 = NULL, y1 = NULL, n1 = NULL), 
                           prior = list(p.alpha = 1, p.beta = 1, 
                                        delta.alpha = 1, delta.beta = 1))
{
  if(missing(CompStat)){
    y0 <- Data.Hist[2]; n0 <- Data.Hist[1]
    y1 <- Data.Cur[2]; n1 <- Data.Cur[1]
  }else{
    y0 <- CompStat$y0
    n0 <- CompStat$n0
    y1 <- CompStat$y1
    n1 <- CompStat$n1
  }
  x <- seq(.Machine$double.eps, 1-.Machine$double.eps, length = npoints)
  propDen <- lbeta(x*y0 + y1 + prior$p.alpha, x*(n0-y0)+n1-y1+prior$p.beta)+
    (prior$delta.alpha-1)*log(x)+(prior$delta.beta-1)*log(1-x)
  return(data.frame(x = x, logden = propDen))
}

###########################################################################################
#### Original Joint Power Prior, Log-Posterior Densities on a Grid for Bernoulli Data
#### Use binomial likelihood (including the (n choose x) constant)
###########################################################################################
LPDeltaBerJPP2 <- function(Data.Cur, Data.Hist, npoints = 1000, 
                           CompStat = list(y0 = NULL, n0 = NULL, y1 = NULL, n1 = NULL), 
                           prior = list(p.alpha = 1, p.beta = 1, 
                                        delta.alpha = 1, delta.beta = 1))
{
  if(missing(CompStat)){
    y0 <- Data.Hist[2]; n0 <- Data.Hist[1]
    y1 <- Data.Cur[2]; n1 <- Data.Cur[1]
  }else{
    y0 <- CompStat$y0
    n0 <- CompStat$n0
    y1 <- CompStat$y1
    n1 <- CompStat$n1
  }
  x <- seq(.Machine$double.eps, 1-.Machine$double.eps, length = npoints)
  propDen <- lbeta(x*y0 + y1 + prior$p.alpha, x*(n0-y0)+n1-y1+prior$p.beta)+
    x*lchoose(n0, y0)+
    (prior$delta.alpha-1)*log(x)+(prior$delta.beta-1)*log(1-x)
  return(data.frame(x = x, logden = propDen))
}


###########################################################################################
#### Normalized Power Prior, Log-Posterior Densities on a Grid for Bernoulli Data
###########################################################################################
LPDeltaBerNPP <- function(Data.Cur, Data.Hist, npoints = 1000, 
                          CompStat = list(y0 = NULL, n0 = NULL, y1 = NULL, n1 = NULL), 
                          prior = list(p.alpha = 1, p.beta = 1, 
                                       delta.alpha = 1, delta.beta = 1))
{
  if(missing(CompStat)){
    y0 <- Data.Hist[2]; n0 <- Data.Hist[1]
    y1 <- Data.Cur[2]; n1 <- Data.Cur[1]
  }else{
    y0 <- CompStat$y0
    n0 <- CompStat$n0
    y1 <- CompStat$y1
    n1 <- CompStat$n1
  }
  x <- seq(.Machine$double.eps, 1-.Machine$double.eps, length = npoints)
  propDen <- lbeta(x*y0 + y1 + prior$p.alpha, x*(n0-y0)+n1-y1+prior$p.beta)-
    lbeta(x*y0 + prior$p.alpha, x*(n0-y0) + prior$p.beta) + 
    (prior$delta.alpha-1)*log(x)+(prior$delta.beta-1)*log(1-x)
  return(data.frame(x = x, logden = propDen))
}


###########################################################################################
#### Start the model fitting ####
###########################################################################################

### Historical Control 1: y01 = 417; N = 576
### Historical Control 2: y02 = 90; N = 111
### Historical Control 3: y03 = 49; N = 62
### Historical Control 4: y04 = 376; N = 487
### Total Historical Control: y = 932; N = 1236
### Current Control: y1 = 426; N = 592
### Current Test: y = 415; N = 558

#### Using Jeffrey's Prior, to get the Posterior of the Treatment Group
set.seed(4321)
Trtgrp <- rbeta(50000, shape1 = 1/2+415, shape2 = 1/2+558-415)


###########################################################################################
#### If no borrowing, with Jeffrey's prior for p
###########################################################################################
set.seed(1234567)
CControl <- rbeta(50000, shape1 = 1/2+426, shape2 = 1/2+592-426)

#### Estimates
round(100*mean(CControl),2)  #71.92%
plot(density(CControl, 0.003)) 

#### 95% HPD for treatment vs control: p_t- p_c
plot(Trtgrp - CControl, type = 'l')
plot(density(Trtgrp - CControl, 0.003))
round(100*quantile(Trtgrp - CControl, c(0.025, 0.975)),2)  # -2.72, 7.49 for equal-tail Credible Interval
round(100*emp.hpd(Trtgrp - CControl, conf = 0.95), 2)      # -2.61, 7.58 for 95% HPD 



###########################################################################################
#### If we use joint power prior, with product of Bernoulli likelihood
###########################################################################################
set.seed(12345)
ControlJPP1 <- BerJPP1_MCMC(Data.Cur = c(592, 426), Data.Hist = c(1236, 932), 
                            prior = list(p.alpha = 0.5, p.beta = 0.5, delta.alpha = 1, delta.beta = 1), 
                            MCMCmethod = 'RW', rw.logit.delta = 1, 
                            ind.delta.alpha= 1, ind.delta.beta= 1, nsample = 50000,
                            control.mcmc = list(delta.ini = NULL, burnin = 2000, thin = 5))

#### Model Fitting Checking; converged well
plot(ControlJPP1$delta, type = 'l')
plot(ControlJPP1$p, type = 'l')
acf(ControlJPP1$p)
acf(ControlJPP1$delta)

#### Estimates
round(100*mean(ControlJPP1$p), 2)         # 71.93%
round(mean(ControlJPP1$delta), 3)         # 0.001
plot(density(ControlJPP1$p, 0.003)) 
plot(density(ControlJPP1$delta, 0.0001))

#### 95% Credible Intervals for treatment vs control: p_t- p_c
plot(Trtgrp - ControlJPP1$p, type = 'l')
plot(density(Trtgrp - ControlJPP1$p, 0.003))
100*mean(Trtgrp - ControlJPP1$p)  
round(100*quantile(Trtgrp - ControlJPP1$p, c(0.025, 0.975)), 2)  # -2.70, 7.51 for equal-tail Credible Interval
round(100*emp.hpd(Trtgrp - ControlJPP1$p, conf = 0.95), 2)       # -2.89, 7.31 for 95% HPD  

#### Posterior Mode 
JPP1postLik <- LPDeltaBerJPP1(Data.Cur = c(592, 426), Data.Hist = c(1236, 932), npoints = 5000, 
                              prior = list(p.alpha = 0.5, p.beta = 0.5, delta.alpha = 1, delta.beta = 1))
round(JPP1postLik[which.max(JPP1postLik$logden), ][1], 3)        # 0




###########################################################################################
#### If we use the original power prior, using binomial likelihood 
#### Including the (N0 choose k0)^delta in the likelihood 
###########################################################################################
set.seed(123456)
ControlJPP2 <- BerJPP2_MCMC(Data.Cur = c(592, 426), Data.Hist = c(1236, 932), 
                      prior = list(p.alpha = 0.5, p.beta = 0.5, delta.alpha = 1, delta.beta = 1), 
                      MCMCmethod = 'RW', rw.logit.delta = 1, nsample = 50000,
                      control.mcmc = list(delta.ini = NULL, burnin = 2000, thin = 5))

#### Model Fitting Checking; converged well
plot(ControlJPP2$delta, type = 'l')
plot(ControlJPP2$p, type = 'l')
acf(ControlJPP2$p)
acf(ControlJPP2$delta)

#### Estimates
round(100*mean(ControlJPP2$p), 2)         # 72.68%
round(mean(ControlJPP2$delta), 3)         # 0.1656
plot(density(ControlJPP2$p, 0.003)) 
plot(density(ControlJPP2$delta, 0.01))

#### 95% Credible Interval for treatment vs control: p_t- p_c
plot(Trtgrp - ControlJPP2$p, type = 'l')
plot(density(Trtgrp - ControlJPP2$p, 0.003))
mean(Trtgrp - ControlJPP2$p)                   #0.165
round(100*quantile(Trtgrp - ControlJPP2$p, c(0.025, 0.975)), 2)  # -3.23, 6.62 for equal-tail Credible Interval
round(100*emp.hpd(Trtgrp - ControlJPP2$p, conf = 0.95), 2)       # -3.26, 6.59 for 95% HPD 

#### Posterior Mode
JPP2postLik <- LPDeltaBerJPP2(Data.Cur = c(592, 426), Data.Hist = c(1236, 932), npoints = 5000, 
                              prior = list(p.alpha = 0.5, p.beta = 0.5, delta.alpha = 1, delta.beta = 1))
round(JPP2postLik[which.max(JPP2postLik$logden), ][1], 3)  # 0




###########################################################################################
#### Using Normalized Power Prior to borrow historical control
###########################################################################################
set.seed(1234)
ControlNPP <- BerNPP_MCMC(Data.Cur = c(592, 426), Data.Hist = c(1236, 932), 
                          prior = list(p.alpha = 0.5, p.beta = 0.5, delta.alpha = 1, delta.beta = 1), 
                          MCMCmethod = 'RW', rw.logit.delta = 1, nsample = 50000,
                          control.mcmc = list(burnin = 2000, thin = 5))

#### Model Fitting Checking; converged well
plot(ControlNPP$delta, type = 'l')
plot(ControlNPP$p, type = 'l')
acf(ControlNPP$p)
acf(ControlNPP$delta)

#### Estimates
round(100*mean(ControlNPP$p), 2)              # 73.50%
round(mean(ControlNPP$delta), 3)              # 0.4817
plot(density(ControlNPP$p, 0.003))
plot(density(ControlNPP$delta, 0.1))

#### 95% Credible Intervals for treatment vs control: p_t- p_c
plot(Trtgrp - ControlNPP$p, type = 'l')
plot(density(Trtgrp - ControlNPP$p, 0.003))
mean(Trtgrp - ControlNPP$p)                   #0.0083
round(100*quantile(Trtgrp - ControlNPP$p, c(0.025, 0.975)), 2)  # -3.76, 5.54 for Equal Tail
round(100*emp.hpd(Trtgrp - ControlNPP$p, conf = 0.95), 2)       #-3.76 5.54 for HPD 

#### Posterior Mode
NPPpostLik <- LPDeltaBerNPP(Data.Cur = c(592, 426), Data.Hist = c(1236, 932), npoints = 5000, 
                            prior = list(p.alpha = 0.5, p.beta = 0.5, delta.alpha = 1, delta.beta = 1))
NPPpostLik[which.max(NPPpostLik$logden), ]
round(NPPpostLik[which.max(NPPpostLik$logden), ][1], 3)  # 0.181
