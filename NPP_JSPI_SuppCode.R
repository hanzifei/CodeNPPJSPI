########################################################################################################
########################################################################################################
#### Example 1: Normalized Power Prior in Borrowing Historical pH Data in Water Quality Assessment  ####
#### Normal Endpoints (Four Sites; using Four Priors)                                             ######
#### Assume Joint Prior of mu and sigmasq is (1/sigmasq)^a                                   ###########
########################################################################################################
########################################################################################################

library(ggplot2)
library(ggpubr)

#### Function to sample from location-scale t distribution
rt_ls <- function(n, df, location, scale) rt(n,df)*scale + location
#### Function to sample from Inverse Gamma
r_igamma <- function(n, shape, rate = 1, scale = 1/rate){
  if(missing(rate) && !missing(scale)) rate <- 1/scale
  1/rgamma(n, shape, rate)
}
#### Function to sample from Inverse Chisq
rinvchisq <- function(n, df, scale=1/df)
{
  df <- rep(df, len=n); scale <- rep(scale, len=n)
  if(any(df <= 0)) stop("The df parameter must be positive.")
  if(any(scale <= 0)) stop("The scale parameter must be positive.")
  z <- rchisq(n, df=df)
  z[which(z == 0)] <- 1e-100
  x <- (df*scale) / z
  return(x)
}


#######################################################################################################
#### Original Joint Power Prior Posterior Sampling Function for Normal Data                     #######
#### Use likelihood as Product of Densities of Sufficient Statistics (Normal and Gamma Density) #######
#### Considering the constant 2pi                                                               #######
#######################################################################################################

NormalOPP_MCMC_suff <- function(Data.Cur, Data.Hist, 
                                CompStat = list(mean0 = NULL, n0 = NULL, var0 = NULL, mean1 = NULL, n1 = NULL, var1 = NULL), 
                                prior = list(joint.a = 1.5, delta.alpha = 1, delta.beta = 1), 
                                MCMCmethod = 'IND', rw.logit.delta = 0.1, 
                                ind.delta.alpha= 1, ind.delta.beta= 1, nsample = 5000, 
                                control.mcmc = list(delta.ini = NULL, burnin = 0, thin = 1))
{
  if(missing(CompStat)){
    mean0 <- mean(Data.Hist)
    n0 <- length(Data.Hist)
    var0 <- var(Data.Hist)*(n0-1)/n0
    mean1 <- mean(Data.Cur)
    n1 <- length(Data.Cur)
    var1 <- var(Data.Cur)*(n1-1)/n1
  }else{
    mean0 <- CompStat$mean0
    n0 <- CompStat$n0
    var0 <- CompStat$var0
    mean1 <- CompStat$mean1
    n1 <- CompStat$n1
    var1 <- CompStat$var1
  }
  LogPostNDelta <- function(x){
    K <- (x*n0*n1*(mean1-mean0)^2)/(2*(x*n0+n1))+
      (x*n0*var0+n1*var1)/2
    s0 <- var0*n0/(n0-1)
    logC <- ((n0-1)/2)*log((n0-1)/2) + ((n0-3)/2)*log(s0) + 
      0.5*(log(n0)-log(2*pi)) - lgamma((n0-1)/2)
    propDen <- (prior$delta.alpha-1)*log(x)+(prior$delta.beta-1)*log(1-x)-
      log(x*n0+n1)/2+lgamma((x*n0+n1-3)/2 +prior$joint.a)-
      ((x*n0+n1-3)/2 + prior$joint.a)*log(K)+x*logC
    return(propDen)
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
      llik.prop <- LogPostNDelta(delta_prop)
      llik.cur <- LogPostNDelta(delta_cur)
      logr <- min(0, (llik.prop-llik.cur+log(delta_prop)+log(1-delta_prop)-log(delta_cur)-log(1-delta_cur)))
    }
    if(MCMCmethod == 'IND'){
      delta_prop <- rbeta(1, shape1 = ind.delta.alpha, shape2 = ind.delta.beta)
      llik.prop <- LogPostNDelta(delta_prop)
      llik.cur <- LogPostNDelta(delta_cur)
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
  #### Generate mu and sigmasq conditional on delta 
  K = (delta*n0*n1*(mean0-mean1)^2/(delta*n0+n1) + delta*n0*var0 + n1*var1 )/2
  mu = rt_ls(nsample, df= delta*n0+n1+2*prior$joint.a-3, location= (delta*n0*mean0+n1*mean1)/(delta*n0 + n1), 
             scale= sqrt(2*K/( (delta*n0+n1+2*prior$joint.a-3)*(delta*n0+n1) )) )
  sigmasq= r_igamma(n = nsample, shape = (delta*n0+n1+2*prior$joint.a-3)/2, rate = K)
  return(list(mu = mu, sigmasq = sigmasq, delta = delta, acceptrate = counter/niter))
}


#### Mode Search 
ModeDeltaNormalOPP_suff <- function(Data.Cur, Data.Hist,
                               CompStat = list(n0 = NULL, mean0 = NULL, sd0 = NULL, n1 = NULL, mean1 = NULL, sd1 = NULL),
                               npoints = 1000,
                               prior = list(joint.a = 1.5, delta.alpha = 1, delta.beta = 1))
{
  if(missing(CompStat)){
    mean0 <- mean(Data.Hist)
    n0 <- length(Data.Hist)
    var0 <- var(Data.Hist)*(n0-1)/n0
    mean1 <- mean(Data.Cur)
    n1 <- length(Data.Cur)
    var1 <- var(Data.Cur)*(n1-1)/n1
  }else{
    mean0 <- CompStat$mean0
    n0 <- CompStat$n0
    var0 <- (CompStat$sd0^2)*(1-1/n0)
    mean1 <- CompStat$mean1
    n1 <- CompStat$n1
    var1 <- (CompStat$sd1^2)*(1-1/n1)
  }
  x <- seq(.Machine$double.eps, 1-.Machine$double.eps, length = npoints)
  K <- (x*n0*n1*(mean1-mean0)^2)/(2*(x*n0+n1))+
    (x*n0*var0+n1*var1)/2
  s0 <- var0*n0/(n0-1)
  logC <- ((n0-1)/2)*log((n0-1)/2) + ((n0-3)/2)*log(s0) + 
    0.5*(log(n0)-log(2*pi)) - lgamma((n0-1)/2)
  propDen <- (prior$delta.alpha-1)*log(x)+(prior$delta.beta-1)*log(1-x)-
    log(x*n0+n1)/2+lgamma((x*n0+n1-3)/2 +prior$joint.a)-
    ((x*n0+n1-3)/2 + prior$joint.a)*log(K)+x*logC
  deltamin <-  max(0, (3-2*prior$a)/n0)
  propDen <-  ifelse(x>deltamin, propDen, log(.Machine$double.xmin))
  NPPpostLik <- data.frame(x = x, logden = propDen)
  modedelta <- NPPpostLik[which.max(NPPpostLik$logden), 1]
  return(modedelta)
}

#######################################################################################################
#### Original Joint Power Prior Posterior Sampling Function for Normal Data                     #######
#### Use likelihood as Product of n normal Densities based on the original observations         #######
#### Considering the constant 2pi                                                               #######
#######################################################################################################

NormalOPP_MCMC_datalik <- function(Data.Cur, Data.Hist, 
                                   CompStat = list(mean0 = NULL, n0 = NULL, var0 = NULL, mean1 = NULL, n1 = NULL, var1 = NULL), 
                                   prior = list(joint.a = 1.5, delta.alpha = 1, delta.beta = 1), 
                                   MCMCmethod = 'IND', rw.logit.delta = 0.1, 
                                   ind.delta.alpha= 1, ind.delta.beta= 1, nsample = 5000, 
                                   control.mcmc = list(delta.ini = NULL, burnin = 0, thin = 1))
{
  if(missing(CompStat)){
    mean0 <- mean(Data.Hist)
    n0 <- length(Data.Hist)
    var0 <- var(Data.Hist)*(n0-1)/n0
    mean1 <- mean(Data.Cur)
    n1 <- length(Data.Cur)
    var1 <- var(Data.Cur)*(n1-1)/n1
  }else{
    mean0 <- CompStat$mean0
    n0 <- CompStat$n0
    var0 <- CompStat$var0
    mean1 <- CompStat$mean1
    n1 <- CompStat$n1
    var1 <- CompStat$var1
  }

  LogPostNDelta <- function(x){
    K <- (x*n0*n1*(mean1-mean0)^2)/(2*(x*n0+n1))+
      (x*n0*var0+n1*var1)/2
    propDen <- (prior$delta.alpha-1)*log(x)+(prior$delta.beta-1)*log(1-x)-
      log(x*n0+n1)/2+lgamma((x*n0+n1-3)/2 +prior$joint.a)-
      ((x*n0+n1-3)/2 + prior$joint.a)*log(K)-n0*x*0.5*log(2*pi)
    return(propDen)
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
      llik.prop <- LogPostNDelta(delta_prop)
      llik.cur <- LogPostNDelta(delta_cur)
      logr <- min(0, (llik.prop-llik.cur+log(delta_prop)+log(1-delta_prop)-log(delta_cur)-log(1-delta_cur)))
    }
    if(MCMCmethod == 'IND'){
      delta_prop <- rbeta(1, shape1 = ind.delta.alpha, shape2 = ind.delta.beta)
      llik.prop <- LogPostNDelta(delta_prop)
      llik.cur <- LogPostNDelta(delta_cur)
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
  #### Generate mu and sigmasq conditional on delta 
  K = (delta*n0*n1*(mean0-mean1)^2/(delta*n0+n1) + delta*n0*var0 + n1*var1 )/2
  mu = rt_ls(nsample, df= delta*n0+n1+2*prior$joint.a-3, location= (delta*n0*mean0+n1*mean1)/(delta*n0 + n1), 
             scale= sqrt(2*K/( (delta*n0+n1+2*prior$joint.a-3)*(delta*n0+n1) )) )
  
  sigmasq= r_igamma(n = nsample, shape = (delta*n0+n1+2*prior$joint.a-3)/2, rate = K)
  return(list(mu = mu, sigmasq = sigmasq, delta = delta, acceptrate = counter/niter))
}




ModeDeltaNormalOPP_datalik <- function(Data.Cur, Data.Hist,
                                    CompStat = list(n0 = NULL, mean0 = NULL, sd0 = NULL, n1 = NULL, mean1 = NULL, sd1 = NULL),
                                    npoints = 1000,
                                    prior = list(joint.a = 1.5, delta.alpha = 1, delta.beta = 1))
{
  if(missing(CompStat)){
    mean0 <- mean(Data.Hist)
    n0 <- length(Data.Hist)
    var0 <- var(Data.Hist)*(n0-1)/n0
    mean1 <- mean(Data.Cur)
    n1 <- length(Data.Cur)
    var1 <- var(Data.Cur)*(n1-1)/n1
  }else{
    mean0 <- CompStat$mean0
    n0 <- CompStat$n0
    var0 <- (CompStat$sd0^2)*(1-1/n0)
    mean1 <- CompStat$mean1
    n1 <- CompStat$n1
    var1 <- (CompStat$sd1^2)*(1-1/n1)
  }
  x <- seq(.Machine$double.eps, 1-.Machine$double.eps, length = npoints)
  K <- (x*n0*n1*(mean1-mean0)^2)/(2*(x*n0+n1))+
    (x*n0*var0+n1*var1)/2
  propDen <- (prior$delta.alpha-1)*log(x)+(prior$delta.beta-1)*log(1-x)-
    log(x*n0+n1)/2+lgamma((x*n0+n1-3)/2 +prior$joint.a)-
    ((x*n0+n1-3)/2 + prior$joint.a)*log(K)-n0*x*0.5*log(2*pi)
  deltamin <-  max(0, (3-2*prior$a)/n0)
  propDen <-  ifelse(x>deltamin, propDen, log(.Machine$double.xmin))
  NPPpostLik <- data.frame(x = x, logden = propDen)
  modedelta <- NPPpostLik[which.max(NPPpostLik$logden), 1]
  return(modedelta)
}


#######################################################################################################
#### Original Joint Power Prior Posterior Sampling Function for Normal Data                     #######
#### Use likelihood as Product of n normal Densities based on the original observations         #######
#### Including an arbitrary constant (2pi)^(n0/2)exp(200)                                       #######
#######################################################################################################

NormalOPP_MCMC_arbitrary <- function(Data.Cur, Data.Hist, 
                                     CompStat = list(mean0 = NULL, n0 = NULL, var0 = NULL, mean1 = NULL, n1 = NULL, var1 = NULL), 
                                     prior = list(joint.a = 1.5, delta.alpha = 1, delta.beta = 1), 
                                     MCMCmethod = 'IND', rw.logit.delta = 0.1, 
                                     ind.delta.alpha= 1, ind.delta.beta= 1, nsample = 5000, 
                                     control.mcmc = list(delta.ini = NULL, burnin = 0, thin = 1))
{
  if(missing(CompStat)){
    mean0 <- mean(Data.Hist)
    n0 <- length(Data.Hist)
    var0 <- var(Data.Hist)*(n0-1)/n0
    mean1 <- mean(Data.Cur)
    n1 <- length(Data.Cur)
    var1 <- var(Data.Cur)*(n1-1)/n1
  }else{
    mean0 <- CompStat$mean0
    n0 <- CompStat$n0
    var0 <- CompStat$var0
    mean1 <- CompStat$mean1
    n1 <- CompStat$n1
    var1 <- CompStat$var1
  }
  LogPostNDelta <- function(x){
    K <- (x*n0*n1*(mean1-mean0)^2)/(2*(x*n0+n1))+
      (x*n0*var0+n1*var1)/2
    propDen <- (prior$delta.alpha-1)*log(x)+(prior$delta.beta-1)*log(1-x)-
      log(x*n0+n1)/2+lgamma((x*n0+n1-3)/2 +prior$joint.a)-
      ((x*n0+n1-3)/2 + prior$joint.a)*log(K)+200*x
    return(propDen)
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
      llik.prop <- LogPostNDelta(delta_prop)
      llik.cur <- LogPostNDelta(delta_cur)
      logr <- min(0, (llik.prop-llik.cur+log(delta_prop)+log(1-delta_prop)-log(delta_cur)-log(1-delta_cur)))
    }
    if(MCMCmethod == 'IND'){
      delta_prop <- rbeta(1, shape1 = ind.delta.alpha, shape2 = ind.delta.beta)
      llik.prop <- LogPostNDelta(delta_prop)
      llik.cur <- LogPostNDelta(delta_cur)
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
  #### vectorized generate mu and sigmasq conditional on delta 
  K = (delta*n0*n1*(mean0-mean1)^2/(delta*n0+n1) + delta*n0*var0 + n1*var1 )/2
  mu = rt_ls(nsample, df= delta*n0+n1+2*prior$joint.a-3, location= (delta*n0*mean0+n1*mean1)/(delta*n0 + n1), 
             scale= sqrt(2*K/( (delta*n0+n1+2*prior$joint.a-3)*(delta*n0+n1) )) )
  
  sigmasq= r_igamma(n = nsample, shape = (delta*n0+n1+2*prior$joint.a-3)/2, rate = K)
  return(list(mu = mu, sigmasq = sigmasq, delta = delta, acceptrate = counter/niter))
}



ModeDeltaNormalOPP_arbitrary <- function(Data.Cur, Data.Hist,
                                       CompStat = list(n0 = NULL, mean0 = NULL, sd0 = NULL, n1 = NULL, mean1 = NULL, sd1 = NULL),
                                       npoints = 1000,
                                       prior = list(joint.a = 1.5, delta.alpha = 1, delta.beta = 1))
{
  if(missing(CompStat)){
    mean0 <- mean(Data.Hist)
    n0 <- length(Data.Hist)
    var0 <- var(Data.Hist)*(n0-1)/n0
    mean1 <- mean(Data.Cur)
    n1 <- length(Data.Cur)
    var1 <- var(Data.Cur)*(n1-1)/n1
  }else{
    mean0 <- CompStat$mean0
    n0 <- CompStat$n0
    var0 <- (CompStat$sd0^2)*(1-1/n0)
    mean1 <- CompStat$mean1
    n1 <- CompStat$n1
    var1 <- (CompStat$sd1^2)*(1-1/n1)
  }
  x <- seq(.Machine$double.eps, 1-.Machine$double.eps, length = npoints)
  K <- (x*n0*n1*(mean1-mean0)^2)/(2*(x*n0+n1))+(x*n0*var0+n1*var1)/2
  propDen <- (prior$delta.alpha-1)*log(x)+(prior$delta.beta-1)*log(1-x)-
    log(x*n0+n1)/2+lgamma((x*n0+n1-3)/2 +prior$joint.a)-
    ((x*n0+n1-3)/2 + prior$joint.a)*log(K)+200*x
  deltamin <-  max(0, (3-2*prior$a)/n0)
  propDen <-  ifelse(x>deltamin, propDen, log(.Machine$double.xmin))
  NPPpostLik <- data.frame(x = x, logden = propDen)
  modedelta <- NPPpostLik[which.max(NPPpostLik$logden), 1]
  return(modedelta)
}


#######################################################################################################
#### Normalized Power Prior Posterior Sampling Function for Normal Data                         #######
#######################################################################################################

NormalNPP_MCMC <- function(Data.Cur, Data.Hist, 
                           CompStat = list(mean0 = NULL, n0 = NULL, var0 = NULL, mean1 = NULL, n1 = NULL, var1 = NULL), 
                           prior = list(joint.a = 1.5, delta.alpha = 1, delta.beta = 1), 
                           MCMCmethod = 'IND', rw.logit.delta = 0.1, 
                           ind.delta.alpha= 1, ind.delta.beta= 1, nsample = 5000, 
                           control.mcmc = list(delta.ini = NULL, burnin = 0, thin = 1))
{
  if(missing(CompStat)){
    mean0 <- mean(Data.Hist)
    n0 <- length(Data.Hist)
    var0 <- var(Data.Hist)*(n0-1)/n0
    mean1 <- mean(Data.Cur)
    n1 <- length(Data.Cur)
    var1 <- var(Data.Cur)*(n1-1)/n1
  }else{
    mean0 <- CompStat$mean0
    n0 <- CompStat$n0
    var0 <- CompStat$var0
    mean1 <- CompStat$mean1
    n1 <- CompStat$n1
    var1 <- CompStat$var1
  }
  #### Prior pi(mu,sigmasq) propto (1/sigmasq)^a
  #### Normalized Power Prior for Normal Log Marginal Posterior (unnormalized) of Delta
  deltamin <- (3-2*prior$joint.a)/n0
  LogPostNDelta <- function(x){
    lden = (x*n0/2+prior$joint.a+prior$delta.alpha-2)*log(x)+(prior$delta.beta-1)*log(1-x)+
      lgamma((x*n0+n1-3)/2+prior$joint.a)-lgamma((x*n0-3)/2+prior$joint.a)-
      ((x*n0+n1-3)/2+prior$joint.a)*log((x*n1*(mean0-mean1)^2)/((x*n0+n1)*var0)+x+(n1*var1)/(n0*var0))-
      log(x*n0+n1)/2
    out = ifelse(x>deltamin, lden, log(.Machine$double.xmin))
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
      llik.prop <- LogPostNDelta(delta_prop)
      llik.cur <- LogPostNDelta(delta_cur)
      logr <- min(0, (llik.prop-llik.cur+log(delta_prop)+log(1-delta_prop)-log(delta_cur)-log(1-delta_cur)))
    }
    if(MCMCmethod == 'IND'){
      delta_prop <- rbeta(1, shape1 = ind.delta.alpha, shape2 = ind.delta.beta)
      llik.prop <- LogPostNDelta(delta_prop)
      llik.cur <- LogPostNDelta(delta_cur)
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
  #### vectorized generate mu and sigmasq conditional on delta 
  K = (delta*n0*n1*(mean0-mean1)^2/(delta*n0+n1) + delta*n0*var0 + n1*var1 )/2
  mu = rt_ls(nsample, df= delta*n0+n1+2*prior$joint.a-3, location= (delta*n0*mean0+n1*mean1)/(delta*n0 + n1), 
             scale= sqrt(2*K/( (delta*n0+n1+2*prior$joint.a-3)*(delta*n0+n1) )) )
  sigmasq= r_igamma(n = nsample, shape = (delta*n0+n1+2*prior$joint.a-3)/2, rate = K)
  meanmu <- mean(mu)
  meansigmasq <- mean(sigmasq)
  ### DIC without constant term
  D <- -2*(-n1*log(sigmasq)/2 -n1*(var1 + (mu-mean1)^2 )/(2*sigmasq))
  Dbar <- -2*(-n1*log(meansigmasq)/2 -n1*(var1 + (meanmu-mean1)^2 )/(2*meansigmasq))
  DIC <- 2*mean(D)-Dbar
  return(list(mu = mu, sigmasq = sigmasq, delta = delta, 
              acceptrate = counter/niter, DIC = DIC))
}


ModeDeltaNormalNPP <- function(Data.Cur, Data.Hist,
                                         CompStat = list(n0 = NULL, mean0 = NULL, sd0 = NULL, n1 = NULL, mean1 = NULL, sd1 = NULL),
                                         npoints = 1000,
                                         prior = list(joint.a = 1.5, delta.alpha = 1, delta.beta = 1))
{
  if(missing(CompStat)){
    mean0 <- mean(Data.Hist)
    n0 <- length(Data.Hist)
    var0 <- var(Data.Hist)*(n0-1)/n0
    mean1 <- mean(Data.Cur)
    n1 <- length(Data.Cur)
    var1 <- var(Data.Cur)*(n1-1)/n1
  }else{
    mean0 <- CompStat$mean0
    n0 <- CompStat$n0
    var0 <- (CompStat$sd0^2)*(1-1/n0)
    mean1 <- CompStat$mean1
    n1 <- CompStat$n1
    var1 <- (CompStat$sd1^2)*(1-1/n1)
  }
  x <- seq(.Machine$double.eps, 1-.Machine$double.eps, length = npoints)
  propDen = (x*n0/2+prior$joint.a+prior$delta.alpha-2)*log(x)+(prior$delta.beta-1)*log(1-x)+
    lgamma((x*n0+n1-3)/2+prior$joint.a)-lgamma((x*n0-3)/2+prior$joint.a)-
    ((x*n0+n1-3)/2+prior$joint.a)*log((x*n1*(mean0-mean1)^2)/((x*n0+n1)*var0)+x+(n1*var1)/(n0*var0))-
    log(x*n0+n1)/2
  deltamin <-  max(0, (3-2*prior$a)/n0)
  propDen <-  ifelse(x>deltamin, propDen, log(.Machine$double.xmin))
  NPPpostLik <- data.frame(x = x, logden = propDen)
  modedelta <- NPPpostLik[which.max(NPPpostLik$logden), 1]
  return(modedelta)
}


#####################################################################################################
##### Functions to calculate the  marginal posterior density of \delta after normalization ##########
##### Used to generate the marginal posterior plot                                         ##########
#####################################################################################################

DeltaNormalOPP_den_suff <- function(Data.Cur, Data.Hist,
                                    CompStat = list(n0 = NULL, mean0 = NULL, sd0 = NULL, n1 = NULL, mean1 = NULL, sd1 = NULL),
                                    npoints = 1000,
                                    prior = list(joint.a = 1.5, delta.alpha = 1, delta.beta = 1))
{
  if(missing(CompStat)){
    mean0 <- mean(Data.Hist)
    n0 <- length(Data.Hist)
    var0 <- var(Data.Hist)*(n0-1)/n0
    mean1 <- mean(Data.Cur)
    n1 <- length(Data.Cur)
    var1 <- var(Data.Cur)*(n1-1)/n1
  }else{
    mean0 <- CompStat$mean0
    n0 <- CompStat$n0
    var0 <- (CompStat$sd0^2)*(1-1/n0)
    mean1 <- CompStat$mean1
    n1 <- CompStat$n1
    var1 <- (CompStat$sd1^2)*(1-1/n1)
  }
  x <- seq(.Machine$double.eps, 1-.Machine$double.eps, length = npoints)
  K <- (x*n0*n1*(mean1-mean0)^2)/(2*(x*n0+n1))+
    (x*n0*var0+n1*var1)/2
  s0 <- var0*n0/(n0-1)
  logC <- ((n0-1)/2)*log((n0-1)/2) + ((n0-3)/2)*log(s0) + 
    0.5*(log(n0)-log(2*pi)) - lgamma((n0-1)/2)
  propDen <- (prior$delta.alpha-1)*log(x)+(prior$delta.beta-1)*log(1-x)-
    log(x*n0+n1)/2+lgamma((x*n0+n1-3)/2 +prior$joint.a)-
    ((x*n0+n1-3)/2 + prior$joint.a)*log(K)+x*logC
  deltamin <-  max(0, (3-2*prior$a)/n0)
  propDen <-  ifelse(x>deltamin, propDen, log(.Machine$double.xmin))
  NPPpostLik <- data.frame(x = x, logden = propDen)
  denum <- sum((1/npoints)*(exp(propDen[-1]) + exp(propDen[-npoints]))/2)  ## Midpoint rule to do normalization
  normden <- data.frame(x = x, den = exp(propDen)/denum)
  return(normden)
}


DeltaNormalOPP_den_datalik <- function(Data.Cur, Data.Hist,
                                       CompStat = list(n0 = NULL, mean0 = NULL, sd0 = NULL, n1 = NULL, mean1 = NULL, sd1 = NULL),
                                       npoints = 1000,
                                       prior = list(joint.a = 1.5, delta.alpha = 1, delta.beta = 1))
{
  if(missing(CompStat)){
    mean0 <- mean(Data.Hist)
    n0 <- length(Data.Hist)
    var0 <- var(Data.Hist)*(n0-1)/n0
    mean1 <- mean(Data.Cur)
    n1 <- length(Data.Cur)
    var1 <- var(Data.Cur)*(n1-1)/n1
  }else{
    mean0 <- CompStat$mean0
    n0 <- CompStat$n0
    var0 <- (CompStat$sd0^2)*(1-1/n0)
    mean1 <- CompStat$mean1
    n1 <- CompStat$n1
    var1 <- (CompStat$sd1^2)*(1-1/n1)
  }
  x <- seq(.Machine$double.eps, 1-.Machine$double.eps, length = npoints)
  K <- (x*n0*n1*(mean1-mean0)^2)/(2*(x*n0+n1))+
    (x*n0*var0+n1*var1)/2
  propDen <- (prior$delta.alpha-1)*log(x)+(prior$delta.beta-1)*log(1-x)-
    log(x*n0+n1)/2+lgamma((x*n0+n1-3)/2 +prior$joint.a)-
    ((x*n0+n1-3)/2 + prior$joint.a)*log(K)-n0*x*0.5*log(2*pi)
  deltamin <-  max(0, (3-2*prior$a)/n0)
  propDen <-  ifelse(x>deltamin, propDen, log(.Machine$double.xmin))
  NPPpostLik <- data.frame(x = x, logden = propDen)
  denum <- sum((1/npoints)*(exp(propDen[-1]) + exp(propDen[-npoints]))/2)  ## Midpoint rule to do normalization
  normden <- data.frame(x = x, den = exp(propDen)/denum)
  return(normden)
}


DeltaNormalOPP_den_arbitrary <- function(Data.Cur, Data.Hist,
                                         CompStat = list(n0 = NULL, mean0 = NULL, sd0 = NULL, n1 = NULL, mean1 = NULL, sd1 = NULL),
                                         npoints = 1000,
                                         prior = list(joint.a = 1.5, delta.alpha = 1, delta.beta = 1))
{
  if(missing(CompStat)){
    mean0 <- mean(Data.Hist)
    n0 <- length(Data.Hist)
    var0 <- var(Data.Hist)*(n0-1)/n0
    mean1 <- mean(Data.Cur)
    n1 <- length(Data.Cur)
    var1 <- var(Data.Cur)*(n1-1)/n1
  }else{
    mean0 <- CompStat$mean0
    n0 <- CompStat$n0
    var0 <- (CompStat$sd0^2)*(1-1/n0)
    mean1 <- CompStat$mean1
    n1 <- CompStat$n1
    var1 <- (CompStat$sd1^2)*(1-1/n1)
  }
  x <- seq(.Machine$double.eps, 1-.Machine$double.eps, length = npoints)
  K <- (x*n0*n1*(mean1-mean0)^2)/(2*(x*n0+n1))+
    (x*n0*var0+n1*var1)/2
  propDen <- (prior$delta.alpha-1)*log(x)+(prior$delta.beta-1)*log(1-x)-
    log(x*n0+n1)/2+lgamma((x*n0+n1-3)/2 +prior$joint.a)-
    ((x*n0+n1-3)/2 + prior$joint.a)*log(K)+200*x
  deltamin <-  max(0, (3-2*prior$a)/n0)
  propDen <-  ifelse(x>deltamin, propDen, log(.Machine$double.xmin))
  NPPpostLik <- data.frame(x = x, logden = propDen)
  denum <- sum((1/npoints)*(exp(propDen[-1]) + exp(propDen[-npoints]))/2)  ## Midpoint rule to do normalization
  normden <- data.frame(x = x, den = exp(propDen)/denum)
  return(normden)
}


DeltaNormalNPP_den <- function(Data.Cur, Data.Hist,
                               CompStat = list(n0 = NULL, mean0 = NULL, sd0 = NULL, n1 = NULL, mean1 = NULL, sd1 = NULL),
                               npoints = 1000,
                               prior = list(joint.a = 1.5, delta.alpha = 1, delta.beta = 1))
{
  if(missing(CompStat)){
    mean0 <- mean(Data.Hist)
    n0 <- length(Data.Hist)
    var0 <- var(Data.Hist)*(n0-1)/n0
    mean1 <- mean(Data.Cur)
    n1 <- length(Data.Cur)
    var1 <- var(Data.Cur)*(n1-1)/n1
  }else{
    mean0 <- CompStat$mean0
    n0 <- CompStat$n0
    var0 <- (CompStat$sd0^2)*(1-1/n0)
    mean1 <- CompStat$mean1
    n1 <- CompStat$n1
    var1 <- (CompStat$sd1^2)*(1-1/n1)
  }
  x <- seq(.Machine$double.eps, 1-.Machine$double.eps, length = npoints)
  propDen <- (x*n0/2+prior$joint.a+prior$delta.alpha-2)*log(x)+(prior$delta.beta-1)*log(1-x)+
    lgamma((x*n0+n1-3)/2+prior$joint.a)-lgamma((x*n0-3)/2+prior$joint.a)-
    ((x*n0+n1-3)/2+prior$joint.a)*log((x*n1*(mean0-mean1)^2)/((x*n0+n1)*var0)+x+(n1*var1)/(n0*var0))-
    log(x*n0+n1)/2
  deltamin <-  max(0, (3-2*prior$a)/n0)
  propDen <-  ifelse(x>deltamin, propDen, log(.Machine$double.xmin))
  NPPpostLik <- data.frame(x = x, logden = propDen)
  denum <- sum((1/npoints)*(exp(propDen[-1]) + exp(propDen[-npoints]))/2)  ## Midpoint rule to do normalization
  normden <- data.frame(x = x, den = exp(propDen)/denum)
  return(normden)
}


#######################################################################################################
#### Preparing some functions to calculate the sampling of 10% quantile                         #######
#######################################################################################################

PHNPPL <- function(mean0, n0, sd0, mean1, n1, sd1, nsample){
  POST <- NormalNPP_MCMC(CompStat = list(mean0 = mean0, n0 = n0, var0 = (sd0^2)*(1-1/n0), 
                                         mean1 = mean1, n1 = n1, var1 = (sd1^2)*(1-1/n1)), 
                 prior = list(joint.a = 1, delta.alpha = 1, delta.beta = 1), MCMCmethod = 'IND', 
                 rw.logit.delta = 0.1, ind.delta.alpha= 1, ind.delta.beta= 1, nsample = nsample, 
                 control.mcmc = list(delta.ini = NULL, burnin = 2000, thin = 5))
  POSTL <- POST$mu+qnorm(0.1)*sqrt(POST$sigmasq)
  return(list(L = POSTL, delta = POST$delta))
}

PHOPPL1 <- function(mean0, n0, sd0, mean1, n1, sd1, nsample){
  POST <- NormalOPP_MCMC_suff(CompStat = list(mean0 = mean0, n0 = n0, var0 = (sd0^2)*(1-1/n0), 
                                         mean1 = mean1, n1 = n1, var1 = (sd1^2)*(1-1/n1)), 
                         prior = list(joint.a = 1, delta.alpha = 1, delta.beta = 1), MCMCmethod = 'RW', 
                         rw.logit.delta = 0.3, ind.delta.alpha= 1, ind.delta.beta= 1, nsample = nsample, 
                         control.mcmc = list(delta.ini = NULL, burnin = 2000, thin = 5))
  POSTL <- POST$mu+qnorm(0.1)*sqrt(POST$sigmasq)
  return(list(L = POSTL, delta = POST$delta))
}

PHOPPL2 <- function(mean0, n0, sd0, mean1, n1, sd1, nsample){
  POST <- NormalOPP_MCMC_datalik(CompStat = list(mean0 = mean0, n0 = n0, var0 = (sd0^2)*(1-1/n0), 
                                         mean1 = mean1, n1 = n1, var1 = (sd1^2)*(1-1/n1)), 
                         prior = list(joint.a = 1, delta.alpha = 1, delta.beta = 1), MCMCmethod = 'RW', 
                         rw.logit.delta = 0.3, ind.delta.alpha= 1, ind.delta.beta= 1, nsample = nsample, 
                         control.mcmc = list(delta.ini = NULL, burnin = 2000, thin = 5))
  POSTL <- POST$mu+qnorm(0.1)*sqrt(POST$sigmasq)
  return(list(L = POSTL, delta = POST$delta))
}

PHOPPL3 <- function(mean0, n0, sd0, mean1, n1, sd1, nsample){
  POST <- NormalOPP_MCMC_arbitrary(CompStat = list(mean0 = mean0, n0 = n0, var0 = (sd0^2)*(1-1/n0), 
                                         mean1 = mean1, n1 = n1, var1 = (sd1^2)*(1-1/n1)), 
                         prior = list(joint.a = 1, delta.alpha = 1, delta.beta = 1), MCMCmethod = 'RW', 
                         rw.logit.delta = 0.3, ind.delta.alpha= 1, ind.delta.beta= 1, nsample = nsample, 
                         control.mcmc = list(delta.ini = NULL, burnin = 2000, thin = 5))
  POSTL <- POST$mu+qnorm(0.1)*sqrt(POST$sigmasq)
  return(list(L = POSTL, delta = POST$delta))
}

PHNBL <- function(mean1, n1, sd1, nsample){
  var1 = (sd1^2)*(1-1/n1)
  mu=  rt(nsample, df = n1-1)*sqrt(var1/(n1-1))+mean1
  sigmasq= rinvchisq(nsample, df = n1-1, scale = sd1^2)
  POSTL <- mu+qnorm(0.1)*sqrt(sigmasq)
  return(POSTL)
}

FunL <- function(x) c(round(length(x[x>6])/length(x),3), round(sd(x),2))


##################################################################################
################ Model Fitting Starts here                        ################
##################################################################################

### OPP1: density based on sufficient statistics; normal times gamma
set.seed(1234)
OPPL1_A <- PHOPPL1(mean0 = 7.05, n0 = 62, sd0 = 0.47, mean1 = 6.91, n1 = 16, sd1 = 0.90, nsample = 20000)
OPPL1_B <- PHOPPL1(mean0 = 6.73, n0 = 31, sd0 = 0.71, mean1 = 6.78, n1 = 12, sd1 = 1.03, nsample = 20000)
OPPL1_C <- PHOPPL1(mean0 = 6.95, n0 = 84, sd0 = 0.49, mean1 = 6.43, n1 = 24, sd1 = 0.88, nsample = 20000)
OPPL1_D <- PHOPPL1(mean0 = 7.88, n0 = 75, sd0 = 0.67, mean1 = 7.87, n1 = 21, sd1 = 1.11, nsample = 20000)

FunL(OPPL1_A$L); FunL(OPPL1_B$L); FunL(OPPL1_C$L); FunL(OPPL1_D$L)
plot(OPPL1_A$L,type ='l');acf(OPPL1_A$L)
plot(OPPL1_B$L,type ='l');acf(OPPL1_B$L)
plot(OPPL1_C$L,type ='l');acf(OPPL1_C$L)
plot(OPPL1_D$L,type ='l');acf(OPPL1_D$L)
#[1] 0.385 0.310
#[1] 0.051 0.300
#[1] 0.003 0.250
#[1] 0.959 0.300

mean(OPPL1_A$delta); mean(OPPL1_B$delta); mean(OPPL1_C$delta); mean(OPPL1_D$delta)
#[1] 0.154
#[1] 0.464
#[1] 0.038
#[1] 0.226


#### OPP2: based on the loglikelihood of original data, including the 2pi
set.seed(1234)
OPPL2_A <- PHOPPL2(mean0 = 7.05, n0 = 62, sd0 = 0.47, mean1 = 6.91, n1 = 16, sd1 = 0.90, nsample = 20000)
OPPL2_B <- PHOPPL2(mean0 = 6.73, n0 = 31, sd0 = 0.71, mean1 = 6.78, n1 = 12, sd1 = 1.03, nsample = 20000)
OPPL2_C <- PHOPPL2(mean0 = 6.95, n0 = 84, sd0 = 0.49, mean1 = 6.43, n1 = 24, sd1 = 0.88, nsample = 20000)
OPPL2_D <- PHOPPL2(mean0 = 7.88, n0 = 75, sd0 = 0.67, mean1 = 7.87, n1 = 21, sd1 = 1.11, nsample = 20000)

FunL(OPPL2_A$L); FunL(OPPL2_B$L); FunL(OPPL2_C$L); FunL(OPPL2_D$L)
plot(OPPL2_A$L,type ='l');acf(OPPL2_A$L)
plot(OPPL2_B$L,type ='l');acf(OPPL2_B$L)
plot(OPPL2_C$L,type ='l');acf(OPPL2_C$L)
plot(OPPL2_D$L,type ='l');acf(OPPL2_D$L)
#[1] 0.201 0.320
#[1] 0.07 0.45
#[1] 0.002 0.250
#[1] 0.886 0.350

round(mean(OPPL2_A$delta), 3); round(mean(OPPL2_B$delta), 3); round(mean(OPPL2_C$delta), 3); round(mean(OPPL2_D$delta), 3)
#[1] 0.016
#[1] 0.026
#[1] 0.01
#[1] 0.011

#### OPP3: with arbitrary constant exp(200)*(2pi)^(n0/2)
set.seed(123)
OPPL3_A <- PHOPPL3(mean0 = 7.05, n0 = 62, sd0 = 0.47, mean1 = 6.91, n1 = 16, sd1 = 0.90, nsample = 20000)
OPPL3_B <- PHOPPL3(mean0 = 6.73, n0 = 31, sd0 = 0.71, mean1 = 6.78, n1 = 12, sd1 = 1.03, nsample = 20000)
OPPL3_C <- PHOPPL3(mean0 = 6.95, n0 = 84, sd0 = 0.49, mean1 = 6.43, n1 = 24, sd1 = 0.88, nsample = 20000)
OPPL3_D <- PHOPPL3(mean0 = 7.88, n0 = 75, sd0 = 0.67, mean1 = 7.87, n1 = 21, sd1 = 1.11, nsample = 20000)

FunL(OPPL3_A$L); FunL(OPPL3_B$L); FunL(OPPL3_C$L); FunL(OPPL3_D$L)
plot(OPPL3_A$L,type ='l');acf(OPPL3_A$L)
plot(OPPL3_B$L,type ='l');acf(OPPL3_B$L)
plot(OPPL3_C$L,type ='l');acf(OPPL3_C$L)
plot(OPPL3_D$L,type ='l');acf(OPPL3_D$L)
#[1] 0.997 0.090
#[1] 0.033 0.170
#[1] 0.592 0.080
#[1] 1.00 0.11

mean(OPPL3_A$delta)
mean(OPPL3_B$delta)
mean(OPPL3_C$delta)
mean(OPPL3_D$delta)



#### NPP: 
set.seed(1234)
NPPL_A <- PHNPPL(mean0 = 7.05, n0 = 62, sd0 = 0.47, mean1 = 6.91, n1 = 16, sd1 = 0.90, nsample = 20000)
NPPL_B <- PHNPPL(mean0 = 6.73, n0 = 31, sd0 = 0.71, mean1 = 6.78, n1 = 12, sd1 = 1.03, nsample = 20000)
NPPL_C <- PHNPPL(mean0 = 6.95, n0 = 84, sd0 = 0.49, mean1 = 6.43, n1 = 24, sd1 = 0.88, nsample = 20000)
NPPL_D <- PHNPPL(mean0 = 7.88, n0 = 75, sd0 = 0.67, mean1 = 7.87, n1 = 21, sd1 = 1.11, nsample = 20000)

FunL(NPPL_A$L); FunL(NPPL_B$L); FunL(NPPL_C$L); FunL(NPPL_D$L)
plot(NPPL_A$L,type ='l');acf(NPPL_A$L)
plot(NPPL_B$L,type ='l');acf(NPPL_B$L)
plot(NPPL_C$L,type ='l');acf(NPPL_C$L)
plot(NPPL_D$L,type ='l');acf(NPPL_D$L)

#[1] 0.488 0.260
#[1] 0.047 0.260
#[1] 0.004 0.240
#[1] 0.986 0.250

round(mean(NPPL_A$delta), 3); round(mean(NPPL_B$delta),3); round(mean(NPPL_C$delta), 3); round(mean(NPPL_D$delta), 3)
#[1] 0.209
#[1] 0.551
#[1] 0.085
#[1] 0.302

round(ModeDeltaNormalNPP(CompStat = list(mean0 = 7.05, n0 = 62, sd0 = 0.47, mean1 = 6.91, n1 = 16, sd1 = 0.90), npoints = 1000), 3)
round(ModeDeltaNormalNPP(CompStat = list(mean0 = 6.73, n0 = 31, sd0 = 0.71, mean1 = 6.78, n1 = 12, sd1 = 1.03), npoints = 1000), 3)
round(ModeDeltaNormalNPP(CompStat = list(mean0 = 6.95, n0 = 84, sd0 = 0.49, mean1 = 6.43, n1 = 24, sd1 = 0.88), npoints = 1000), 3)
round(ModeDeltaNormalNPP(CompStat = list(mean0 = 7.88, n0 = 75, sd0 = 0.67, mean1 = 7.87, n1 = 21, sd1 = 1.11), npoints = 1000), 3)

#[1] 0.06
#[1] 0.482
#[1] 0.031
#[1] 0.081


#### Reference Prior Only, no historical data
set.seed(123)
NBL_A <- PHNBL(mean1 = 6.91, n1 = 16, sd1 = 0.90, nsample = 20000)
NBL_B <- PHNBL(mean1 = 6.78, n1 = 12, sd1 = 1.03, nsample = 20000)
NBL_C <- PHNBL(mean1 = 6.43, n1 = 24, sd1 = 0.88, nsample = 20000)
NBL_D <- PHNBL(mean1 = 7.87, n1 = 21, sd1 = 1.11, nsample = 20000)

FunL(NBL_A); FunL(NBL_B); FunL(NBL_C); FunL(NBL_D)
plot(NBL_A,type ='l');acf(NBL_A)
plot(NBL_B,type ='l');acf(NBL_B)
plot(NBL_C,type ='l');acf(NBL_C)
plot(NBL_D,type ='l');acf(NBL_D)

#[1] 0.177 0.340
#[1] 0.069 0.470
#[1] 0.001 0.260
#[1] 0.865 0.360



#######################################################################################
#######################################################################################
#### The Plot of the Marginal Posterior of delta for four sites using four priors  ####
#######################################################################################
#######################################################################################


Delta_A_opp1 = DeltaNormalOPP_den_suff(CompStat = list(mean0 = 7.05, n0 = 62, sd0 = 0.47, mean1 = 6.91, n1 = 16, sd1 = 0.90), npoints = 1000)
Delta_A_opp2 = DeltaNormalOPP_den_datalik(CompStat = list(mean0 = 7.05, n0 = 62, sd0 = 0.47, mean1 = 6.91, n1 = 16, sd1 = 0.90), npoints = 1000)
Delta_A_opp3 = DeltaNormalOPP_den_arbitrary(CompStat = list(mean0 = 7.05, n0 = 62, sd0 = 0.47, mean1 = 6.91, n1 = 16, sd1 = 0.90), npoints = 1000)
Delta_A_npp = DeltaNormalNPP_den(CompStat = list(mean0 = 7.05, n0 = 62, sd0 = 0.47, mean1 = 6.91, n1 = 16, sd1 = 0.90), npoints = 1000)

Delta_B_opp1 = DeltaNormalOPP_den_suff(CompStat = list(mean0 = 6.73, n0 = 31, sd0 = 0.71, mean1 = 6.78, n1 = 12, sd1 = 1.03), npoints = 1000)
Delta_B_opp2 = DeltaNormalOPP_den_datalik(CompStat = list(mean0 = 6.73, n0 = 31, sd0 = 0.71, mean1 = 6.78, n1 = 12, sd1 = 1.03), npoints = 1000)
Delta_B_opp3 = DeltaNormalOPP_den_arbitrary(CompStat = list(mean0 = 6.73, n0 = 31, sd0 = 0.71, mean1 = 6.78, n1 = 12, sd1 = 1.03), npoints = 1000)
Delta_B_npp = DeltaNormalNPP_den(CompStat = list(mean0 = 6.73, n0 = 31, sd0 = 0.71, mean1 = 6.78, n1 = 12, sd1 = 1.03), npoints = 1000)

Delta_C_opp1 = DeltaNormalOPP_den_suff(CompStat = list(mean0 = 6.95, n0 = 84, sd0 = 0.49, mean1 = 6.43, n1 = 24, sd1 = 0.88), npoints = 1000)
Delta_C_opp2 = DeltaNormalOPP_den_datalik(CompStat = list(mean0 = 6.95, n0 = 84, sd0 = 0.49, mean1 = 6.43, n1 = 24, sd1 = 0.88), npoints = 1000)
Delta_C_opp3 = DeltaNormalOPP_den_arbitrary(CompStat = list(mean0 = 6.95, n0 = 84, sd0 = 0.49, mean1 = 6.43, n1 = 24, sd1 = 0.88), npoints = 1000)
Delta_C_npp = DeltaNormalNPP_den(CompStat = list(mean0 = 6.95, n0 = 84, sd0 = 0.49, mean1 = 6.43, n1 = 24, sd1 = 0.88), npoints = 1000)

Delta_D_opp1 = DeltaNormalOPP_den_suff(CompStat = list(mean0 = 7.88, n0 = 75, sd0 = 0.67, mean1 = 7.87, n1 = 21, sd1 = 1.11), npoints = 1000)
Delta_D_opp2 = DeltaNormalOPP_den_datalik(CompStat = list(mean0 = 7.88, n0 = 75, sd0 = 0.67, mean1 = 7.87, n1 = 21, sd1 = 1.11), npoints = 1000)
Delta_D_opp3 = DeltaNormalOPP_den_arbitrary(CompStat = list(mean0 = 7.88, n0 = 75, sd0 = 0.67, mean1 = 7.87, n1 = 21, sd1 = 1.11), npoints = 1000)
Delta_D_npp = DeltaNormalNPP_den(CompStat = list(mean0 = 7.88, n0 = 75, sd0 = 0.67, mean1 = 7.87, n1 = 21, sd1 = 1.11), npoints = 1000)


### Start the ggplot 
Delta_A_df = data.frame(rbind(Delta_A_npp, Delta_A_opp1, Delta_A_opp2, Delta_A_opp3), 
                        variable = rep(c("ANPP", "JPP1", "JPP2", "JPP3"), each = 1000))
Delta_B_df = data.frame(rbind(Delta_B_npp, Delta_B_opp1, Delta_B_opp2, Delta_B_opp3), 
                        variable = rep(c("ANPP", "JPP1", "JPP2", "JPP3"), each = 1000))
Delta_C_df = data.frame(rbind(Delta_C_npp, Delta_C_opp1, Delta_C_opp2, Delta_C_opp3), 
                        variable = rep(c("ANPP", "JPP1", "JPP2", "JPP3"), each = 1000))
Delta_D_df = data.frame(rbind(Delta_D_npp, Delta_D_opp1, Delta_D_opp2, Delta_D_opp3), 
                        variable = rep(c("ANPP", "JPP1", "JPP2", "JPP3"), each = 1000))

FADelta <- ggplot(Delta_A_df, aes(x = x, y = den, col = variable, linetype = variable)) +
  geom_line(lwd = 0.6)+
  scale_linetype_manual(values = c("solid", "dashed", "dotdash", "dotdash"), 
                        labels = c("NPP", "JPP1", "JPP2", "JPP3")) +
  scale_color_manual(values = c("red", "blue", "forestgreen", "brown"), 
                     labels = c("NPP", "JPP1", "JPP2", "JPP3")) +
  theme_bw()+theme(legend.title=element_blank(), 
                   axis.text=element_text(size=9),
                   axis.title=element_text(size=8))+
  labs(x=expression(delta), y= 'Density')+
  scale_y_continuous(breaks=c(0, 1, 2, 3, 4, 5),limits = c(0, 5.1))+
  scale_x_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0, 1.02))+
  theme(legend.key.width = unit(2.5, "line"), legend.position="bottom")+
  guides(col = guide_legend(nrow=2.5), linetype = guide_legend(nrow = 1))

FBDelta <- ggplot(Delta_B_df, aes(x = x, y = den, col = variable, linetype = variable)) +
  geom_line(lwd = 0.6)+
  scale_linetype_manual(values = c("solid", "dashed", "dotdash", "dotdash"), 
                        labels = c("NPP", "JPP1", "JPP2", "JPP3")) +
  scale_color_manual(values = c("red", "blue", "forestgreen", "brown"), 
                     labels = c("NPP", "JPP1", "JPP2", "JPP3")) +
  theme_bw()+theme(legend.title=element_blank(), 
                   axis.text=element_text(size=9),
                   axis.title=element_text(size=8))+
  labs(x=expression(delta), y= 'Density')+
  scale_y_continuous(breaks=c(0, 0.5, 1, 1.5, 2),limits = c(0, 2.1))+
  scale_x_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0, 1.02))+
  theme(legend.key.width = unit(2.5, "line"), legend.position="bottom")+
  guides(col = guide_legend(nrow=2.5), linetype = guide_legend(nrow = 1))


FCDelta <- ggplot(Delta_C_df, aes(x = x, y = den, col = variable, linetype = variable)) +
  geom_line(lwd = 0.6)+
  scale_linetype_manual(values = c("solid", "dashed", "dotdash", "dotdash"), 
                        labels = c("NPP", "JPP1", "JPP2", "JPP3")) +
  scale_color_manual(values = c("red", "blue", "forestgreen", "brown"), 
                     labels = c("NPP", "JPP1", "JPP2", "JPP3")) +
  theme_bw()+theme(legend.title=element_blank(), 
                   axis.text=element_text(size=9),
                   axis.title=element_text(size=8))+
  labs(x=expression(delta), y= 'Density')+
  scale_y_continuous(breaks=c(0, 3, 6, 9, 12),limits = c(0, 14))+
  scale_x_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0, 1.02))+
  theme(legend.key.width = unit(2.5, "line"), legend.position="bottom")+
  guides(col = guide_legend(nrow=2.5), linetype = guide_legend(nrow = 1))



FDDelta <- ggplot(Delta_D_df, aes(x = x, y = den, col = variable, linetype = variable)) +
  geom_line(lwd = 0.6)+
  scale_linetype_manual(values = c("solid", "dashed", "dotdash", "dotdash"), 
                        labels = c("NPP", "JPP1", "JPP2", "JPP3")) +
  scale_color_manual(values = c("red", "blue", "forestgreen", "brown"), 
                     labels = c("NPP", "JPP1", "JPP2", "JPP3")) +
  theme_bw()+theme(legend.title=element_blank(), 
                   axis.text=element_text(size=9),
                   axis.title=element_text(size=8))+
  labs(x=expression(delta), y= 'Density')+
  scale_y_continuous(breaks=c(0, 1, 2, 3, 4),limits = c(0, 4.2))+
  scale_x_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0, 1.02))+
  theme(legend.key.width = unit(2.5, "line"), legend.position="bottom")+
  guides(col = guide_legend(nrow=2.5), linetype = guide_legend(nrow = 1))

PH_delta_4plots = ggarrange(FADelta, FBDelta, FCDelta, FDDelta, ncol=4, nrow=1, common.legend = TRUE, legend="bottom")






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














#################################################################################################
#################################################################################################
#### Example 3: Normalized Power Prior in Diagnostic Tests for Medical Device                ####
#### Multinomial Endpoints (4-cells in total)                                                ####
#################################################################################################
#################################################################################################

## generate n random deviates from the Dirichlet function with shape
## parameters alpha is a vector of shapes
## source: gtools
r_dir<-function(n,alpha){
  l<-length(alpha);
  x<-matrix(rgamma(l*n,alpha),ncol=l,byrow=TRUE);
  sm<-x%*%rep(1,l);
  x/as.vector(sm);
}

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
#### Original Joint Power Prior Posterior Sampling Function for Multinomial Data
#### Use likelihood without constant (N!/(n1!...nk!))
###########################################################################################
MultinomialJPP_MCMC <- function(Data.Cur = c(10, 10, 10), Data.Hist = c(10, 10, 10),  
                                CompStat = list(n0 = NULL, n1 = NULL), 
                                prior = list(theta.dir.alpha = c(0.5, 0.5, 0.5), 
                                             delta.alpha = 1, delta.beta = 1), 
                                MCMCmethod = 'IND', rw.logit.delta = 0.1, 
                                ind.delta.alpha = 1, ind.delta.beta = 1, nsample = 5000, 
                                control.mcmc = list(delta.ini = NULL, burnin = 0, thin = 1))
{
  if(missing(CompStat)){
    if(length(Data.Cur) != length(Data.Cur)) stop("Number of Categories Must Be Equal!")
    k <- length(Data.Cur)
    n0 <- Data.Hist
    n1 <- Data.Cur
  }else{
    if(length(CompStat$n0) != length(CompStat$n1)) stop("Number of Categories Must Be Equal!")
    k <- length(CompStat$n1)
    n0 <- CompStat$n0
    n1 <- CompStat$n1
  }
  
  n0sum <- sum(n0); n1sum <- sum(n1)
  alphasum <- sum(prior$theta.dir.alpha)
  
  #### Normalized Power Prior for Dirichlet Log Marginal Posterior (unnormalized) of Delta
  LogPostDirDelta <- function(x){
    s1 <- sum(lgamma(n0*x+n1+prior$theta.dir.alpha))
    result <- ((prior$delta.alpha-1)*log(x)+(prior$delta.beta-1)*log(1-x)+
                 x*(lgamma(n0sum+1) - sum(lgamma(n0+1)))+
                 s1-lgamma(alphasum+x*n0sum+n1sum))
    return(result)
  }
  
  if(is.null(control.mcmc$delta.ini)) delta.ini = 0.5
  delta_cur <- delta.ini
  delta <- rep(delta.ini, nsample)
  theta <- matrix(0, ncol = k, nrow = nsample)
  counter <- 0
  niter <- nsample*control.mcmc$thin + control.mcmc$burnin
  
  for (i in 1:niter){
    ### Update delta with RW MH for Logit delta
    if(MCMCmethod == 'RW'){
      lgdelta_cur <- log(delta_cur/(1-delta_cur))
      lgdelta_prop <- rnorm(1, mean = lgdelta_cur, sd = sqrt(rw.logit.delta))
      delta_prop <- exp(lgdelta_prop)/(1+exp(lgdelta_prop))
      llik.prop <- LogPostDirDelta(delta_prop)
      llik.cur <- LogPostDirDelta(delta_cur)
      logr <- min(0, (llik.prop-llik.cur+log(delta_prop)+log(1-delta_prop)-log(delta_cur)-log(1-delta_cur)))
    }
    
    if(MCMCmethod == 'IND'){
      delta_prop <- rbeta(1, shape1 = ind.delta.alpha, shape2 = ind.delta.beta)
      llik.prop <- LogPostDirDelta(delta_prop)
      llik.cur <- LogPostDirDelta(delta_cur)
      logr <- min(0, (llik.prop-llik.cur+
                        dbeta(delta_cur, shape1 = ind.delta.alpha, shape2 = ind.delta.beta, log = TRUE) -
                        dbeta(delta_prop, shape1 = ind.delta.alpha, shape2 = ind.delta.beta, log = TRUE)))
    }
    
    if(runif(1) <= exp(logr)){
      delta_cur = delta_prop; counter = counter+1
    }
    if( i > control.mcmc$burnin & (i-control.mcmc$burnin)%%control.mcmc$thin==0)
    {
      delta[(i-control.mcmc$burnin)/control.mcmc$thin] <- delta_cur
      theta[(i-control.mcmc$burnin)/control.mcmc$thin, ] <- r_dir(1, alpha = n0*delta_cur+n1+prior$theta.dir.alpha)  
    }
  }
  # Calculate DIC
  meantheta <- apply(theta, 2, mean)
  D <- -2*( log(theta)%*%as.vector(n1))
  Dpbar <- -2*(sum(n1*log(meantheta)))
  DIC <- 2*mean(D)-Dpbar
  return(list(theta = theta, delta = delta, acceptrate = counter/niter, DIC = DIC))
}



###########################################################################################
#### Normalized Power Prior Posterior Sampling Function for Multinomial Data
###########################################################################################
MultinomialNPP_MCMC <- function(Data.Cur = c(10, 10, 10), Data.Hist = c(10, 10, 10),  
                                CompStat = list(n0 = NULL, n1 = NULL), 
                                prior = list(theta.dir.alpha = c(0.5, 0.5, 0.5), 
                                             delta.alpha = 1, delta.beta = 1), 
                                MCMCmethod = 'IND', rw.logit.delta = 0.1, 
                                ind.delta.alpha = 1, ind.delta.beta = 1, nsample = 5000, 
                                control.mcmc = list(delta.ini = NULL, burnin = 0, thin = 1))
{
  if(missing(CompStat)){
    if(length(Data.Cur) != length(Data.Cur)) stop("Number of Categories Must Be Equal!")
    k <- length(Data.Cur)
    n0 <- Data.Hist
    n1 <- Data.Cur
  }else{
    if(length(CompStat$n0) != length(CompStat$n1)) stop("Number of Categories Must Be Equal!")
    k <- length(CompStat$n1)
    n0 <- CompStat$n0
    n1 <- CompStat$n1
  }
  n0sum <- sum(n0); n1sum <- sum(n1)
  alphasum <- sum(prior$theta.dir.alpha)
  
  LogPostDirDelta <- function(x){
    s1 <- sum(lgamma(n0*x+n1+prior$theta.dir.alpha))
    s2 <- sum(lgamma(n0*x+prior$theta.dir.alpha))
    result <- ((prior$delta.alpha-1)*log(x)+(prior$delta.beta-1)*log(1-x)+
                 lgamma(alphasum+x*n0sum)+s1-s2-
                 lgamma(alphasum+x*n0sum+n1sum))
    return(result)
  }
  
  if(is.null(control.mcmc$delta.ini)) delta.ini = 0.5
  delta_cur <- delta.ini
  delta <- rep(delta.ini, nsample)
  theta <- matrix(0, ncol = k, nrow = nsample)
  counter <- 0
  niter <- nsample*control.mcmc$thin + control.mcmc$burnin
  
  for (i in 1:niter){
    ### Update delta with RW MH for Logit delta
    if(MCMCmethod == 'RW'){
      lgdelta_cur <- log(delta_cur/(1-delta_cur))
      lgdelta_prop <- rnorm(1, mean = lgdelta_cur, sd = sqrt(rw.logit.delta))
      delta_prop <- exp(lgdelta_prop)/(1+exp(lgdelta_prop))
      llik.prop <- LogPostDirDelta(delta_prop)
      llik.cur <- LogPostDirDelta(delta_cur)
      logr <- min(0, (llik.prop-llik.cur+log(delta_prop)+log(1-delta_prop)-log(delta_cur)-log(1-delta_cur)))
    }
    
    if(MCMCmethod == 'IND'){
      delta_prop <- rbeta(1, shape1 = ind.delta.alpha, shape2 = ind.delta.beta)
      llik.prop <- LogPostDirDelta(delta_prop)
      llik.cur <- LogPostDirDelta(delta_cur)
      logr <- min(0, (llik.prop-llik.cur+
                        dbeta(delta_cur, shape1 = ind.delta.alpha, shape2 = ind.delta.beta, log = TRUE) -
                        dbeta(delta_prop, shape1 = ind.delta.alpha, shape2 = ind.delta.beta, log = TRUE)))
    }
    if(runif(1) <= exp(logr)){
      delta_cur = delta_prop; counter = counter+1
    }
    if( i > control.mcmc$burnin & (i-control.mcmc$burnin)%%control.mcmc$thin==0)
    {
      delta[(i-control.mcmc$burnin)/control.mcmc$thin] <- delta_cur
      theta[(i-control.mcmc$burnin)/control.mcmc$thin, ] <- r_dir(1, alpha = n0*delta_cur+n1+prior$theta.dir.alpha)  
    }
  }
  # Calculate DIC
  meantheta <- apply(theta, 2, mean)
  D <- -2*( log(theta)%*%as.vector(n1))
  Dpbar <- -2*(sum(n1*log(meantheta)))
  DIC <- 2*mean(D)-Dpbar
  return(list(theta = theta, delta = delta, acceptrate = counter/niter, DIC = DIC))
}




###########################################################################################
#### Original Joint Power Prior, Log-Posterior Densities on a Grid for Multinomial Data
###########################################################################################
LPDeltaMultinomialJPP <- function(Data.Cur, Data.Hist, npoints = 1000, 
                                  CompStat = list(n0 = NULL, n1 = NULL), 
                                  prior = list(theta.dir.alpha = c(0.5, 0.5, 0.5), 
                                               delta.alpha = 1, delta.beta = 1))
{
  if(missing(CompStat)){
    if(length(Data.Cur) != length(Data.Cur)) stop("Number of Categories Must Be Equal!")
    k <- length(Data.Cur)
    n0 <- Data.Hist
    n1 <- Data.Cur
  }else{
    if(length(CompStat$n0) != length(CompStat$n1)) stop("Number of Categories Must Be Equal!")
    k <- length(CompStat$n1)
    n0 <- CompStat$n0
    n1 <- CompStat$n1
  }
  n0sum <- sum(n0); n1sum <- sum(n1)
  alphasum <- sum(prior$theta.dir.alpha)
  x <- seq(.Machine$double.eps, 1-.Machine$double.eps, length = npoints)
  
  s1 <- c()
  for(i in 1:length(x)){
    s1[i] <- sum(lgamma(n0*x[i]+n1+prior$theta.dir.alpha))
  }
  propDen <- ((prior$delta.alpha-1)*log(x)+(prior$delta.beta-1)*log(1-x)+
                x*(lgamma(n0sum+1) - sum(lgamma(n0+1)))+
                s1-lgamma(alphasum+x*n0sum+n1sum))
  return(data.frame(x = x, logden = propDen))
}

###########################################################################################
#### Normalized Power Prior, Log-Posterior Densities on a Grid for Multinomial Data
###########################################################################################
LPDeltaMultinomialNPP <- function(Data.Cur, Data.Hist, npoints = 1000, 
                                  CompStat = list(n0 = NULL, n1 = NULL), 
                                  prior = list(theta.dir.alpha = c(0.5, 0.5, 0.5), 
                                               delta.alpha = 1, delta.beta = 1))
{
  if(missing(CompStat)){
    if(length(Data.Cur) != length(Data.Cur)) stop("Number of Categories Must Be Equal!")
    k <- length(Data.Cur)
    n0 <- Data.Hist
    n1 <- Data.Cur
  }else{
    if(length(CompStat$n0) != length(CompStat$n1)) stop("Number of Categories Must Be Equal!")
    k <- length(CompStat$n1)
    n0 <- CompStat$n0
    n1 <- CompStat$n1
  }
  n0sum <- sum(n0); n1sum <- sum(n1)
  alphasum <- sum(prior$theta.dir.alpha)
  x <- seq(.Machine$double.eps, 1-.Machine$double.eps, length = npoints)
  s1 <- c(); s2 <- c()
  for(i in 1:length(x)){
    s1[i] <- sum(lgamma(n0*x[i]+n1+prior$theta.dir.alpha))
    s2[i] <- sum(lgamma(n0*x[i]+prior$theta.dir.alpha))
  }
  propDen <- ((prior$delta.alpha-1)*log(x)+(prior$delta.beta-1)*log(1-x)+
                lgamma(alphasum+x*n0sum)+s1-s2-
                lgamma(alphasum+x*n0sum+n1sum))
  return(data.frame(x = x, logden = propDen))
}




###########################################################################################
#### Start the model fitting ####
###########################################################################################

n0 <- c(9, 20, 9, 473)
n1 <- c(3, 11, 3, 669)


#### Using Joint Power Prior 
set.seed(4321)
Test_JPP <- MultinomialJPP_MCMC(Data.Cur = n1, Data.Hist = n0,
                                prior = list(theta.dir.alpha = c(0.5,0.5,0.5,0.5), delta.alpha = 1, delta.beta = 1), 
                                MCMCmethod = 'IND', ind.delta.alpha = 1, ind.delta.beta = 1, nsample = 20000, 
                                control.mcmc = list(burnin = 2000, thin = 5))

#### Using Normalized Power Prior 
set.seed(4321)
Test_NPP <- MultinomialNPP_MCMC(Data.Cur = n1, Data.Hist = n0,
                                prior = list(theta.dir.alpha = c(0.5,0.5,0.5,0.5), delta.alpha = 1, delta.beta = 1), 
                                MCMCmethod = 'IND', rw.logit.delta = 0.1, 
                                ind.delta.alpha = 1, ind.delta.beta = 1, nsample = 20000, 
                                control.mcmc = list(burnin = 2000, thin = 5))


#### No Borrowing
set.seed(123)
Test_zero<-r_dir(10000, alpha = n1+c(0.5,0.5,0.5,0.5))

#### Full borrowing
set.seed(101)
Test_full<-r_dir(10000, alpha = n0+n1+c(0.5,0.5,0.5,0.5))



#### Diagnostics
plot(Test_JPP$delta, type = 'l')
plot(Test_JPP$theta[,1], type = 'l')
plot(Test_JPP$theta[,2], type = 'l')
plot(Test_JPP$theta[,3], type = 'l')
plot(Test_JPP$theta[,4], type = 'l')
acf(Test_JPP$delta)

plot(Test_NPP$delta, type = 'l')
plot(Test_NPP$theta[,1], type = 'l')
plot(Test_NPP$theta[,2], type = 'l')
plot(Test_NPP$theta[,3], type = 'l')
plot(Test_NPP$theta[,4], type = 'l')
acf(Test_NPP$delta)


###########################################################################################
#### Model Fitting Results ####
###########################################################################################
#### Zero 
round(100*mean(Test_zero[,1]/(Test_zero[,1]+Test_zero[,3])),2)  #  50.04
round(100*emp.hpd(Test_zero[,1]/(Test_zero[,1]+Test_zero[,3]), conf = 0.95),2) # (16.67, 82.80)

round(100*mean(Test_zero[,4]/(Test_zero[,2]+Test_zero[,4])),2)  #  98.31
round(100*emp.hpd(Test_zero[,4]/(Test_zero[,2]+Test_zero[,4]), conf = 0.95),2) # (97.32, 99.22)

#### Full
round(100*mean(Test_full[,1]/(Test_full[,1]+Test_full[,3])),2)  #  49.85
round(100*emp.hpd(Test_full[,1]/(Test_full[,1]+Test_full[,3]), conf = 0.95),2) # (31.40, 68.70)

round(100*mean(Test_full[,4]/(Test_full[,2]+Test_full[,4])),2)  #  97.32
round(100*emp.hpd(Test_full[,4]/(Test_full[,2]+Test_full[,4]), conf = 0.95),2) # (96.38, 98.17)

#### JPP
round(100*mean(Test_JPP$theta[,1]/(Test_JPP$theta[,1]+Test_JPP$theta[,3])),2) # 49.98
round(100*emp.hpd(Test_JPP$theta[,1]/(Test_JPP$theta[,1]+Test_JPP$theta[,3]), conf = 0.95),2) # (18.94, 83.05)

round(100*mean(Test_JPP$theta[,4]/(Test_JPP$theta[,2]+Test_JPP$theta[,4])),2) #  98.24
round(100*emp.hpd(Test_JPP$theta[,4]/(Test_JPP$theta[,2]+Test_JPP$theta[,4]), conf = 0.95),2) # (97.27, 99.18)

round(mean(Test_JPP$delta),3)  #### 0.044
JPPpostLik <- LPDeltaMultinomialJPP(Data.Cur = c(3,11,3,669), Data.Hist = c(9,20,9,473), npoints = 5000, 
                                    prior = list(theta.dir.alpha = c(0.5,0.5,0.5,0.5), delta.alpha = 1, delta.beta = 1))
round(JPPpostLik[which.max(JPPpostLik$logden), ][1], 3) #### 0

#### NPP
round(100*mean(Test_NPP$theta[,1]/(Test_NPP$theta[,1]+Test_NPP$theta[,3])),2) # 49.88
round(100*emp.hpd(Test_NPP$theta[,1]/(Test_NPP$theta[,1]+Test_NPP$theta[,3]), conf = 0.95),2) ## (21.60, 78.84)

round(100*mean(Test_NPP$theta[,4]/(Test_NPP$theta[,2]+Test_NPP$theta[,4])),2) # 98.02
round(100*emp.hpd(Test_NPP$theta[,4]/(Test_NPP$theta[,2]+Test_NPP$theta[,4]), conf = 0.95),2) ## (96.93, 99.00)


round(mean(Test_NPP$delta),3)  #### 0.216
NPPpostLik <- LPDeltaMultinomialNPP(Data.Cur = c(3,11,3,669), Data.Hist = c(9,20,9,473), npoints = 5000, 
                                    prior = list(theta.dir.alpha = c(0.5,0.5,0.5,0.5), delta.alpha = 1, delta.beta = 1))
round(NPPpostLik[which.max(NPPpostLik$logden), ][1],3) #### 0.085


