########################################################################################################
########################################################################################################
#### Example 1: Normalized Power Prior in Borrowing Historical pH Data in Water Quality Assessment  ####
#### Normal Endpoints (response, non-response)                                               ###########
#### Assume Joint Prior of mu and sigmasq is (1/sigmasq)^a                                   ###########
#### a = 1.5 Jeffreys; a = 1 reference                                                       ###########
########################################################################################################
########################################################################################################

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
  ## a = 1 for Reference Prior; a = 1.5 for Jeffrey's Prior
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
}

PHOPPL1 <- function(mean0, n0, sd0, mean1, n1, sd1, nsample){
  POST <- NormalOPP_MCMC_suff(CompStat = list(mean0 = mean0, n0 = n0, var0 = (sd0^2)*(1-1/n0), 
                                         mean1 = mean1, n1 = n1, var1 = (sd1^2)*(1-1/n1)), 
                         prior = list(joint.a = 1, delta.alpha = 1, delta.beta = 1), MCMCmethod = 'RW', 
                         rw.logit.delta = 0.3, ind.delta.alpha= 1, ind.delta.beta= 1, nsample = nsample, 
                         control.mcmc = list(delta.ini = NULL, burnin = 2000, thin = 5))
  POSTL <- POST$mu+qnorm(0.1)*sqrt(POST$sigmasq)
}

PHOPPL2 <- function(mean0, n0, sd0, mean1, n1, sd1, nsample){
  POST <- NormalOPP_MCMC_datalik(CompStat = list(mean0 = mean0, n0 = n0, var0 = (sd0^2)*(1-1/n0), 
                                         mean1 = mean1, n1 = n1, var1 = (sd1^2)*(1-1/n1)), 
                         prior = list(joint.a = 1, delta.alpha = 1, delta.beta = 1), MCMCmethod = 'RW', 
                         rw.logit.delta = 0.3, ind.delta.alpha= 1, ind.delta.beta= 1, nsample = nsample, 
                         control.mcmc = list(delta.ini = NULL, burnin = 2000, thin = 5))
  POSTL <- POST$mu+qnorm(0.1)*sqrt(POST$sigmasq)
}

PHOPPL3 <- function(mean0, n0, sd0, mean1, n1, sd1, nsample){
  POST <- NormalOPP_MCMC_arbitrary(CompStat = list(mean0 = mean0, n0 = n0, var0 = (sd0^2)*(1-1/n0), 
                                         mean1 = mean1, n1 = n1, var1 = (sd1^2)*(1-1/n1)), 
                         prior = list(joint.a = 1, delta.alpha = 1, delta.beta = 1), MCMCmethod = 'RW', 
                         rw.logit.delta = 0.3, ind.delta.alpha= 1, ind.delta.beta= 1, nsample = nsample, 
                         control.mcmc = list(delta.ini = NULL, burnin = 2000, thin = 5))
  POSTL <- POST$mu+qnorm(0.1)*sqrt(POST$sigmasq)
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

FunL(OPPL1_A); FunL(OPPL1_B); FunL(OPPL1_C); FunL(OPPL1_D)
plot(OPPL1_A,type ='l');acf(OPPL1_A)
plot(OPPL1_B,type ='l');acf(OPPL1_B)
plot(OPPL1_C,type ='l');acf(OPPL1_C)
plot(OPPL1_D,type ='l');acf(OPPL1_D)
#[1] 0.385 0.310
#[1] 0.051 0.300
#[1] 0.003 0.250
#[1] 0.959 0.300



#### OPP2: based on the loglikelihood of original data, including the 2pi
set.seed(1234)
OPPL2_A <- PHOPPL2(mean0 = 7.05, n0 = 62, sd0 = 0.47, mean1 = 6.91, n1 = 16, sd1 = 0.90, nsample = 20000)
OPPL2_B <- PHOPPL2(mean0 = 6.73, n0 = 31, sd0 = 0.71, mean1 = 6.78, n1 = 12, sd1 = 1.03, nsample = 20000)
OPPL2_C <- PHOPPL2(mean0 = 6.95, n0 = 84, sd0 = 0.49, mean1 = 6.43, n1 = 24, sd1 = 0.88, nsample = 20000)
OPPL2_D <- PHOPPL2(mean0 = 7.88, n0 = 75, sd0 = 0.67, mean1 = 7.87, n1 = 21, sd1 = 1.11, nsample = 20000)

FunL(OPPL2_A); FunL(OPPL2_B); FunL(OPPL2_C); FunL(OPPL2_D)
plot(OPPL2_A,type ='l');acf(OPPL2_A)
plot(OPPL2_B,type ='l');acf(OPPL2_B)
plot(OPPL2_C,type ='l');acf(OPPL2_C)
plot(OPPL2_D,type ='l');acf(OPPL2_D)
#[1] 0.201 0.320
#[1] 0.07 0.45
#[1] 0.002 0.250
#[1] 0.886 0.350



#### OPP3: with arbitrary constant exp(200)*(2pi)^(n0/2)
set.seed(123)
OPPL3_A <- PHOPPL3(mean0 = 7.05, n0 = 62, sd0 = 0.47, mean1 = 6.91, n1 = 16, sd1 = 0.90, nsample = 20000)
OPPL3_B <- PHOPPL3(mean0 = 6.73, n0 = 31, sd0 = 0.71, mean1 = 6.78, n1 = 12, sd1 = 1.03, nsample = 20000)
OPPL3_C <- PHOPPL3(mean0 = 6.95, n0 = 84, sd0 = 0.49, mean1 = 6.43, n1 = 24, sd1 = 0.88, nsample = 20000)
OPPL3_D <- PHOPPL3(mean0 = 7.88, n0 = 75, sd0 = 0.67, mean1 = 7.87, n1 = 21, sd1 = 1.11, nsample = 20000)

FunL(OPPL3_A); FunL(OPPL3_B); FunL(OPPL3_C); FunL(OPPL3_D)
plot(OPPL3_A,type ='l');acf(OPPL3_A)
plot(OPPL3_B,type ='l');acf(OPPL3_B)
plot(OPPL3_C,type ='l');acf(OPPL3_C)
plot(OPPL3_D,type ='l');acf(OPPL3_D)
#[1] 0.184 0.330
#[1] 0.068 0.470
#[1] 0.001 0.260
#[1] 0.874 0.350



#### NPP: 
set.seed(1234)
NPPL_A <- PHNPPL(mean0 = 7.05, n0 = 62, sd0 = 0.47, mean1 = 6.91, n1 = 16, sd1 = 0.90, nsample = 20000)
NPPL_B <- PHNPPL(mean0 = 6.73, n0 = 31, sd0 = 0.71, mean1 = 6.78, n1 = 12, sd1 = 1.03, nsample = 20000)
NPPL_C <- PHNPPL(mean0 = 6.95, n0 = 84, sd0 = 0.49, mean1 = 6.43, n1 = 24, sd1 = 0.88, nsample = 20000)
NPPL_D <- PHNPPL(mean0 = 7.88, n0 = 75, sd0 = 0.67, mean1 = 7.87, n1 = 21, sd1 = 1.11, nsample = 20000)

FunL(NPPL_A); FunL(NPPL_B); FunL(NPPL_C); FunL(NPPL_D)
plot(NPPL_A,type ='l');acf(NPPL_A)
plot(NPPL_B,type ='l');acf(NPPL_B)
plot(NPPL_C,type ='l');acf(NPPL_C)
plot(NPPL_D,type ='l');acf(NPPL_D)

#[1] 0.488 0.260
#[1] 0.047 0.260
#[1] 0.004 0.240
#[1] 0.986 0.250


#### Reference Prior Only, no historical data
set.seed(123)
NBL_A <- PHNBL(mean1 = 6.91, n1 = 16, sd1 = 0.90, nsample = 20000)
NBL_B <- PHNBL(mean1 = 6.78, n1 = 12, sd1 = 1.03, nsample = 20000)
NBL_C <- PHNBL(mean1 = 6.43, n1 = 24, sd1 = 0.88, nsample = 20000)
NBL_D <- PHNBL(mean1 = 7.87, n1 = 21, sd1 = 1.11, nsample = 20000)

FunL(NBL_A); FunL(NBL_B); FunL(NBL_C); FunL(NBL_D)
plot(OPPL3_A,type ='l');acf(OPPL3_A)
plot(OPPL3_B,type ='l');acf(OPPL3_B)
plot(OPPL3_C,type ='l');acf(OPPL3_C)
plot(OPPL3_D,type ='l');acf(OPPL3_D)



#### How much they can borrow?
# mean0 = 7.05; n0 = 62; sd0 = 0.47; mean1 = 6.91; n1 = 16; sd1 = 0.90; nsample = 20000
# OPP1 = NormalOPP_MCMC_suff(CompStat = list(mean0 = mean0, n0 = n0, var0 = (sd0^2)*(1-1/n0), 
#                                 mean1 = mean1, n1 = n1, var1 = (sd1^2)*(1-1/n1)), 
#                 prior = list(joint.a = 1, delta.alpha = 1, delta.beta = 1), MCMCmethod = 'RW', 
#                 rw.logit.delta = 0.3, ind.delta.alpha= 1, ind.delta.beta= 1, nsample = nsample, 
#                 control.mcmc = list(delta.ini = NULL, burnin = 2000, thin = 5))
# 
# OPP2 = NormalOPP_MCMC_datalik(CompStat = list(mean0 = mean0, n0 = n0, var0 = (sd0^2)*(1-1/n0), 
#                                        mean1 = mean1, n1 = n1, var1 = (sd1^2)*(1-1/n1)), 
#                        prior = list(joint.a = 1, delta.alpha = 1, delta.beta = 1), MCMCmethod = 'RW', 
#                        rw.logit.delta = 0.3, ind.delta.alpha= 1, ind.delta.beta= 1, nsample = nsample, 
#                        control.mcmc = list(delta.ini = NULL, burnin = 2000, thin = 5))
# 
# OPP3 = NormalOPP_MCMC_arbitrary(CompStat = list(mean0 = mean0, n0 = n0, var0 = (sd0^2)*(1-1/n0), 
#                                        mean1 = mean1, n1 = n1, var1 = (sd1^2)*(1-1/n1)), 
#                        prior = list(joint.a = 1, delta.alpha = 1, delta.beta = 1), MCMCmethod = 'RW', 
#                        rw.logit.delta = 0.3, ind.delta.alpha= 1, ind.delta.beta= 1, nsample = nsample, 
#                        control.mcmc = list(delta.ini = NULL, burnin = 2000, thin = 5))
# 
# mean(OPP1$delta) 0.1466039
# mean(OPP2$delta) 0.01620353
# mean(OPP3$delta) 0.9952083











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
















