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

n0 <- c(9,20,9,473)
n1 <- c(3,11,3,669)


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
