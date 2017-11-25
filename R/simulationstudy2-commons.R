alpha2=-5.92
nCases=16
var.alpha=c(seq(3,6,length.out=nCases))

# Generates random sample of (Z,delta) from 
# z~UNI(0,1) and delta~BERN
# depending on input parameter alpha1
generateRandomSample<-function(alpha1){
  Z = runif(samplesize,0,1)
  delta = sapply(lapply(Z, true.m, theta=c(alpha1,alpha2)), rbinom, n=1, size=1)
  return(data.frame(Z,delta))
}

# True Survival function in this case
S<-function(z){
  return(1-punif(z,0,1))
}

# True model for SRCM in this case (two-parameter complementary log-log model)
true.m<-function(z,theta){
  return(1-exp(-exp(theta[1]+theta[2]*z)))
}

# Misspecified model for SRCM in this case (one-parameter complementary log-log model with alpha1 = 4)
m<-function(z,theta){
  return(1-exp(-exp(4+theta*z)))
}

# Log likelihood function for m 
# data = sample with Z and delta, par = theta
m.LL<-function(data,par){
  R = m(data$Z,par)
  R[!is.finite(log(R))] = 0.0001
  R[!is.finite(log(1-R))] = 0.9999
  return(-1/length(data$Z) * sum(data$delta*log(R)+(1-data$delta)*log(1-R)))
}

calc.mle<-function(sample){
  return(optim(par=1,fn=m.LL,data=sample,method="L-BFGS-B"))  
}