alpha1=2
beta1=3
beta2=4.5

nCases=10

var.alpha=c(seq(1.1,5.5,length.out=nCases))

# Generates random sample of (Z,delta) from X,Y~weibull 
# depending on input parameter alpha2
generateRandomSample<-function(alpha2){
  X = rweibull(samplesize, alpha1, beta1)
  Y = rweibull(samplesize, alpha2, beta2)
  
  Z = mapply(min, X, Y)
  delta = ifelse(X<=Y,1,0)
  return(data.frame(Z,delta))
}

# True Survival function in this case
S<-function(z){
  return(1-pweibull(z,alpha1,beta1))
}

# True model for SRCM in this case
m<-function(z,theta){
  return(theta[1]/(theta[1]+z^theta[2]))
}

# Log likelihood function for m 
# data = sample with Z and delta, par = c(theta[1],theta[2])
m.LL<-function(data,par){
  R = m(data$Z,par)
  R[!is.finite(log(R))] = 0.0001
  R[!is.finite(log(1-R))] = 0.9999
#   print(data)
#   print(par)
#   print(-1/length(data$Z) * sum(data$delta*log(R)+(1-data$delta)*log(1-R)))
  return(-1/length(data$Z) * sum(data$delta*log(R)+(1-data$delta)*log(1-R)))
}

calc.mle<-function(sample){
  return(optim(par=c(1,1),fn=m.LL,data=sample,lower=c(0.0001,-Inf),upper=c(Inf,Inf),method="L-BFGS-B"))  
}

calculateTheta<-function(alpha2){
  theta1 = (alpha1*beta1^(-alpha1))/(alpha2*beta2^(-alpha2))
  theta2 = alpha2-alpha1
  return(c(theta1, theta2))
}
