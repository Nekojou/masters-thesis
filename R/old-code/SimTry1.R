# Obtain one Bootstrap sample of size n from Z, m and thetaD
twoStageBootstrap<-function(Z, mf, thetaD){
  n = length(Z)
  Zs = sample(Z, replace = TRUE)
  deltas = sapply(lapply(Zs, mf, theta=thetaD), rbinom, n=1, size=1)
  return(list(Zs, deltas))
}

#
m3.1<-function(x,theta){
  return(theta[1]/(theta[1]+x^theta[2]))
}

#
my.rweibull<-function(n,alpha,beta){
  return(rweibull(n,alpha,1/beta))
}

#
generateData3.1<-function(samplesize, alpha, beta){
  X = my.rweibull(samplesize, alpha[1], beta[1])
  Y = my.rweibull(samplesize, alpha[2], beta[2])
  
  Z = mapply(min, X, Y)
  delta = ifelse(X<=Y,1,0)
  return(list(Z,delta))
}
  