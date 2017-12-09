
estimators.kaplan_meyer <- function(sample)
{
  return(survfit(Surv(sample$Z, sample$delta)~1))
}

estimators.dikta_2 <- function(sample, modelFunction, mleTheta)
{
  surv = c()
  temp = 1 - modelFunction(sample$Z, mleTheta) / (length(sample$Z) - rank(sample$Z) + 1)
  orderedZ = sort(sample$Z)
  for(i in 1:length(sample$Z)){
    newvalue=prod(temp[sample$Z<=orderedZ[i]])
    surv=c(surv,newvalue)
  }
  return(list(time=orderedZ,surv=surv))
}

estimators.dikta_3 <- function(sample, modelFunction, mleTheta)
{
  surv = c()
  temp = 1 - modelFunction(sample$Z, mleTheta) / (length(sample$Z) - rank(sample$Z) + modelFunction(sample$Z, mleTheta))
  orderedZ = sort(sample$Z)
  for(i in 1:length(sample$Z)){
    newvalue=prod(temp[sample$Z<=orderedZ[i]])
    surv=c(surv,newvalue)
  }
  return(list(time=orderedZ,surv=surv))
}

estimators.likelihoodFunction <- function(sample, theta, modelFunction)
{
  R = sapply(sample$Z, modelFunction, theta = theta)
  # R[!is.finite(log(R))] = 0.0001
  # R[!is.finite(log(1-R))] = 0.9999
  #   print(sample)
  #   print(theta)
  #   print(-1/length(sample$Z) * sum(data$sample*log(R)+(1-sample$delta)*log(1-R)))
  result = -1/length(sample$Z) * 
            sum(sample$delta * log(R) + 
            (1 - sample$delta) * log(1 - R))
  if(!is.finite(result))
  {
    print("here")
  }
  return(result)
}

estimators.calculateMaximumLikelihoodEstimator <- function(sample, modelFunction)
{
  mle = optim(par = c(1,1), fn = estimators.likelihoodFunction,
              sample = sample, modelFunction = modelFunction,
              lower = c(0.0001, -Inf), upper = c(Inf, Inf), method = "L-BFGS-B")
  return(mle$par)
}