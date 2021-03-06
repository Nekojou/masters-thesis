
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
  R[R == 0] = 1e-250
  R[R == 1] = 1-1e-10
  loglikelihood = 1/length(sample$Z) * 
            sum(sample$delta * log(R) + 
            (1 - sample$delta) * log(1 - R))
  if(!is.finite(loglikelihood))
  {
    print("here")
  }
  return(-loglikelihood)
}

estimators.calculateMaximumLikelihoodEstimator <- function(sample, modelFunction, parameterLimits)
{
  mle = optim(par = c(1,1), fn = estimators.likelihoodFunction,
              sample = sample, modelFunction = modelFunction,
              lower = parameterLimits[["lower"]], upper = parameterLimits[["upper"]], method = "L-BFGS-B")
  return(mle$par)
}

# function to return the estimators value at time t
estimators.getEstimatorFromT <- function(t, estimator)
{
  return(min(1, estimator$surv[estimator$time <= t]))
}