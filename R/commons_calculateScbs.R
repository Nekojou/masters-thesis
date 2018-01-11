calculateScb_hall_wellner <- function(estimator)
{
  scb = km.ci(estimator, conf.level = 1 - confidencelevel, method = "hall-wellner")
  scb$estimator = estimator
  return(scb)
}

calculateScb_akritas <- function(estimator)
{
  scb = estimator
  scb$lower = scb$surv - 0.15
  scb$upper = scb$surv + 0.15
  scb$estimator = estimator
  return(scb)
}

calculateScb_nairs_equal_precision <- function(estimator)
{
  scb = km.ci(estimator, conf.level = 1 - confidencelevel, method = "epband")
  scb$estimator = estimator
  return(scb)
}

calculateScb_proposed_I <- function(estimator)
{
  scb = estimator
  scb$lower = scb$surv - 0.15
  scb$upper = scb$surv + 0.15
  scb$estimator = estimator
  return(scb)
}

calculateScb_proposed_II <- function(estimator)
{
  scb = estimator
  scb$lower = scb$surv - 0.15
  scb$upper = scb$surv + 0.15
  scb$estimator = estimator
  return(scb)
}

calculateScb_new <- function(estimator, sample, modelFunction, mleTheta, parameterLimits)
{
  W.hat.new = c()
  
  for(bootIterator in 1:bootstrapRepetitions){
    bootstrapSample = generateBootstrapSample(sample, modelFunction, mleTheta)
    bootstrapMleTheta = estimators.calculateMaximumLikelihoodEstimator(bootstrapSample, modelFunction, parameterLimits)
    
    # TODO: check if this calculation is done correctly
    bootstrapEstimator = estimators.dikta_3(bootstrapSample, modelFunction, bootstrapMleTheta)
    bootstrapEstimatorAtOriginalTimes = sapply(estimator$time, estimators.getEstimatorFromT, estimator = bootstrapEstimator)
    W.hat.new = c(W.hat.new, 
                  sqrt(samplesize) * max(abs(bootstrapEstimatorAtOriginalTimes - estimator$surv)))
  }
  
  W.hat.new = sort(W.hat.new)
  q.alpha.new = W.hat.new[floor(bootstrapRepetitions * (1 - confidencelevel))]

  upper.new = estimator$surv + 1 / sqrt(samplesize) * q.alpha.new
  lower.new = estimator$surv - 1 / sqrt(samplesize) * q.alpha.new
  scb.new = list(time = estimator$time, surv = estimator$surv,
               upper = upper.new, lower = lower.new, estimator = estimator)
  
  return(scb.new)
}

calculateScb_transformed_hall_wellner <- function(estimator)
{
  scb = km.ci(estimator, conf.level = 1 - confidencelevel, method = "loghall")
  scb$estimator = estimator
  return(scb)
}

calculateScb_transformed_akritas <- function(estimator)
{
  scb = estimator
  scb$lower = scb$surv - 0.15
  scb$upper = scb$surv + 0.15
  scb$estimator = estimator
  return(scb)
}

calculateScb_transformed_nairs_equal_precision <- function(estimator)
{
  scb = km.ci(estimator, conf.level = 1 - confidencelevel, method = "logep")
  scb$estimator = estimator
  return(scb)
}

calculateScb_proposed_III <- function(estimator)
{
  scb = estimator
  scb$lower = scb$surv - 0.15
  scb$upper = scb$surv + 0.15
  scb$estimator = estimator
  return(scb)
}

calculateScb_proposed_IV <- function(estimator)
{
  scb = estimator
  scb$lower = scb$surv - 0.15
  scb$upper = scb$surv + 0.15
  scb$estimator = estimator
  return(scb)
}

calculateScb_transformed_new <- function(estimator)
{
  scb = estimator
  scb$lower = scb$surv - 0.15
  scb$upper = scb$surv + 0.15
  scb$estimator = estimator
  return(scb)
}