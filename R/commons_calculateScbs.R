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

calculateScb_proposed_I <- function(estimator, sample, modelFunction, mleTheta, parameterLimits)
{
  W.stars = replicate(bootstrapRepetitions, getWFromOneBootstrapSample(estimator, sample, modelFunction, mleTheta, parameterLimits, 
                                                                           estimators.dikta_2))
  W.stars = sort(W.stars)
  
  q.alpha = W.stars[floor(bootstrapRepetitions * (1 - confidencelevel))]
  
  upper = estimator$surv + 1 / sqrt(samplesize) * q.alpha
  lower = estimator$surv - 1 / sqrt(samplesize) * q.alpha
  scb = list(time = estimator$time, surv = estimator$surv,
                 upper = upper, lower = lower, estimator = estimator)
  
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

calculateScb_new_proposed_I <- function(estimator, sample, modelFunction, mleTheta, parameterLimits)
{
  W.stars = replicate(bootstrapRepetitions, getWFromOneBootstrapSample(estimator, sample, modelFunction, mleTheta, parameterLimits, 
                                                                 estimators.dikta_3))
  W.stars = sort(W.stars)
  
  q.alpha = W.stars[floor(bootstrapRepetitions * (1 - confidencelevel))]

  upper = estimator$surv + 1 / sqrt(samplesize) * q.alpha
  lower = estimator$surv - 1 / sqrt(samplesize) * q.alpha
  scb = list(time = estimator$time, surv = estimator$surv,
               upper = upper, lower = lower, estimator = estimator)
  
  return(scb)
}

calculateScb_new_proposed_II <- function(estimator)
{
  scb = estimator
  scb$lower = scb$surv - 0.15
  scb$upper = scb$surv + 0.15
  scb$estimator = estimator
  return(scb)
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

calculateScb_proposed_III <- function(estimator, sample, modelFunction, mleTheta, parameterLimits)
{
  W.stars = replicate(bootstrapRepetitions, getWFromOneBootstrapSample(estimator, sample, modelFunction, mleTheta, parameterLimits, 
                                                                       estimators.dikta_2))
  W.stars = sort(W.stars)
  
  q.alpha = W.stars[floor(bootstrapRepetitions * (1 - confidencelevel))]
  
  exponent = 1/sqrt(samplesize) * q.alpha / (estimator$surv * log(estimator$surv))
  upper = estimator$surv ^ (exp(-exponent))
  lower = estimator$surv ^ (exp(exponent))
  scb = list(time = estimator$time, surv = estimator$surv,
                 upper = upper, lower = lower, estimator = estimator)
  
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

calculateScb_new_proposed_III <- function(estimator, sample, modelFunction, mleTheta, parameterLimits)
{
  W.stars = replicate(bootstrapRepetitions, getWFromOneBootstrapSample(estimator, sample, modelFunction, mleTheta, parameterLimits, 
                                                                       estimators.dikta_3))
  W.stars = sort(W.stars)
  
  q.alpha = W.stars[floor(bootstrapRepetitions * (1 - confidencelevel))]
  
  exponent = 1/sqrt(samplesize) * q.alpha / (estimator$surv * log(estimator$surv))
  upper = estimator$surv ^ (exp(-exponent))
  lower = estimator$surv ^ (exp(exponent))
  scb = list(time = estimator$time, surv = estimator$surv,
                 upper = upper, lower = lower, estimator = estimator)
  
  return(scb)
}

calculateScb_new_proposed_IV <- function(estimator)
{
  scb = estimator
  scb$lower = scb$surv - 0.15
  scb$upper = scb$surv + 0.15
  scb$estimator = estimator
  return(scb)
}


getWFromOneBootstrapSample <- function(estimator, sample, modelFunction, mleTheta, parameterLimits, estimatorFunction)
{
  bootstrapSample = generateBootstrapSample(sample, modelFunction, mleTheta)
  bootstrapMleTheta = estimators.calculateMaximumLikelihoodEstimator(bootstrapSample, modelFunction, parameterLimits)
  
  bootstrapEstimator = estimatorFunction(bootstrapSample, modelFunction, bootstrapMleTheta)
  bootstrapEstimatorAtOriginalTimes = sapply(estimator$time, estimators.getEstimatorFromT, estimator = bootstrapEstimator)
  
  W = sqrt(samplesize) * max(abs(bootstrapEstimatorAtOriginalTimes - estimator$surv))
  
  return(W)
}