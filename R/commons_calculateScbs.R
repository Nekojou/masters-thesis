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

calculateScb_new <- function(estimator)
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