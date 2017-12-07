calculateScb_hall_wellner <- function(estimator)
{
  return(km.ci(estimator, conf.level = 1 - confidencelevel, method = "hall-wellner"))
}

calculateScb_akritas <- function(estimator)
{
  return(estimator)
}

calculateScb_nairs_equal_precision <- function(estimator)
{
  return(km.ci(estimator, conf.level = 1 - confidencelevel, method = "epband"))
}

calculateScb_proposed_I <- function(estimator)
{
  return(estimator)
}

calculateScb_proposed_II <- function(estimator)
{
  return(estimator)
}

calculateScb_new <- function(estimator)
{
  return(estimator)
}

calculateScb_transformed_hall_wellner <- function(estimator)
{
  return(km.ci(estimator, conf.level = 1 - confidencelevel, method = "loghall"))
}

calculateScb_transformed_akritas <- function(estimator)
{
  return(estimator)
}

calculateScb_transformed_nairs_equal_precision <- function(estimator)
{
  return(km.ci(estimator, conf.level = 1 - confidencelevel, method = "logep"))
}

calculateScb_proposed_III <- function(estimator)
{
  return(estimator)
}

calculateScb_proposed_IV <- function(estimator)
{
  return(estimator)
}

calculateScb_transformed_new <- function(estimator)
{
  return(estimator)
}