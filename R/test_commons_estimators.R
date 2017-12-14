# test function estimators.dikta_2
# As sample, the example 
# Z =  1.7, 4.9, 1.4, 8.5, 3.4, 0.6, 2.4, 3.5, 3.7, 2.0
# delta = 1, 1, 1, 1, 1, 1, 1, 0, 1, 1
# modelFunction(Z, theta) = theta
# theta = 0.9
# Second example:
# same sample, different model function
# modelFunction(Z, theta) = theta^Z
# theta = 0.9
# all results for expectedSurv were calculated by hand
test.commons.estimators.dikta_2 <- function()
{
  # example 1
  modelFunction = function(Z, theta){ return(theta) }
  theta = 0.9
  
  Z = c(1.7, 4.9, 1.4, 8.5, 3.4, 0.6, 2.4, 3.5, 3.7, 2.0)
  delta = c(1, 1, 1, 1, 1, 1, 1, 0, 1, 1)
  sample = data.frame(Z = Z, delta = delta)
  
  estimator = estimators.dikta_2(sample, modelFunction, theta)
  
  expectedTime = sort(Z)
  expectedSurv = c(0.91, 0.819, 0.7268625, 0.63340875, 0.5383974375, 0.441485899875, 
                   0.34215157153125, 0.239506100071875, 0.13172835503953125, 0.013172835503953125)
  stopifnot(all(floatCompare(expectedTime, estimator$time)))
  stopifnot(all(floatCompare(expectedSurv, estimator$surv)))
  
  # example 2
  modelFunction = function(Z, theta){ return(theta^Z) }
  theta = 0.9
  
  estimator = estimators.dikta_2(sample, modelFunction, theta)
  
  expectedSurv = c(0.9061259606, 0.8192528091, 0.7336397100, 0.6487471150, 0.5647805838,
                   0.4858338284, 0.4018343590, 0.3111310567, 0.2182982108, 0.1291502309)
  stopifnot(all(floatCompare(expectedTime, estimator$time)))
  stopifnot(all(floatCompare(expectedSurv, estimator$surv)))
}

# test function estimators.dikta_3
# First example:
# Z =  1.7, 4.9, 1.4, 8.5, 3.4, 0.6, 2.4, 3.5, 3.7, 2.0
# delta = 1, 1, 1, 1, 1, 1, 1, 0, 1, 1
# modelFunction(Z, theta) = theta
# theta = 0.9
# Second example:
# same sample, different model function
# modelFunction(Z, theta) = theta^Z
# theta = 0.9
# all results for expectedSurv were calculated by hand
test.commons.estimators.dikta_3 <- function()
{
  # example 1
  modelFunction = function(Z, theta){ return(theta) }
  theta = 0.9
  
  Z = c(1.7, 4.9, 1.4, 8.5, 3.4, 0.6, 2.4, 3.5, 3.7, 2.0)
  delta = c(1, 1, 1, 1, 1, 1, 1, 0, 1, 1)
  sample = data.frame(Z = Z, delta = delta)
  
  estimator = estimators.dikta_3(sample, modelFunction, theta)
  
  expectedTime = sort(Z)
  expectedSurv = c(0.9090909090, 0.8171603677, 0.7240661486, 0.6296227379, 0.5335785914,
                   0.4355743603, 0.3350572002, 0.2310739312, 0.1216178585, 0)
  stopifnot(all(floatCompare(expectedTime, estimator$time)))
  stopifnot(all(floatCompare(expectedSurv, estimator$surv)))
  
  # example 2
  modelFunction = function(Z, theta){ return(theta^Z) }
  theta = 0.9
  
  estimator = estimators.dikta_3(sample, modelFunction, theta)
  
  expectedSurv = c(0.9055473474, 0.8173862744, 0.7301806411, 0.6433309613, 0.5568448749,
                   0.4740199321, 0.3852160582, 0.2877786271, 0.1802283702, 0)
  stopifnot(all(floatCompare(expectedTime, estimator$time)))
  stopifnot(all(floatCompare(expectedSurv, estimator$surv)))
}

test.commons.estimators.likelihoodFunction <- function()
{
  estimators.likelihoodFunction(sample, theta, modelFunction)
}

test.commons.estimators.calculateMaximumLikelihoodEstimator <- function()
{
  set.seed(22431873)
  alpha1 = 2
  beta1 = 3
  alpha2 = 3.5
  beta2 = 4.5
  X = rweibull(100, alpha1, beta1)
  Y = rweibull(100, alpha2, beta2)
  
  Z = mapply(min, X, Y)
  delta = ifelse(X<=Y,1,0)
  sample = data.frame(Z,delta)
  
  modelFunction = function(Z, theta){ return(theta[1]/(theta[1]+Z^theta[2])) }
  theta1 = (alpha1*beta1^(-alpha1))/(alpha2*beta2^(-alpha2))
  theta2 = alpha2-alpha1
  theta = c(theta1, theta2)
  
  mleTheta = estimators.calculateMaximumLikelihoodEstimator(sample, modelFunction)
  
  print(theta)
  print(mleTheta)
  
  stopifnot(mleTheta == c(8.387158, 1.340541))
}
