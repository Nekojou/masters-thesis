setwd(dirname(parent.frame(2)$ofile))
source("commons.R")

# set fixed distribution parameters
study1.alpha1 = 2
study1.beta1 = 3
study1.beta2 = 4.5

# generate variable parameters alpha2
study1.numberOfCases = 10
study1.alpha2List = c(seq(1.1,5.5,length.out=study1.numberOfCases))

# Generates random sample of (Z,delta) from X,Y~weibull 
# depending on input parameter alpha2
study1.generateRandomSample <- function(alpha2, samplesize=100)
{
  X = rweibull(samplesize, study1.alpha1, study1.beta1)
  Y = rweibull(samplesize, alpha2, study1.beta2)
  
  Z = mapply(min, X, Y)
  delta = ifelse(X<=Y,1,0)
  return(data.frame(Z,delta))
}

# true Survival function in this study
study1.survivalfunction <- function(t,alpha=study1.alpha1,beta=study1.beta1)
{
  return(1-pweibull(t,alpha,beta))
}

# true model function for SRCM in this study
study1.modelfunction <- function(z,theta)
{
  return(theta[1]/(theta[1]+z^theta[2]))
}

# calculate parameters theta = (theta1,theta2) from this studys parameters
# depending on the variable parameter alpha2
study1.calculateNewParametrizationTheta <- function(alpha2)
{
  theta1 = (study1.alpha1*study1.beta1^(-study1.alpha1))/(alpha2*study1.beta2^(-alpha2))
  theta2 = alpha2-study1.alpha1
  return(c(theta1, theta2))
}

study1.run <- function()
{
  
  # generate all samples
  allSamples = generateAllRandomSamples(study1.generateRandomSample, study1.alpha2List)
  
  # calculate global time interval 
  globalTimeInterval = calculateGlobalTimeInterval(allSamples)
  
  # set up results list
  resultsList = setUpResultsVariables()

  runAllStudyCases(allSamples, study1.alpha2List)
  
}