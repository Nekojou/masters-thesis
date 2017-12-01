# used packages


# common parameters
bootstrapRepetitions=500
samplesize=100
montecarloRepetitions=1000
confidencelevel=0.05


runAllStudyCases <- function(allSamples, variableParameterList)
{
  for(parameter in variableParameterList)
  {
    runOneCaseStudy(parameter, allSamples[[parameter]])
  }
}

runOneCaseStudy <- function(caseParameter, caseSamples)
{
  print(caseParameter)
}

# TODO: Functions to be transferred:
# S.estimator2<-function(z,theta.hat){}
# S.estimator2.t<-function(t,z,theta.hat){}
# S.estimator3<-function(z,theta.hat){}
# S.estimator3.t<-function(t,z,theta.hat){}
# calculateCoverageAndEnclosedAreaAndWidth<-function(scb,globalTimeInterval){}
# save.ergs<-function(prefix){}
# load.ergs<-function(prefix){}
# plot.ergs<-function(prefix){}

# generate all samples
# arguments are:
# the function to generate one random sample
# the list of parameters for the sample generation 
# and the number of montecarlo repetitions to be generated for each parameter
generateAllRandomSamples <- function(generateRandomSample, variableParameterList, numberOfMontecarloRepetitions = montecarloRepetitions)
{
  samples <- list()
  for(parameter in variableParameterList){
    samples[[parameter]] = replicate(n = numberOfMontecarloRepetitions,
                                     generateRandomSample(parameter),
                                     simplify = FALSE)
  }
  return(samples)
}

# Obtain one bootstrap sample from original sample using the two-stage bootstrap method
# ZStar are generated through classical bootstrap from original Z
# deltaStar are generated from bernoulli variables 
# with success parameter modelFunction(ZStar_i, thetaMLE)
generateBootstrapSample <- function(originalSample, modelFunction, thetaMLE)
{
  ZStar = sample(originalSample$Z, replace = TRUE)
  deltaStar = sapply(sapply(ZStar, modelFunction, theta=thetaMLE), rbinom, n=1, size=1)
  return(data.frame(Z=ZStar, delta=deltaStar))
}


# calculates the censoring rate from a given sample
calculateCensoringRate <- function(sample){
  return(1-mean(sample$delta))
}

# utility function to calculate the max of all uncensored Z within a sample 
maxUncensoredZ <- function(sample)
{
  return(max(sample$Z[sample$delta==1]))
}
# utility function to calculate the min of all uncensored Z within a sample 
minUncensoredZ <- function(sample)
{
  return(min(sample$Z[sample$delta==1]))
}

# calculates the limits t1 and t2 of the global time intervall for all samples
# t1 is slightly greater than the maximum of all uncensored minimums per sample 
# t2 is slightly smaller than the minimum of all uncensored maximums per sample
calculateGlobalTimeInterval <- function(allSamples)
{
  deltaT = 0.0001
  
  t1 = max(sapply(allSamples, sapply, minUncensoredZ)) + deltaT
  t2 = min(sapply(allSamples, sapply, maxUncensoredZ)) - deltaT
  
  return(c(t1,t2))
}

# calculates m1 and m2 for a specific (ordered) sample Z so that 
# (Z_(m1+1),Z_(m2-1)) is subset of (t1,t2) is subset of (z_m1, z_m2)
# where t1 and t2 are the limits of the globalTimeInterval 
# (see above function calculateGlobalTimeInterval).
# The index limits m1 and m2 are used to calculate the eaea and eaw 
# consistantly for all scb types 
calculateIndexLimitsForStatistics<-function(sample,globalTimeInterval)
{
  m1 = sum(sample$Z<globalTimeInterval[1])
  m2 = sum(sample$Z<globalTimeInterval[2]) + 1
  return(c(m1,m2))
}

# generate empty datasets for results for each scb type
# define a list to store all results
# and give names to the list elements
setUpResultsList <- function()
{
  resultsList = replicate(n = 10,
                          expr = {data.frame(censoringrateOrParameter=double(),
                                    ecp=double(),
                                    eaea=double(),
                                    eaw=double(),
                                    stringsAsFactors=FALSE)},
                          simplify = FALSE)

  names(resultsList) = c("hall-wellner",
                         "nairs-equal-precision",
                         "akritas",
                         "proposed-I",
                         "new",
                         "transformed-hall-wellner",
                         "transformed-nairs-equal-precision",
                         "transformed-akritas",
                         "proposed-III",
                         "transformed-new")
  return(resultsList)
}
