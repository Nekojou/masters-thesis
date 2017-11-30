# used packages


# common parameters
bootstrapRepetitions=500
samplesize=100
montecarloRepetitions=1000
confidencelevel=0.05


# TODO: Functions to be transferred:
# S.estimator2<-function(z,theta.hat){}
# S.estimator2.t<-function(t,z,theta.hat){}
# S.estimator3<-function(z,theta.hat){}
# S.estimator3.t<-function(t,z,theta.hat){}
# calculateCoverageAndEnclosedAreaAndWidth<-function(scb,globalTimeInterval){}
# save.ergs<-function(prefix){}
# load.ergs<-function(prefix){}
# plot.ergs<-function(prefix){}

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
calculateGlobalTimeInterval <- function(allSamplesMatrix)
{
  deltaT = 0.0001
  
  t1 = max(sapply(allSamplesMatrix, minUncensoredZ)) + deltaT
  t2 = min(sapply(allSamplesMatrix, maxUncensoredZ)) - deltaT
  
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

