setwd(dirname(parent.frame(2)$ofile))
source("commons.R")


# test function generateRandomSamplesMatrix
# this test is a little special because it only tests 
# if the functions itself are applied correctly to each element of the result matrix
# the function generateRandomSample is artificial and produces no useful data
# the parameterlist is of size 5
# numberOfRepisitions is chosen to be 100
test.commons.generateRandomSamplesMatrix <- function()
{
  generateRandomSample = function(parameter){return(data.frame(Z=sample(parameter,100,replace=TRUE),delta=rep(1,100)))}
  variableParameterList = c(1,2,3,4,5)
  numberOfRepetitions = 100
  
  set.seed(456167561)
  allSamplesByFunction <- generateRandomSamplesMatrix(generateRandomSample, variableParameterList, numberOfRepetitions)
  stopifnot(nrow(allSamplesByFunction)*ncol(allSamplesByFunction) == length(variableParameterList)*numberOfRepetitions)
  
  set.seed(456167561)
  for(parameterIterator in 1:length(variableParameterList)){
    for(repetitionsIterator in 1:numberOfRepetitions){
        stopifnot(
          allSamplesByFunction[[parameterIterator,repetitionsIterator]] == 
            generateRandomSample(variableParameterList[parameterIterator]))
    }
  }
}

# test function generateBootstrapSample
# this test is a little special because it only tests if my apply functions work correctly
# this is done by comparing the results to those computed by a for loop
# a fixed random seed is used to ensure the same random results
# Z is generated randomply from numbers between 1 and 20
# deltas are 10% censored and 90% uncensored
# modelFunction is chosen as m(t,theta) = theta*z
# thetaMLE is chosen as 0.05 to let the values of m be between 0 and 1
test.commons.generateBootstrapSample <- function()
{
  Z = sample(20,100,replace=TRUE)
  delta = c(rep(0,10), rep(1,90))
  sample = data.frame(Z,delta)
  modelFunction = function(z,theta){return(theta*z)}
  thetaMLE = 0.05
  
  set.seed(76469874)
  bootstrapSample = generateBootstrapSample(sample,modelFunction,thetaMLE)
  
  set.seed(76469874)
  ZStar = sample(sample$Z, replace = TRUE)
  deltaStar = c()
  for(i in seq(1:100))
  {
    successparameter = modelFunction(ZStar[i], thetaMLE)
    deltaStar = c(deltaStar, rbinom(1,1,successparameter))
  }

  stopifnot(bootstrapSample$Z == ZStar)
  stopifnot(bootstrapSample$delta == deltaStar)
}

# test function calculateCensoringRate
# sample is artificial with 25 censored and 75 uncensored observations
# censoring rate should then be 0.25
test.commons.calculateCensoringRate <- function()
{
  Z = sample(20,100,replace=TRUE)
  delta = c(rep(0,25), rep(1,75))
  sample = data.frame(Z,delta)
  censoringrate = calculateCensoringRate(sample)
  stopifnot(censoringrate == 0.25)
}

# test function maxUncensoredZ
# Z is containing all values from 1 to 100, so index matches value
# delta is 1 for all Z but the last (100)
# maximumUncensoredZ should be 99 because 100 is censored 
test.commons.maxUncensoredZ <- function()
{
  Z = seq(1,100)
  delta = c(rep(1,99), 0)
  sample = data.frame(Z,delta)
  maximumUncensoredZ = maxUncensoredZ(sample)
  stopifnot(maximumUncensoredZ == 99)
}

# test function minUncensoredZ
# Z is containing all values from 1 to 100, so index matches value
# delta is 1 for all Z but the first (1)
# maximumUncensoredZ should be 2 because 1 is censored 
test.commons.minUncensoredZ <- function()
{
  Z = seq(1,100)
  delta = c(0, rep(1,99))
  sample = data.frame(Z,delta)
  minimumUncensoredZ = minUncensoredZ(sample)
  stopifnot(minimumUncensoredZ == 2)
}

# test function calculateGlobalTimeInterval
# generate 4 sample dataframes (2 times 2 in matrix)
# sample 1 has values Z from 1 to 25
# sample 2 has values Z from 26 to 50
# sample 3 has values Z from 51 to 75
# sample 4 has values Z from 76 to 100
# all delta values are set to 1
# t1 should be slightly greater than 76
# t2 should be slightly smaller than 25 
test.commons.calculateGlobalTimeInterval <- function()
{
  samples<-matrix(data.frame(), nrow=2, ncol=2)
  for(casesIterator in 1:2){
    for(repetitionsIterator in 1:2){
      {
        Z = seq(1,25) + 25 * (2 * (casesIterator - 1) + (repetitionsIterator-1))
        delta = rep(1,25)
        samples[[casesIterator,repetitionsIterator]] = data.frame(Z,delta)
      }
    }
  }
  
  deltaT = 0.0001
  t = calculateGlobalTimeInterval(samples)
  
  stopifnot(t[1] == 76 + deltaT)
  stopifnot(t[2] == 25 - deltaT)
}

# test function calculateIndexLimitsForStatistics
# Z is containing all values from 1 to 100, so index matches value (for comparison)
# delta is 1 for all Z
# globalTimeInterval was set to (t1,t2) = (1.5,88.9)
# indexLimit m1 should be 1
# indexLimit m2 should be 89
test.commons.calculateIndexLimitsForStatistics <- function()
{
  Z = seq(1,100)
  delta = rep(1,100)
  sample = data.frame(Z,delta)
  globalTimeInterval = c(1.5,88.9)
  
  indexLimits = calculateIndexLimitsForStatistics(sample, globalTimeInterval)
  stopifnot(indexLimits[1] == 1)
  stopifnot(indexLimits[2] == 89)
}


# function to run all tests listed in functionnames
test.commons.runAll <- function()
{
  functionnames = c("generateRandomSamplesMatrix",
                    "generateBootstrapSample",
                    "calculateCensoringRate", 
                    "maxUncensoredZ", 
                    "minUncensoredZ", 
                    "calculateGlobalTimeInterval",
                    "calculateIndexLimitsForStatistics")
  
  print("Run tests for commons: ")
  for(fi in functionnames)
  {
    print(paste("Test function", fi, "..."))
    get(paste("test","commons",fi,sep = "."))()
    print(paste("Test function", fi, "successfull!"))
  }
  
  print("All tests successfull!")
}