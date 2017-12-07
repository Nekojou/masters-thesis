setwd(dirname(parent.frame(2)$ofile))
source("commons.R")


test.commons.runAllStudyCases <- function()
{
}

test.commons.applyRunOneCaseStudy <- function()
{
}

test.commons.runOneCaseStudy <- function()
{
}

test.commons.applyRowMeans <- function()
{
}

test.commons.calculateResultsForOneSample <- function()
{
}

# tests function calculateStatistics
# the used scb is constructed 
test.commons.calculateStatistics <- function()
{
  # test case 1
  Z = 1/10*seq(1, 100)
  delta = rbinom(100, 1, 0.9)
  
  #exclude first and last observation
  indexLimits = c(1,100)
  
  survivalfunction = function(t){ return(1-t/10)}
  survfit = survfit(Surv(Z, delta)~1)
  
  scb = list()
  scb[["time"]] = survfit$time[indexLimits[1]:indexLimits[2]]
  scb[["surv"]] = survfit$surv[indexLimits[1]:indexLimits[2]]
  scb[["upper"]] = scb$surv + 0.025
  scb[["lower"]] = scb$surv - 0.025
  
  results = calculateStatistics(scb, survfit, indexLimits, survivalfunction)

  stopifnot(!is.na(results[["coverage"]]))
  stopifnot(results[["coverage"]] == 1)
  
  stopifnot(!is.na(results[["enclosedArea"]]))
  stopifnot(!is.na(results[["width"]]))
  
  print(results)
  
  # test case 2
  #exclude first and last observations
  indexLimits = c(3,97)
  
  survivalfunction = function(t){ return(1-t)}
  
  scb = list()
  scb[["time"]] = survfit$time[indexLimits[1]:indexLimits[2]]
  scb[["surv"]] = survfit$surv[indexLimits[1]:indexLimits[2]]
  scb[["upper"]] = scb$surv + 0.025
  scb[["lower"]] = scb$surv - 0.025
  
  results = calculateStatistics(scb, survfit, indexLimits, survivalfunction)
  
  stopifnot(!is.na(results[["coverage"]]))
  stopifnot(results[["coverage"]] == 0)
  
  stopifnot(!is.na(results[["enclosedArea"]]))
  stopifnot(!is.na(results[["width"]]))
  
  print(results)
}

# test function generateAllRandomSamples
# this test is a little special because it only tests 
# if the functions itself are applied correctly to each element of the result list of lists
# and if the sizes of the lists are correct
# the function generateRandomSample is artificial and produces no useful data
# the parameterlist is of size 5
# numberOfRepisitions is chosen to be 100
test.commons.generateAllRandomSamples <- function()
{
  generateRandomSample = function(parameter){return(data.frame(Z=sample(parameter,100,replace=TRUE),delta=rep(1,100)))}
  variableParameterList = c(1.3,2.9,3,4.1,5.99)
  numberOfRepetitions = 100
  
  set.seed(456167561)
  allSamplesByFunction <- generateAllRandomSamples(generateRandomSample, variableParameterList, numberOfRepetitions)
  stopifnot(length(allSamplesByFunction) == length(variableParameterList))
  stopifnot(length(allSamplesByFunction[[which(variableParameterList == variableParameterList[1])]]) == numberOfRepetitions)
  
  set.seed(456167561)
  for(parameterIterator in 1:length(variableParameterList)){
    for(repetitionsIterator in 1:numberOfRepetitions){
        stopifnot(allSamplesByFunction[[parameterIterator]][[repetitionsIterator]] == 
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
# generate 4 sample dataframes (2 times 2 in lists)
# sample 1 has values Z from 1 to 25
# sample 2 has values Z from 26 to 50
# sample 3 has values Z from 51 to 75
# sample 4 has values Z from 76 to 100
# all delta values are set to 1
# t1 should be slightly greater than 76
# t2 should be slightly smaller than 25 
test.commons.calculateGlobalTimeInterval <- function()
{
  samples = list()
  for(casesIterator in 1:2){
    innerList = list()
    for(repetitionsIterator in 1:2){
      {
        Z = seq(1,25) + 25 * (2 * (casesIterator - 1) + (repetitionsIterator-1))
        delta = rep(1,25)
        innerList[[repetitionsIterator]] = data.frame(Z,delta)
      }
    }
    samples[[casesIterator]] = innerList
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

# test function reorderResults
# generates random resultsByCase
# and tests if they reordered results have the correct format
test.commons.reorderResults <- function(resultsByCase)
{
  subList = list()
  subList[["caseParameter"]] = 1
  subList[["censoringrate"]] = 0.2
  for(scbName in scbNames)
  {
    subList[[scbName]] = c(coverage=1, enclosedArea=2, width=3)
  }

  resultsByCase = list(subList)
  
  resultsByStatistic = reorderResults(resultsByCase)
  
  stopifnot(resultsByStatistic[["caseParameter"]] == c(1))
  stopifnot(resultsByStatistic[["censoringrate"]] == c(0.2))
  for(scbName in scbNames)
  {
    stopifnot(resultsByStatistic[["coverage"]][[scbName]] == c(1))
    stopifnot(resultsByStatistic[["enclosedArea"]][[scbName]] == c(2))
    stopifnot(resultsByStatistic[["width"]][[scbName]] == c(3))
  }
}

# test function saveResults and loadResults
# generates random data, saves it to a file
# then loads the same data again and compares it
# after that the temporary file is removed
test.commons.saveAndLoadResults <- function(resultsByStatistic, studyId)
{
  resultsByStatistic = list()
  resultsByStatistic[["censoringrate"]] = rnorm(10)
  resultsByStatistic[["caseParameter"]] = rnorm(10)
  
  for (type in statisticTypes)
  {
    subList = list()
    for (name in scbNames)
    {
      subList[[name]] = rnorm(10)
    }
    resultsByStatistic[[type]] = subList
  }
  
  saveResults(resultsByStatistic, 99)
  loadedResults = loadResults(paste(Sys.Date(), "study99_results.txt", sep="_"))
  
  stopifnot(all(floatCompare(resultsByStatistic[["censoringrate"]], loadedResults[["censoringrate"]])))
  stopifnot(all(floatCompare(resultsByStatistic[["caseParameter"]], loadedResults[["caseParameter"]])))
  
  for (type in statisticTypes)
  {
    for (name in scbNames)
    {
      stopifnot(all(floatCompare(resultsByStatistic[[type]][[name]], loadedResults[[type]][[name]])))
    }
  }
  
  invisible(file.remove(paste(Sys.Date(), "study99_results.txt", sep="_")))
}

# test function floatCompare
# tests two small floats to be equal 
# and two unequal numbers to be unequal
test.commons.floatCompare <- function()
{
  x1 = 0.00000001
  x2 = 0.00000001
  stopifnot(floatCompare(x1, x2) == TRUE)
  stopifnot(floatCompare(1, 2) == FALSE)
}

# function to run all tests listed in functionnames
test.commons.runAll <- function()
{
  functionnames = c("runAllStudyCases",
                    "applyRunOneCaseStudy",
                    "runOneCaseStudy",
                    "applyRowMeans",
                    "calculateResultsForOneSample",
                    "calculateStatistics",
                    "generateAllRandomSamples",
                    "generateBootstrapSample",
                    "calculateCensoringRate", 
                    "maxUncensoredZ", 
                    "minUncensoredZ", 
                    "calculateGlobalTimeInterval",
                    "calculateIndexLimitsForStatistics",
                    "reorderResults",
                    "saveAndLoadResults",
                    "floatCompare")
  
  print("Run tests for commons: ")
  for(fi in functionnames)
  {
    print(paste("Test function", fi, "..."))
    get(paste("test","commons",fi,sep = "."))()
    print(paste("Test function", fi, "successfull!"))
  }
  
  print("All tests successfull!")
}