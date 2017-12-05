# used packages
library(survival)
library(km.ci)

# common parameters
bootstrapRepetitions=500
samplesize=100
montecarloRepetitions=1000
confidencelevel=0.05

classicalUntransformedScbNames = c("hall-wellner",
                                   #"akritas",
                                   "nairs-equal-precision")
classicalTransformedScbNames = c("transformed-hall-wellner",
                                 #"transformed-akritas",
                                 "transformed-nairs-equal-precision")
classicalScbNames = c(classicalUntransformedScbNames,
                      classicalTransformedScbNames)

semiparametricUntransformedScbNames = c("proposed-I",
                                        "new")
semiparametricTransformedScbNames = c("proposed-III",
                                      "transformed-new")
semiparametricScbNames = c(semiparametricUntransformedScbNames,
                           semiparametricTransformedScbNames)
colors = c("blue", "red", "green", "violet", "darkgreen")

# Run the studies for each parameter in variableParameterList
# by calling the wrapper function applyRunOneCaseStudy
runAllStudyCases <- function(allSamples, variableParameterList, globalTimeInterval, trueSurvivalFunction)
{
  allResultsByCase = lapply(variableParameterList, applyRunOneCaseStudy,
                            variableParameterList = variableParameterList,
                            allSamples = allSamples, globalTimeInterval = globalTimeInterval,
                            trueSurvivalFunction = trueSurvivalFunction)
  
  # TODO: reorderResults(allResultsByCase)
  return(allResultsByCase)
}

# This function is a helper function to make the call of runOneCaseStudy possible with lapply
# the main problem is that both variableParameterList as well as allSamples need to be variated
applyRunOneCaseStudy <- function(parameter, variableParameterList, allSamples, globalTimeInterval, trueSurvivalFunction)
{
  return(runOneCaseStudy(parameter, allSamples[[which(variableParameterList == parameter)]], 
                         globalTimeInterval, trueSurvivalFunction))
}

# Runs the study for one specific case parameter
# Calculates the restults for all types of scbs
# groups results together in list with 
# caseParameter and censoringrate
runOneCaseStudy <- function(caseParameter, caseSamples, globalTimeInterval, trueSurvivalFunction)
{
  censoringrate = mean(sapply(caseSamples, calculateCensoringRate))
  
  resultsFromClassicalScbs = calculateResultsForClassicalScbs(caseSamples, globalTimeInterval, trueSurvivalFunction)
  
  resultsFromSubramanianScbs = list("subramanian")
  
  resultsFromNewScbs = list("new")
  
  return(c(caseParameter = caseParameter, censoringrate = censoringrate, resultsFromClassicalScbs, resultsFromSubramanianScbs, resultsFromNewScbs))
}

# Calculates empirical coverage probability,
# estimated average enclosed area and 
# estimated average width 
# for the 6 classical scb types based on 
# the Kaplan-Meier estimator
calculateResultsForClassicalScbs <- function(monteCarloSamples, globalTimeInterval, trueSurvivalFunction)
{
  coveragesEnclosedAreasWidths = lapply(monteCarloSamples, calculateClassicalResultsForOneSample, 
                                        globalTimeInterval = globalTimeInterval, 
                                        trueSurvivalFunction = trueSurvivalFunction)
  
  resultsForClassicalSCBs = lapply(classicalScbNames, applyRowMeans, coveragesEnclosedAreasWidths = coveragesEnclosedAreasWidths)
  names(resultsForClassicalSCBs) = classicalScbNames
  
  return(resultsForClassicalSCBs)
}

# This function is a helper function to make is possible to use lapply
# for calculating the rowMeans for each scb-type
applyRowMeans <- function(name, coveragesEnclosedAreasWidths)
{
  return(rowMeans(sapply(coveragesEnclosedAreasWidths, "[[", name)))
}

# Calculates the classical scbs for one specific sample.
# Returns a list of results, containing
# coverage (0 or 1), enclosed area
# and width of the scb for each type
calculateClassicalResultsForOneSample <- function(sample, globalTimeInterval, trueSurvivalFunction)
{
  # calculate Kaplan-Meier estimator
  sample.survfit = survfit(Surv(sample$Z, sample$delta)~1)
  
  scbList = list()
  #TODO; AKRITAS?
  # untransformed bands based on KM
  scbList[["hall-wellner"]]          = km.ci(sample.survfit, conf.level=1-confidencelevel, method="hall-wellner")
  scbList[["nairs-equal-precision"]] = km.ci(sample.survfit, conf.level=1-confidencelevel, method="epband")
  #scb.akritas=km.ci(sample.survfit, conf.level=1-confidencelevel, method="akritas")
  
  # transformed bands based on KM
  scbList[["transformed-hall-wellner"]]          = km.ci(sample.survfit, conf.level=1-confidencelevel, method="loghall")
  scbList[["transformed-nairs-equal-precision"]] = km.ci(sample.survfit, conf.level=1-confidencelevel, method="logep")
  #scb.akritas.transformed=km.ci(sample.survfit, conf.level=1-confidencelevel, method="logakritas")
  
  indexLimitsForStatistics = calculateIndexLimitsForStatistics(sample, globalTimeInterval)
  
  return(lapply(scbList, calculateStatistics, kmEstimator = sample.survfit,
                indexLimitsForStatistics = indexLimitsForStatistics, 
                trueSurvivalFunction = trueSurvivalFunction))
}

# Calculates the coverage, enclosed area and width 
# for a given scb within the given index limits
# TODO: open questions:
# 1. for coverage: is it correct to only evaluate at discrete times? 
# and only for m1:m2 and not t1:t2 (because upper and lower might not be defined there)?
# 2. what to do if m2 = 100 because z_j+1-z_j is not defined then
# 3. why can indexlimits exceed scbs limits?
# this was observed when hall-wellner dropped the 2 last samples. but the last one wasn't censored...
calculateStatistics <- function(scb, kmEstimator, indexLimitsForStatistics, trueSurvivalFunction)
{
  indexOffset = which(kmEstimator$time == scb$time[1])
  
  # I thought this was not necessary, but it happened that the indizes where out of range of the scb
  if(indexLimitsForStatistics[1] < indexOffset)
  {
    indexLimitsForStatistics[1] = indexOffset
  }
  if(indexLimitsForStatistics[2] > (length(scb$time) + indexOffset - 1))
  {
    indexLimitsForStatistics[2] = (length(scb$time) + indexOffset - 1)
  }
  
  indexRange = indexLimitsForStatistics[1]:indexLimitsForStatistics[2]
  indexRangeForScb = indexRange + (1 - indexOffset)
  
  trueSurvivalFunctionData = trueSurvivalFunction(kmEstimator$time[indexRange])
  
  coverage = 0
  if(all(scb$lower[indexRangeForScb]<=trueSurvivalFunctionData 
         && scb$upper[indexRangeForScb]>=trueSurvivalFunctionData))
  {
    coverage = 1
  }
  
  indexRange = indexLimitsForStatistics[1]:(indexLimitsForStatistics[2]+1)
  if(max(indexRange) > samplesize)
  {
    indexRange = indexLimitsForStatistics[1]:indexLimitsForStatistics[2]
    indexRangeForScb = indexRangeForScb[1:(length(indexRangeForScb)-1)]
  }
  
  bandwidths = scb$upper[indexRangeForScb] - scb$lower[indexRangeForScb]
  
  enclosedArea = sum(
    diff(kmEstimator$time[indexRange], 1)
    * bandwidths)
  
#  if(is.na(enclosedArea))
#  {
#    print("d")
#  }
  
  indexRange = indexLimitsForStatistics[1]:indexLimitsForStatistics[2]
  indexRangeForScb = indexRange + (1 - indexOffset)
  survforDiff = kmEstimator$surv[indexRange]
  if(min(indexRange) == 1)
  {
    survforDiff = c(1, survforDiff)
  }
  else
  {
    survforDiff = c(kmEstimator$surv[indexLimitsForStatistics[1]-1], survforDiff)
  }
  
  bandwidths = scb$upper[indexRangeForScb] - scb$lower[indexRangeForScb]
  
  width = sum(
    diff(-survforDiff, 1)
    * bandwidths)
  
#  if(is.na(width))
#  {
#    print("d")
#  }
  
  return(c(coverage=coverage, enclosedArea=enclosedArea, width=width))
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
    samples[[which(variableParameterList == parameter)]] = replicate(n = numberOfMontecarloRepetitions,
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

# reorder the results list to group results by statistic type and not by case
# returns a list with following structure
# list(censoringrates,coverage,enclosedArea,width)
# where coverage/enclosedArea/width are lists of results for each scb type
reorderResults <- function(resultsByCase, censoringrateAsParameter)
{
  coverage = list()
  for(scbName in classicalScbNames)
  {
    coverage[[scbName]] = sapply(resultsByCase, "[[", scbName)["coverage",]
  }
  
  enclosedArea = list()
  for(scbName in classicalScbNames)
  {
    enclosedArea[[scbName]] = sapply(resultsByCase, "[[", scbName)["enclosedArea",]
  }
  
  width = list()
  for(scbName in classicalScbNames)
  {
    width[[scbName]] = sapply(resultsByCase, "[[", scbName)["width",]
  }
  
  if(censoringrateAsParameter == TRUE)
  {
    censoringrates = sapply(resultsByCase, "[[", "censoringrate")
    return(list(censoringrates = censoringrates, coverage = coverage, enclosedArea = enclosedArea, width = width))
  }
  
  caseParameters = sapply(resultsByCase, "[[", "caseParameter")
  return(list(caseParameters = caseParameters, coverage = coverage, enclosedArea = enclosedArea, width = width))
}

saveResults <- function(resultsByStatistic, studyId)
{
  write.table(resultsByStatistic, paste(Sys.Date(), "_study", studyId, "_results.txt", sep = ""))
}

plotAllResults <- function(resultsByStatistic, plotLimits, censoringrateAsParameter)
{
  if(censoringrateAsParameter == TRUE)
  {
    x = resultsByStatistic[["censoringrates"]]
  }
  else
  {
    x = resultsByStatistic[["caseParameters"]]
  }
  
  untransformedNames = c(classicalUntransformedScbNames)
  generateOnePlot(x, resultsByStatistic[["coverage"]], untransformedNames, 
                  plotLimits[["coverage"]], "Untransformed Bands", "censoringrate", "ECP")
  generateOnePlot(x, resultsByStatistic[["enclosedArea"]], untransformedNames, 
                  plotLimits[["enclosedArea"]], "Untransformed Bands", "censoringrate", "EAEA")
  generateOnePlot(x, resultsByStatistic[["width"]], untransformedNames, 
                  plotLimits[["width"]], "Untransformed Bands", "censoringrate", "EAW")

  transformedNames = c(classicalTransformedScbNames)
  generateOnePlot(x, resultsByStatistic[["coverage"]], transformedNames, 
                  plotLimits[["coverage"]], "Transformed Bands", "censoringrate", "ECP")
  generateOnePlot(x, resultsByStatistic[["enclosedArea"]], transformedNames, 
                  plotLimits[["enclosedArea"]], "Transformed Bands", "censoringrate", "EAEA")
  generateOnePlot(x, resultsByStatistic[["width"]], transformedNames, 
                  plotLimits[["width"]], "Transformed Bands", "censoringrate", "EAW")
}

generateOnePlot <- function(x, listOfY, listOfYNames, yLimits, maintitle, xlabel, ylabel)
{
  plot(x=x,y=listOfY[[listOfYNames[[1]]]],type="l",xlab=xlabel,ylab=ylabel,col="blue",ylim=yLimits)
  for (i in 2:length(listOfYNames))
  {
    lines(x=x,y=listOfY[[listOfYNames[[i]]]],col=colors[i])
  }
  title(main=maintitle)
  #legend("topleft",listOfYNames,col=colors) # Does not display color
}