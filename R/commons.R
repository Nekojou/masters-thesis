setwd(dirname(parent.frame(2)$ofile))
# used packages
library(survival)
library(km.ci)

# common parameters
bootstrapRepetitions = 500
samplesize = 100
montecarloRepetitions = 1000
confidencelevel = 0.05

statisticTypes = c("coverage", "enclosedArea", "width")

untransformedScbNames = c("hall_wellner"
                          ,"akritas"
                          ,"nairs_equal_precision"
                          ,"proposed_I"
                          ,"proposed_II"
                          ,"new_proposed_I"
                          ,"new_proposed_II"
                          )
transformedScbNames = c("transformed_hall_wellner"
                        ,"transformed_akritas"
                        ,"transformed_nairs_equal_precision"
                        ,"proposed_III"
                        ,"proposed_IV"
                        ,"new_proposed_III"
                        ,"new_proposed_IV"
                        )
scbNames = c(untransformedScbNames, transformedScbNames)

colors = c("blue", "green", "orange", "violet", "darkgreen", "red")

# load functions to calculate the different estimators for S
source("commons_estimators.R")

# load functions to calculate the different scbs
source("commons_calculateScbs.R")

# Run the studies for each parameter in variableParameterList
# by calling the wrapper function applyRunOneCaseStudy
runAllStudyCases <- function(allSamples, variableParameterList, globalTimeInterval, trueSurvivalFunction, modelFunction, parameterLimits)
{
  print(sprintf("Run the following cases: %s", paste(variableParameterList, collapse=" ")))
  resultsByCase = lapply(variableParameterList, applyRunOneCaseStudy,
                            variableParameterList = variableParameterList,
                            allSamples = allSamples, 
                            globalTimeInterval = globalTimeInterval,
                            trueSurvivalFunction = trueSurvivalFunction,
                            modelFunction = modelFunction,
                            parameterLimits = parameterLimits)
  
  # reorder results to make them plottable
  resultsByStatistic = reorderResults(resultsByCase)
  return(resultsByStatistic)
}

# This function is a helper function to make the call of runOneCaseStudy possible with lapply
# the main problem is that both variableParameterList as well as allSamples need to be variated
applyRunOneCaseStudy <- function(parameter, variableParameterList, allSamples, globalTimeInterval, trueSurvivalFunction, modelFunction, parameterLimits)
{
  return(runOneCaseStudy(parameter, allSamples[[which(variableParameterList == parameter)]], 
                         globalTimeInterval, trueSurvivalFunction, modelFunction, parameterLimits))
}

# Runs the study for one specific case parameter
# Calculates empirical coverage probability,
# estimated average enclosed area and 
# estimated average width 
# for all scb types based on montecarlo simulations
# groups results together in list with 
# caseParameter and censoringrate
runOneCaseStudy <- function(caseParameter, caseSamples, globalTimeInterval, trueSurvivalFunction, modelFunction, parameterLimits)
{
  print(sprintf("Start case %f", caseParameter))

  censoringrate = mean(sapply(caseSamples, calculateCensoringRate))

  coveragesEnclosedAreasWidths = parLapply(cl, caseSamples, calculateResultsForOneSample, 
                                        globalTimeInterval = globalTimeInterval, 
                                        trueSurvivalFunction = trueSurvivalFunction,
                                        modelFunction = modelFunction,
                                        parameterLimits = parameterLimits)
  
  results = lapply(scbNames, applyRowMeans, coveragesEnclosedAreasWidths = coveragesEnclosedAreasWidths)
  names(results) = scbNames
  
  return(c(caseParameter = caseParameter, censoringrate = censoringrate, results))
}

# This function is a helper function to make is possible to use lapply
# for calculating the rowMeans for each scb-type
applyRowMeans <- function(name, coveragesEnclosedAreasWidths)
{
  return(rowMeans(sapply(coveragesEnclosedAreasWidths, "[[", name)))
}

# Calculates the scbs for one specific sample.
# Returns a list of results, containing
# coverage (0 or 1), enclosed area
# and width of the scb for each type
calculateResultsForOneSample <- function(sample, globalTimeInterval, trueSurvivalFunction, modelFunction, parameterLimits)
{
  # calculate Kaplan-Meier estimator
  sample.kaplan_meier_estimator = survfit(Surv(sample$Z, sample$delta)~1)
  
  # calculate semiparametric estimators
  mleTheta = estimators.calculateMaximumLikelihoodEstimator(sample, modelFunction, parameterLimits)
  
  sample.dikta_2_estimator = estimators.dikta_2(sample, modelFunction, mleTheta)
  sample.dikta_3_estimator = estimators.dikta_3(sample, modelFunction, mleTheta)
  
  scbList = list()
  # untransformed bands
  scbList[["hall_wellner"]]          = calculateScb_hall_wellner(sample.kaplan_meier_estimator)
  scbList[["akritas"]]               = calculateScb_akritas(sample.kaplan_meier_estimator)
  scbList[["nairs_equal_precision"]] = calculateScb_nairs_equal_precision(sample.kaplan_meier_estimator)
  scbList[["proposed_I"]]            = calculateScb_proposed_I(sample.dikta_2_estimator, sample, modelFunction, mleTheta, parameterLimits)
  scbList[["proposed_II"]]           = calculateScb_proposed_II(sample.dikta_2_estimator)
  scbList[["new_proposed_I"]]        = calculateScb_new_proposed_I(sample.dikta_3_estimator, sample, modelFunction, mleTheta, parameterLimits)
  scbList[["new_proposed_II"]]       = calculateScb_new_proposed_II(sample.dikta_3_estimator)

  # transformed bands
  scbList[["transformed_hall_wellner"]]           = calculateScb_transformed_hall_wellner(sample.kaplan_meier_estimator)
  scbList[["transformed_akritas"]]                = calculateScb_akritas(sample.kaplan_meier_estimator)
  scbList[["transformed_nairs_equal_precision"]]  = calculateScb_transformed_nairs_equal_precision(sample.kaplan_meier_estimator)
  scbList[["proposed_III"]]                       = calculateScb_proposed_III(sample.dikta_2_estimator, sample, modelFunction, mleTheta, parameterLimits)
  scbList[["proposed_IV"]]                        = calculateScb_proposed_IV(sample.dikta_2_estimator)
  scbList[["new_proposed_III"]]                   = calculateScb_new_proposed_III(sample.dikta_3_estimator, sample, modelFunction, mleTheta, parameterLimits)
  scbList[["new_proposed_IV"]]                    = calculateScb_new_proposed_IV(sample.dikta_3_estimator)

  indexLimitsForStatistics = calculateIndexLimitsForStatistics(sample, globalTimeInterval)

  return(lapply(scbList, calculateStatistics,
                indexLimitsForStatistics = indexLimitsForStatistics, 
                trueSurvivalFunction = trueSurvivalFunction))
}

# Calculates the coverage, enclosed area and width 
# for a given scb within the given index limits
# TODO: why can indexlimits exceed scbs limits?
# this was observed when hall_wellner dropped the 2 last samples. but the last one wasn't censored...
calculateStatistics <- function(scb, indexLimitsForStatistics, trueSurvivalFunction)
{
  indexOffset = which(scb$estimator$time == scb$time[1])
  
  # I thought this was not necessary, but it happened that the indizes where out of range of the scb
  # TODO: Might not be necessary anymore when KM based estimators are manually calculated
  if(indexLimitsForStatistics[1] < indexOffset)
  {
    indexLimitsForStatistics[1] = indexOffset
  }
  if(indexLimitsForStatistics[2] > (length(scb$time) + indexOffset - 1))
  {
    indexLimitsForStatistics[2] = (length(scb$time) + indexOffset - 1)
  }
  
  coverage = calculateCoverage(scb, scb$estimator, indexLimitsForStatistics, indexOffset, trueSurvivalFunction)
  
  enclosedArea = calculateEnclosedArea(scb, scb$estimator, indexLimitsForStatistics, indexOffset)

  width = calculateWidth(scb, scb$estimator, indexLimitsForStatistics, indexOffset)
  
  return(c(coverage = coverage, enclosedArea = enclosedArea, width = width))
}

# function to check if an SCB completely includes the true survival function
# returns 1 if so, and 0 if not
# TODO: calculate on t1,t2 not m1,m2? Might be possible when KM based estimators are manually calculated
# TODO: Do not evaluate on descrete times only but develop a strategy to check this correctly
calculateCoverage <- function(scb, estimator, indexLimitsForStatistics, indexOffset, trueSurvivalFunction)
{
  indexRange = indexLimitsForStatistics[1]:indexLimitsForStatistics[2]
  indexRangeForScb = indexRange + (1 - indexOffset)
  
  trueSurvivalFunctionData = trueSurvivalFunction(estimator$time[indexRange])
  
  coverage = 0
  if(all(scb$lower[indexRangeForScb] <= trueSurvivalFunctionData 
         && scb$upper[indexRangeForScb] >= trueSurvivalFunctionData
         && scb$lower[indexRangeForScb][1:(length(trueSurvivalFunctionData)-1)]
                <= trueSurvivalFunctionData[2:length(trueSurvivalFunctionData)]))
  {
    coverage = 1
  }
  
  return(coverage)
}

# calculate the enclosed area of a specific scb
# TODO: what to do if m2 = 100 because z_j+1-z_j is not defined then
# and this automatically makes the enclosed area rationally smaller, because one square less is included
calculateEnclosedArea <- function(scb, estimator, indexLimitsForStatistics, indexOffset)
{
  indexRange = indexLimitsForStatistics[1]:(indexLimitsForStatistics[2]+1)
  indexRangeForScb = indexLimitsForStatistics[1]:indexLimitsForStatistics[2] + (1 - indexOffset)
  if(max(indexRange) > samplesize)
  {
    indexRange = indexLimitsForStatistics[1]:indexLimitsForStatistics[2]
    indexRangeForScb = indexRangeForScb[1:(length(indexRangeForScb)-1)]
  }
  
  bandwidths = scb$upper[indexRangeForScb] - scb$lower[indexRangeForScb]
  
  enclosedArea = sum(
    diff(estimator$time[indexRange], 1)
    * bandwidths)
  
  return(enclosedArea)
}

# calculate width of a specific scb
calculateWidth <- function(scb, estimator, indexLimitsForStatistics, indexOffset)
{
  indexRange = indexLimitsForStatistics[1]:indexLimitsForStatistics[2]
  indexRangeForScb = indexRange + (1 - indexOffset)
  survforDiff = estimator$surv[indexRange]
  if(min(indexRange) == 1)
  {
    survforDiff = c(1, survforDiff)
  }
  else
  {
    survforDiff = c(estimator$surv[indexLimitsForStatistics[1]-1], survforDiff)
  }
  
  bandwidths = scb$upper[indexRangeForScb] - scb$lower[indexRangeForScb]
  
  width = sum(
    diff(-survforDiff, 1)
    * bandwidths)
  
  return(width)
}

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
# with success parameter modelFunction(ZStar_i, mleTheta)
generateBootstrapSample <- function(originalSample, modelFunction, mleTheta)
{
  ZStar = sample(originalSample$Z, replace = TRUE)
  deltaStar = sapply(sapply(ZStar, modelFunction, theta=mleTheta), rbinom, n=1, size=1)
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

# reorder the results list to group results by statistic type and not by case
# returns a list with following structure
# list(censoringrate,caseParameter,coverage,enclosedArea,width)
# where coverage/enclosedArea/width are lists of results for each scb type
reorderResults <- function(resultsByCase)
{
  coverage = list()
  enclosedArea = list()
  width = list()
  
  for(scbName in scbNames)
  {
    coverage[[scbName]] = sapply(resultsByCase, "[[", scbName)["coverage",]
    enclosedArea[[scbName]] = sapply(resultsByCase, "[[", scbName)["enclosedArea",]
    width[[scbName]] = sapply(resultsByCase, "[[", scbName)["width",]
  }
  
  censoringrate = sapply(resultsByCase, "[[", "censoringrate")
  caseParameter = sapply(resultsByCase, "[[", "caseParameter")
  return(list(censoringrate = censoringrate, caseParameter = caseParameter, coverage = coverage, enclosedArea = enclosedArea, width = width))

}

saveResults <- function(resultsByStatistic, studyId)
{
  write.table(resultsByStatistic, paste(Sys.Date(), "_study", studyId, "_results.txt", sep = ""))
}

loadResults <- function(filename)
{
  resultsInOtherFormat = read.table(filename)
  
  resultsByStatistic = list()
  resultsByStatistic[["censoringrate"]] = resultsInOtherFormat$censoringrate
  resultsByStatistic[["caseParameter"]] = resultsInOtherFormat$caseParameter
  
  for (type in statisticTypes)
  {
    subList = list()
    for (name in scbNames)
    {
      subList[[name]] = resultsInOtherFormat[[gsub("_", ".", paste(type,name,sep = "."))]]
    }
    resultsByStatistic[[type]] = subList
  }
  
  return(resultsByStatistic)
}

plotAllResults <- function(resultsByStatistic, plotLimits, censoringrateAsParameter)
{
  if(censoringrateAsParameter == TRUE)
  {
    x = resultsByStatistic[["censoringrate"]]
    xLabel = "censoringrate"
  }
  else
  {
    x = resultsByStatistic[["caseParameter"]]
    xLabel = "caseParameter"
  }
  
  generateOnePlot(x, resultsByStatistic[["coverage"]], untransformedScbNames, 
                  plotLimits[["coverage"]], "Untransformed Bands", xLabel, "ECP")
  generateOnePlot(x, resultsByStatistic[["enclosedArea"]], untransformedScbNames, 
                  plotLimits[["enclosedArea"]], "Untransformed Bands", xLabel, "EAEA")
  generateOnePlot(x, resultsByStatistic[["width"]], untransformedScbNames, 
                  plotLimits[["width"]], "Untransformed Bands", xLabel, "EAW")

  generateOnePlot(x, resultsByStatistic[["coverage"]], transformedScbNames, 
                  plotLimits[["coverage"]], "Transformed Bands", xLabel, "ECP")
  generateOnePlot(x, resultsByStatistic[["enclosedArea"]], transformedScbNames, 
                  plotLimits[["enclosedArea"]], "Transformed Bands", xLabel, "EAEA")
  generateOnePlot(x, resultsByStatistic[["width"]], transformedScbNames, 
                  plotLimits[["width"]], "Transformed Bands", xLabel, "EAW")
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

floatCompare <- function(in1, in2)
{
  return(abs(in1-in2) <= 1e-8)
}