# generate some test samples
# according to simulation study 1 (weibull distributions)
# but without varying the parameter alpha2
test.generateTestSamples <- function(numberOfCases=2,montecarloRepetitions=10)
{
  samples<-matrix(data.frame(), nrow=numberOfCases, ncol=montecarloRepetitions)
  for(casesIterator in 1:numberOfCases){
    for(repetitionsIterator in 1:montecarloRepetitions){
      {
        X = rweibull(100, 2, 3)
        Y = rweibull(100, 1.5, 4.5)
        
        Z = mapply(min, X, Y)
        delta = ifelse(X<=Y,1,0)
        samples[[casesIterator,repetitionsIterator]] = data.frame(Z,delta)
      }
    }
  }
  return(samples)
}