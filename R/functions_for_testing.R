# generate some test samples
# according to simulation study 1 (weibull distributions)
# but without varying the parameters
test.generateTestSamples <- function(numberOfCases=2,montecarlo.repetitions=10)
{
  samples<-matrix(data.frame(), nrow=numberOfCases, ncol=montecarlo.repetitions)
  for(casesIterator in 1:numberOfCases){
    for(repetitionsIterator in 1:montecarlo.repetitions){
      {
        X = rweibull(100, 3)
        Y = rweibull(100, 1.5, 4.5)
        
        Z = mapply(min, X, Y)
        delta = ifelse(X<=Y,1,0)
        samples[[casesIterator,repetitionsIterator]] = data.frame(Z,delta)
      }
    }
  }
  return(samples)
}