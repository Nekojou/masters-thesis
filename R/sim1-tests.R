source("simulation-commons.R")
source("simulationstudy1-commons.R")
library(survival)
library(stats4)
library(km.ci)

set.seed(17102016)

pb<-txtProgressBar(max=nCases*nRepititions,title="Progress")

alpha2=c(seq(1.1,5.5,length.out=nCases))

# generate all samples
samples<-matrix(data.frame(), nrow=nCases, ncol=nRepititions)

for(alpha2Iterator in 1:nCases){
 for(repIterator in 1:nRepititions){
   samples[[alpha2Iterator,repIterator]]=generateRandomSample(alpha2[alpha2Iterator])
 }
}
theta.hats<-c()
boot.theta.hats<-c()
for(alpha2Iterator in 1:nCases){
    true.theta<-calculateTheta(alpha2[alpha2Iterator])

    for(repIterator in 1:nRepititions){
      aktSample<-samples[[alpha2Iterator,repIterator]]
      
      mle<-nlm(m.LL, theta<-c(1,1), x=aktSample)
      theta.hats<-c(theta.hats,mle$estimate)
      theta.hat=theta.hats[(2*repIterator-1):(2*repIterator)]
      
      for(bootIterator in 1:nBoot){
        bootSample<-twoStageBootstrap(aktSample$Z,theta.hat)
        mle.boot<-nlm(m.LL, theta<-theta.hat, x=bootSample)
        boot.theta.hats=c(boot.theta.hats,mle.boot$estimate)
      }
      
      
      setTxtProgressBar(pb, alpha2Iterator*nRepititions+repIterator)
  }
}

plot(theta.hats[seq(1,length(theta.hats),by=2)], theta.hats[seq(2,length(theta.hats),by=2)])

close(pb)