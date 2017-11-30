source("simulation-commons.R")
source("simulationstudy1-commons.R")
library(survival)
library(stats4)
library(km.ci)

set.seed(17102016)

nCases=1
nBoot=25
nRepititions=10

pb<-txtProgressBar(max=nCases*nRepititions,title="Progress")

var.alpha=c(seq(1.1,5.5,length.out=nCases))

# generate all samples
samples<-matrix(data.frame(), nrow=nCases, ncol=nRepititions)

for(var.alphaIterator in 1:nCases){
  for(repIterator in 1:nRepititions){
    samples[[var.alphaIterator,repIterator]]=generateRandomSample(var.alpha[var.alphaIterator])
  }
}

for(var.alphaIterator in 1:nCases){
  
  for(repIterator in 1:nRepititions){
    aktSample<-samples[[var.alphaIterator,repIterator]]
    
    mle<-calc.mle(aktSample)
    theta.hat<-mle$par
    
    aktS.estimator2<-S.estimator2(aktSample$Z,theta.hat)
    aktS.estimator3<-S.estimator3(aktSample$Z,theta.hat)
    
    W.hat.S=c()
    W.hat.new=c()
    for(bootIterator in 1:nBoot){
      bootSample=twoStageBootstrap(aktSample$Z,theta.hat)
      boot.mle=calc.mle(bootSample)
      boot.theta.hat=boot.mle$par
      #bootS.estimator2=S.estimator2(bootSample$Z,boot.theta.hat)
      
      help=S.estimator2.t(aktS.estimator2$time,bootSample$Z,boot.theta.hat)
      help.new=S.estimator3.t(aktS.estimator3$time,bootSample$Z,boot.theta.hat)
      W.hat.S=c(W.hat.S,sqrt(samplesize)*max(abs(help-aktS.estimator2$surv)))
      W.hat.new=c(W.hat.new,sqrt(samplesize)*max(abs(help.new-aktS.estimator3$surv)))
      
    }
    
    W.hat.S=sort(W.hat.S)
    W.hat.new=sort(W.hat.new)
    q.alpha=W.hat.S[floor(nBoot*(1-alpha))]
    q.alpha.new=W.hat.new[floor(nBoot*(1-alpha))]
    
    upper.p1=aktS.estimator2$surv+1/sqrt(samplesize)*q.alpha
    lower.p1=aktS.estimator2$surv-1/sqrt(samplesize)*q.alpha
    scb.p1=list(time=aktS.estimator2$time, surv=aktS.estimator2$surv,
                upper=upper.p1, lower=lower.p1)
    
    upper.new=aktS.estimator3$surv+1/sqrt(samplesize)*q.alpha.new
    lower.new=aktS.estimator3$surv-1/sqrt(samplesize)*q.alpha.new
    scb.new=list(time=aktS.estimator3$time, surv=aktS.estimator3$surv,
                 upper=upper.new, lower=lower.new)
    
    tmp.exp=aktS.estimator2$surv*log(aktS.estimator2$surv)
    tmp.bounds=1/sqrt(samplesize)*q.alpha/tmp.exp
    upper.p3=aktS.estimator2$surv^(exp(tmp.bounds))
    lower.p3=aktS.estimator2$surv^(exp(-tmp.bounds))
    scb.p3=list(time=aktS.estimator2$time, surv=aktS.estimator2$surv,
                upper=upper.p3, lower=lower.p3)
    
    tmp.exp.new=aktS.estimator3$surv*log(aktS.estimator3$surv)
    tmp.bounds.new=1/sqrt(samplesize)*q.alpha.new/tmp.exp.new
    upper.new.t=aktS.estimator3$surv^(exp(tmp.bounds.new))
    lower.new.t=aktS.estimator3$surv^(exp(-tmp.bounds.new))
    scb.new.t=list(time=aktS.estimator3$time, surv=aktS.estimator3$surv,
                   upper=upper.new.t, lower=lower.new.t)
    
    
    setTxtProgressBar(pb, var.alphaIterator*nRepititions+repIterator)
  }
}

close(pb)

plot(scb.p1$time, scb.p1$surv)
lines(scb.p1$time, scb.p1$upper)
lines(scb.p1$time, scb.p1$lower)

plot(scb.p3$time, scb.p3$surv)
lines(scb.p3$time, scb.p3$upper)
lines(scb.p3$time, scb.p3$lower)

plot(scb.new$time, scb.new$surv)
lines(scb.new$time, scb.new$upper)
lines(scb.new$time, scb.new$lower)

plot(scb.new.t$time, scb.new.t$surv)
lines(scb.new.t$time, scb.new.t$upper)
lines(scb.new.t$time, scb.new.t$lower)