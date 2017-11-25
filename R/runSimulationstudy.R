library(survival)
library(stats4)
library(km.ci)

#set.seed(17102016)

run<-function(studyId, shortRun=TRUE, enableProgress=TRUE){
  source("simulation-commons.R")
  source(paste("simulationstudy",studyId,"-commons.R",sep=""))
  
  if(shortRun)
  {
    nBoot=250
    nRepititions=100
  }
  
  if(enableProgress)
  {
    pb<-txtProgressBar(max=nCases*nRepititions,title="Progress")
  }
  
  # generate all samples
  samples<-matrix(data.frame(), nrow=nCases, ncol=nRepititions)
  
  for(var.alphaIterator in 1:nCases){
    for(repIterator in 1:nRepititions){
      samples[[var.alphaIterator,repIterator]]=generateRandomSample(var.alpha[var.alphaIterator])
    }
  }
  
  # calculate t1 and t2
  t.interval<<-calculateT(samples)
  
  # declare variables for results
  # for 5 types of confidence bands, each untransformed and transformed
  ergs.hw<-data.frame(censoringrate=double(),
                      ecp=double(),
                      eaea=double(),
                      eaw=double(),
                      stringsAsFactors=FALSE)
  ergs.nep<-data.frame(censoringrate=double(),
                       ecp=double(),
                       eaea=double(),
                       eaw=double(),
                       stringsAsFactors=FALSE)
  ergs.ak<-data.frame(censoringrate=double(),
                      ecp=double(),
                      eaea=double(),
                      eaw=double(),
                      stringsAsFactors=FALSE)
  ergs.p1<-data.frame(censoringrate=double(),
                      ecp=double(),
                      eaea=double(),
                      eaw=double(),
                      stringsAsFactors=FALSE)
  ergs.new<-data.frame(censoringrate=double(),
                       ecp=double(),
                       eaea=double(),
                       eaw=double(),
                       stringsAsFactors=FALSE)
  
  ergs.hw.t<-data.frame(censoringrate=double(),
                        ecp=double(),
                        eaea=double(),
                        eaw=double(),
                        stringsAsFactors=FALSE)
  ergs.nep.t<-data.frame(censoringrate=double(),
                         ecp=double(),
                         eaea=double(),
                         eaw=double(),
                         stringsAsFactors=FALSE)
  ergs.ak.t<-data.frame(censoringrate=double(),
                        ecp=double(),
                        eaea=double(),
                        eaw=double(),
                        stringsAsFactors=FALSE)
  ergs.p3<-data.frame(censoringrate=double(),
                      ecp=double(),
                      eaea=double(),
                      eaw=double(),
                      stringsAsFactors=FALSE)
  ergs.new.t<-data.frame(censoringrate=double(),
                         ecp=double(),
                         eaea=double(),
                         eaw=double(),
                         stringsAsFactors=FALSE)
  
  # for each parameter var.alpha
  for(var.alphaIterator in 1:nCases){
    
    censoringrate=0
    tmp.ergs=matrix(0.,nrow=10,ncol=3)
    
    # calculate estimators and ecp, eaea and eaw by doing nRepitition times 
    for(repIterator in 1:nRepititions){
      aktSample=samples[[var.alphaIterator,repIterator]]
      censoringrate=censoringrate+calculateCensoringRate(aktSample$delta)
      
      aktSample.survfit=survfit(Surv(aktSample$Z, aktSample$delta)~1)
      
      #untransformed bands based on KM
      scb.hw=km.ci(aktSample.survfit, conf.level=1-alpha, method="hall-wellner")
      scb.nep=km.ci(aktSample.survfit, conf.level=1-alpha, method="epband")
      #scb.ak=km.ci(aktSample.survfit, conf.level=1-alpha, method="akritas")
      
      #transformed bands based on KM
      scb.hw.t=km.ci(aktSample.survfit, conf.level=1-alpha, method="loghall")
      scb.nep.t=km.ci(aktSample.survfit, conf.level=1-alpha, method="logep")
      #scb.ak.t=km.ci(aktSample.survfit, conf.level=1-alpha, method="logakritas")
      
      #bootstrap based scbs
      #theta.actual=calculateTheta(var.alpha[var.alphaIterator])
      mle=calc.mle(aktSample)
      theta.hat=mle$par
      
      aktS.estimator2=S.estimator2(aktSample$Z,theta.hat)
      aktS.estimator3=S.estimator3(aktSample$Z,theta.hat)
      
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
      
      tmp.ergs[1,] = tmp.ergs[1,] + calculateCoverageAndEnclosedAreaAndWidth(scb.hw,t.interval)
      tmp.ergs[2,] = tmp.ergs[2,] + calculateCoverageAndEnclosedAreaAndWidth(scb.nep,t.interval)
      #tmp.ergs[3,] = tmp.ergs[3,] + calculateCoverageAndEnclosedAreaAndWidth(scb.ak,t.interval)
      tmp.ergs[4,] = tmp.ergs[4,] + calculateCoverageAndEnclosedAreaAndWidth(scb.p1,t.interval)
      tmp.ergs[5,] = tmp.ergs[5,] + calculateCoverageAndEnclosedAreaAndWidth(scb.new,t.interval)
      
      tmp.ergs[6,] = tmp.ergs[6,] + calculateCoverageAndEnclosedAreaAndWidth(scb.hw.t,t.interval)
      tmp.ergs[7,] = tmp.ergs[7,] + calculateCoverageAndEnclosedAreaAndWidth(scb.nep.t,t.interval)
      #tmp.ergs[8,] = tmp.ergs[8,] + calculateCoverageAndEnclosedAreaAndWidth(scb.ak.t,t.interval)
      tmp.ergs[9,] = tmp.ergs[9,] + calculateCoverageAndEnclosedAreaAndWidth(scb.p3,t.interval)
      tmp.ergs[10,] = tmp.ergs[10,] + calculateCoverageAndEnclosedAreaAndWidth(scb.new.t,t.interval)
      
      if(is.na(tmp.ergs[9,1])||is.na(tmp.ergs[9,2] || is.na(tmp.ergs[9,3])))
      {
        print(calculateCoverageAndEnclosedAreaAndWidth(scb.p3,t.interval))
      }
      
      if(enableProgress)
      {
        setTxtProgressBar(pb, var.alphaIterator*nRepititions+repIterator)
      }
      
    }
    
    censoringrate=censoringrate/nRepititions
    tmp.ergs=1/nRepititions*tmp.ergs
    
    if(studyId==1)
    {
      ergs.hw<-rbind(ergs.hw,list("censoringrate"=censoringrate,"ecp"=tmp.ergs[[1,1]],"eaea"=tmp.ergs[[1,2]],"eaw"=tmp.ergs[[1,3]]))
      ergs.nep<-rbind(ergs.nep,list("censoringrate"=censoringrate,"ecp"=tmp.ergs[[2,1]],"eaea"=tmp.ergs[[2,2]],"eaw"=tmp.ergs[[2,3]]))
      # ergs.ak<-rbind(ergs.ak,list("censoringrate"=censoringrate,"ecp"=tmp.ergs[[3,1]],"eaea"=tmp.ergs[[3,2]],"eaw"=tmp.ergs[[3,3]]))
      ergs.p1<-rbind(ergs.p1,list("censoringrate"=censoringrate,"ecp"=tmp.ergs[[4,1]],"eaea"=tmp.ergs[[4,2]],"eaw"=tmp.ergs[[4,3]]))
      ergs.new<-rbind(ergs.new,list("censoringrate"=censoringrate,"ecp"=tmp.ergs[[5,1]],"eaea"=tmp.ergs[[5,2]],"eaw"=tmp.ergs[[5,3]]))
      
      ergs.hw.t<-rbind(ergs.hw.t,list("censoringrate"=censoringrate,"ecp"=tmp.ergs[[6,1]],"eaea"=tmp.ergs[[6,2]],"eaw"=tmp.ergs[[6,3]]))
      ergs.nep.t<-rbind(ergs.nep.t,list("censoringrate"=censoringrate,"ecp"=tmp.ergs[[7,1]],"eaea"=tmp.ergs[[7,2]],"eaw"=tmp.ergs[[7,3]]))
      # ergs.ak.t<-rbind(ergs.ak.t,list("censoringrate"=censoringrate,"ecp"=tmp.ergs[[8,1]],"eaea"=tmp.ergs[[8,2]],"eaw"=tmp.ergs[[8,3]]))
      ergs.p3<-rbind(ergs.p3,list("censoringrate"=censoringrate,"ecp"=tmp.ergs[[9,1]],"eaea"=tmp.ergs[[9,2]],"eaw"=tmp.ergs[[9,3]]))
      ergs.new.t<-rbind(ergs.new.t,list("censoringrate"=censoringrate,"ecp"=tmp.ergs[[10,1]],"eaea"=tmp.ergs[[10,2]],"eaw"=tmp.ergs[[10,3]]))
    }
    else
    {
      ergs.hw<-rbind(ergs.hw,list("alpha"=var.alpha[var.alphaIterator],"ecp"=tmp.ergs[[1,1]],"eaea"=tmp.ergs[[1,2]],"eaw"=tmp.ergs[[1,3]]))
      ergs.nep<-rbind(ergs.nep,list("alpha"=var.alpha[var.alphaIterator],"ecp"=tmp.ergs[[2,1]],"eaea"=tmp.ergs[[2,2]],"eaw"=tmp.ergs[[2,3]]))
      # ergs.ak<-rbind(ergs.ak,list("alpha"=var.alpha[var.alphaIterator],"ecp"=tmp.ergs[[3,1]],"eaea"=tmp.ergs[[3,2]],"eaw"=tmp.ergs[[3,3]]))
      ergs.p1<-rbind(ergs.p1,list("alpha"=var.alpha[var.alphaIterator],"ecp"=tmp.ergs[[4,1]],"eaea"=tmp.ergs[[4,2]],"eaw"=tmp.ergs[[4,3]]))
      ergs.new<-rbind(ergs.new,list("alpha"=var.alpha[var.alphaIterator],"ecp"=tmp.ergs[[5,1]],"eaea"=tmp.ergs[[5,2]],"eaw"=tmp.ergs[[5,3]]))
      
      ergs.hw.t<-rbind(ergs.hw.t,list("alpha"=var.alpha[var.alphaIterator],"ecp"=tmp.ergs[[6,1]],"eaea"=tmp.ergs[[6,2]],"eaw"=tmp.ergs[[6,3]]))
      ergs.nep.t<-rbind(ergs.nep.t,list("alpha"=var.alpha[var.alphaIterator],"ecp"=tmp.ergs[[7,1]],"eaea"=tmp.ergs[[7,2]],"eaw"=tmp.ergs[[7,3]]))
      # ergs.ak.t<-rbind(ergs.ak.t,list("alpha"=var.alpha[var.alphaIterator],"ecp"=tmp.ergs[[8,1]],"eaea"=tmp.ergs[[8,2]],"eaw"=tmp.ergs[[8,3]]))
      ergs.p3<-rbind(ergs.p3,list("alpha"=var.alpha[var.alphaIterator],"ecp"=tmp.ergs[[9,1]],"eaea"=tmp.ergs[[9,2]],"eaw"=tmp.ergs[[9,3]]))
      ergs.new.t<-rbind(ergs.new.t,list("alpha"=var.alpha[var.alphaIterator],"ecp"=tmp.ergs[[10,1]],"eaea"=tmp.ergs[[10,2]],"eaw"=tmp.ergs[[10,3]]))
    }
    
    
    
  }
  
  if(enableProgress)
  {
    close(pb)
  }
  
  ergs.hw<<-ergs.hw[order(ergs.hw[,1]),]
  ergs.nep<<-ergs.nep[order(ergs.nep[,1]),]  
  # ergs.ak<<-ergs.ak[order(ergs.ak[,1]),]
  ergs.p1<<-ergs.p1[order(ergs.p1[,1]),]
  ergs.new<<-ergs.new[order(ergs.new[,1]),]
  
  ergs.hw.t<<-ergs.hw.t[order(ergs.hw.t[,1]),]
  ergs.nep.t<<-ergs.nep.t[order(ergs.nep.t[,1]),]
  # ergs.ak.t<<-ergs.ak.t[order(ergs.ak.t[,1]),]
  ergs.p3<<-ergs.p3[order(ergs.p3[,1]),]
  ergs.new.t<<-ergs.new.t[order(ergs.new.t[,1]),]
  
}

#run()