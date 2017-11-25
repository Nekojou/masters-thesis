nBoot=500
samplesize=100
nRepititions=1000
alpha=0.05
deltaT=0.001

# Obtain one Bootstrap sample of size n from Z, m and theta.hat
# m has to be a declared funktion of type m(x,theta)
twoStageBootstrap<-function(Z, theta.hat){
  Zs = sample(Z, replace = TRUE)
  deltas = sapply(lapply(Zs, m, theta=theta.hat), rbinom, n=1, size=1)
  return(data.frame(Z=Zs, delta=deltas))
}

# semiparametric estimator for S from Dikta (1999)
S.estimator2<-function(z,theta.hat){
  zs=c()
  z=sort(z)
  temp=1-m(z,theta.hat)/(samplesize-rank(z)+1)
  for(i in 1:samplesize){
    newvalue=prod(temp[z<=z[i]])
    zs=c(zs,newvalue)
  }
  return(list(time=z,surv=zs))
}

S.estimator2.t<-function(t,z,theta.hat){
  ts=c()
  z=sort(z)
  temp=1-m(z,theta.hat)/(samplesize-rank(z)+1)
  for(i in 1:length(t)){
    newvalue=prod(temp[z<=t[i]])
    ts=c(ts,newvalue)
  }
  return(ts)
}

# semiparametric estimator for S from Dikta (2016)
S.estimator3<-function(z,theta.hat){
  zs=c()
  z=sort(z)
  temp=1-m(z,theta.hat)/(samplesize-rank(z)+m(z,theta.hat))
  for(i in 1:samplesize){
    newvalue=prod(temp[z<=z[i]])
    zs=c(zs,newvalue)
  }
  return(list(time=z,surv=zs))
}

S.estimator3.t<-function(t,z,theta.hat){
  ts=c()
  z=sort(z)
  temp=1-m(z,theta.hat)/(samplesize-rank(z)+m(z,theta.hat))
  for(i in 1:samplesize){
    newvalue=prod(temp[z<=t[i]])
    ts=c(ts,newvalue)
  }
  return(ts)
}

# calculates t1 and t2 
calculateT<-function(samples){
  t=c(.Machine$double.xmin,.Machine$double.xmax)
  for(i in 1:(dim(samples)[1])){
    for(j in 1:(dim(samples)[2])){
      t[1]=max(t[1],min(samples[[i,j]]$Z))
      t[2]=min(t[2],max(samples[[i,j]]$Z))
    }
  }
  t[1]=t[1]+deltaT
  t[2]=t[2]-deltaT
  return(t)
}

calculateMs<-function(times,t.interval){
  m1=max(1,sum(times<t.interval[1]))
  m2=min(length(times),length(times)-sum(times>t.interval[2]))
  return(c(m1,m2))
}

calculateCensoringRate<-function(delta){
  return(1-mean(delta))
}

# calculate ecp, eaea and eaw for one scb
calculateCoverageAndEnclosedAreaAndWidth<-function(scb,t.interval){
  ms=calculateMs(scb$time,t.interval)
  
  Sd=S(scb$time[ms[1]:ms[2]])
  covering=all(scb$lower[ms[1]:ms[2]]<=Sd & scb$upper[ms[1]:ms[2]]>=Sd)
  
  bandwidths=(scb$upper[ms[1]:ms[2]]-scb$lower[ms[1]:ms[2]])[1:(length(scb$time[ms[1]:ms[2]])-1)]
  enclosedArea=sum(
          diff(scb$time[ms[1]:ms[2]],1)
          * bandwidths)
  
  width=sum(
    diff(-scb$surv[ms[1]:ms[2]],1)
    * bandwidths)

  return(c(covering,enclosedArea,width))
}

cropToT<-function(estimator,t)
{
  return(estimator$surv[estimator$time>=t[1] & estimator$time<=t[2]])
}

save.ergs<-function(studyID){
  write.table(ergs.hw, paste("study",studyID,"-ergs.hw.txt",sep=""))
  write.table(ergs.nep, paste("study",studyID,"-ergs.nep.txt",sep=""))
  #write.table(ergs.ak, paste("study",studyID,"-ergs.ak.txt",sep=""))
  write.table(ergs.p1, paste("study",studyID,"-ergs.p1.txt",sep=""))
  write.table(ergs.new, paste("study",studyID,"-ergs.new.txt",sep=""))
  
  write.table(ergs.hw.t, paste("study",studyID,"-ergs.hw.t.txt",sep=""))
  write.table(ergs.nep.t, paste("study",studyID,"-ergs.nep.t.txt",sep=""))
  #write.table(ergs.ak.t, paste("study",studyID,"-ergs.ak.t.txt",sep=""))
  write.table(ergs.p3, paste("study",studyID,"-ergs.p3.txt",sep=""))
  write.table(ergs.new.t, paste("study",studyID,"-ergs.new.t.txt",sep=""))
}

load.ergs<-function(studyID){
  ergs.hw<<-read.table(paste("study",studyID,"-ergs.hw.txt",sep=""))
  ergs.nep<<-read.table(paste("study",studyID,"-ergs.nep.txt",sep=""))
  #ergs.ak<<-read.table(paste("study",studyID,"-ergs.ak.txt",sep=""))
  ergs.p1<<-read.table(paste("study",studyID,"-ergs.p1.txt",sep=""))
  ergs.new<<-read.table(paste("study",studyID,"-ergs.new.txt",sep=""))
  
  ergs.hw.t<<-read.table(paste("study",studyID,"-ergs.hw.t.txt",sep=""))
  ergs.nep.t<<-read.table(paste("study",studyID,"-ergs.nep.t.txt",sep=""))
  #ergs.ak.t<<-read.table(paste("study",studyID,"-ergs.ak.t.txt",sep=""))
  ergs.p3<<-read.table(paste("study",studyID,"-ergs.p3.txt",sep=""))
  ergs.new.t<<-read.table(paste("study",studyID,"-ergs.new.t.txt",sep=""))
}

# All ergs need to be calculated before using the method run()
plot.ergs<-function(studyId)
{
  plotlims.ecp = rbind(c(0.8,1.0), c(0.85,1.0), c(0.6,1.0))
  plotlims.eaw = rbind(c(0.15,0.35), c(0.1,0.16), c(0.12,0.24))
  plotlims.eaea = rbind(c(0.6,1.1), c(0.11,0.16), c(0.35,0.5))
  
  # plot untransformed
  plot(x=ergs.hw$censoringrate,y=ergs.hw$ecp,type="l",xlab="Censoring Rate",ylab="ECP",col="blue",ylim=plotlims.ecp[studyId,])
  lines(x=ergs.nep$censoringrate,y=ergs.nep$ecp,col="green")
  # lines(x=ergs.ak$censoringrate,y=ergs.ak$ecp,col="yellow")
  lines(x=ergs.p1$censoringrate,y=ergs.p1$ecp,col="violet")
  lines(x=ergs.new$censoringrate,y=ergs.new$ecp,col="darkgreen")
  abline(0.95,0,col="red")
  title(main="Untransformed Bands")
  #legend("topleft",c("KM","95% Line"),col=c("black","red"))
  
  plot(x=ergs.hw$censoringrate,y=ergs.hw$eaea,type="l",xlab="Censoring Rate",ylab="EAEA",col="blue",ylim=plotlims.eaea[studyId,])
  lines(x=ergs.nep$censoringrate,y=ergs.nep$eaea,col="green")
  # lines(x=ergs.ak$censoringrate,y=ergs.ak$eaea,col="yellow")
  lines(x=ergs.p1$censoringrate,y=ergs.p1$eaea,col="violet")
  lines(x=ergs.new$censoringrate,y=ergs.new$eaea,col="darkgreen")
  title(main="Untransformed Bands")
  
  plot(x=ergs.hw$censoringrate,y=ergs.hw$eaw,type="l",xlab="Censoring Rate",ylab="EAW",col="blue",ylim=plotlims.eaw[studyId,])
  lines(x=ergs.nep$censoringrate,y=ergs.nep$eaw,col="green")
  # lines(x=ergs.ak$censoringrate,y=ergs.ak$eaw,col="yellow")
  lines(x=ergs.p1$censoringrate,y=ergs.p1$eaw,col="violet") 
  lines(x=ergs.new$censoringrate,y=ergs.new$eaw,col="darkgreen")
  title(main="Untransformed Bands")
  
  
  
  #plot transformed
  plot(x=ergs.hw.t$censoringrate,y=ergs.hw.t$ecp,type="l",xlab="Censoring Rate",ylab="ECP",col="blue",ylim=plotlims.ecp[studyId,])
  lines(x=ergs.nep.t$censoringrate,y=ergs.nep.t$ecp,col="green")
  # lines(x=ergs.ak.t$censoringrate,y=ergs.ak.t$ecp,col="yellow")
  lines(x=ergs.p3$censoringrate,y=ergs.p3$ecp,col="violet")
  lines(x=ergs.new.t$censoringrate,y=ergs.new.t$ecp,col="darkgreen")
  abline(0.95,0,col="red")
  title(main="Transformed Bands")
  #legend("topleft",c("KM","95% Line"),col=c("black","red"))
  
  plot(x=ergs.hw.t$censoringrate,y=ergs.hw.t$eaea,type="l",xlab="Censoring Rate",ylab="EAEA",col="blue",ylim=plotlims.eaea[studyId,])
  lines(x=ergs.nep.t$censoringrate,y=ergs.nep.t$eaea,col="green")
  # lines(x=ergs.ak.t$censoringrate,y=ergs.ak.t$eaea,col="yellow")
  lines(x=ergs.p3$censoringrate,y=ergs.p3$eaea,col="violet")
  lines(x=ergs.new.t$censoringrate,y=ergs.new.t$eaea,col="darkgreen")
  title(main="Transformed Bands")
  
  plot(x=ergs.hw.t$censoringrate,y=ergs.hw.t$eaw,type="l",xlab="Censoring Rate",ylab="EAW",col="blue",ylim=plotlims.eaw[studyId,])
  lines(x=ergs.nep.t$censoringrate,y=ergs.nep.t$eaw,col="green")
  # lines(x=ergs.ak.t$censoringrate,y=ergs.ak.t$eaw,col="yellow")
  lines(x=ergs.p3$censoringrate,y=ergs.p3$eaw,col="violet")
  lines(x=ergs.new.t$censoringrate,y=ergs.new.t$eaw,col="darkgreen")
  title(main="Transformed Bands")
}