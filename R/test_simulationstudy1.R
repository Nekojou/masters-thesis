this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
source("simulationstudy1.R")

# test function study1.generateRandomSample
# this mainly tests if Z and delta are constructed correctly via mapply/ifelse
# alpha2 is arbitrary chosen as 1.1
test.study1.generateRandomSample <- function()
{
  alpha2 = 1.1
  
  set.seed(27564688)
  sample = study1.generateRandomSample(alpha2)
  
  set.seed(27564688)
  X = rweibull(samplesize, study1.alpha1, study1.beta1)
  Y = rweibull(samplesize, alpha2, study1.beta2)

  Z = c()
  delta = c()
  for(i in seq(1:100))
  {
    Z = c(Z, min(X[i], Y[i]))
    if(Z[i] == X[i]) # uncensored 
    {
      delta = c(delta, 1)
    }
    else # censored
    {
      delta = c(delta, 0)
    }
  }
  
  stopifnot(sample$Z == Z)
  stopifnot(sample$delta == delta)
}

# test function study1.survivalfunction
# this mainly tests if the function works for vectors and arrays
# t is a random sample from 1 to 25 of size 100 
# alpha1 and beta1 were chosen as 1 and 2
test.study1.survivalfunction<-function(){
  alpha1 = 1
  beta1 = 2
  t = sample(25,100,replace=TRUE)
  
  ergs1 = study1.survivalfunction(t,alpha1,beta1)
  
  ergs2 = c()
  for(i in seq(1:100))
  {
    ergs2 = c(ergs2, study1.survivalfunction(t[i],alpha1,beta1))
  }
  
  ergs3 = sapply(t, study1.survivalfunction, alpha=alpha1, beta=beta1)
  
  stopifnot(ergs1 == ergs2)
  stopifnot(ergs2 == ergs3)
}

# test function study1.modelfunction
# this mainly tests if the function works for vectors and arrays
# z is a random sample from 1 to 25 of size 100 
# theta = (theta1,theta2) is chosen as (1,2)
test.study1.modelfunction<-function(){
  theta = c(1,2)
  Z = sample(25,100,replace=TRUE)
  
  ergs1 = study1.modelfunction(Z,theta)
  
  ergs2 = c()
  for(i in seq(1:100))
  {
    ergs2 = c(ergs2, study1.modelfunction(Z[i],theta))
  }
  
  ergs3 = sapply(Z, study1.modelfunction, theta=theta)
  
  stopifnot(ergs1 == ergs2)
  stopifnot(ergs2 == ergs3)
}

# function to run all tests listed in functionnames
test.study1.runAll <- function()
{
  functionnames = c("generateRandomSample",
                    "survivalfunction",
                    "modelfunction")
  
  print("Run tests for simulationstudy1: ")
  for(fi in functionnames)
  {
    print(paste("Test function", fi, "..."))
    get(paste("test","study1",fi,sep = "."))()
    print(paste("Test function", fi, "successfull!"))
  }
  
  print("All tests successfull!")
}