
library(rstudioapi)
library(foreach)
library(doParallel)
library(beepr)
library(pryr)
library(mc2d)

setwd(dirname(getActiveDocumentContext()$path))

source("Functions_Rabies_Model.R") ##load predefined function and parameters


set.seed(1)
nruns <- 1000
betas <- runif(nruns,0.07, 0.09) #.5,1.5 #300 numbers between X and X #this is our beta range ie from .5 to 30
inf.Radius <- runif(nruns,.15,.25) #.5, .8

cores = detectCores()
cl <- makeCluster(cores[1]-1, type = "FORK")

registerDoParallel(cl)

prevalenceV6 <- foreach(i=1:nruns, .combine=cbind) %dopar%{
  
  sim <- Simulation$new()
  sim$beta = betas[i] 
  sim$inf.Radius = inf.Radius[i]
  sim$nDogs =  1727
  sim$sizeM =  14
  
  sim$initialize() 
  sim$burnin(7305) # 1825 4 year burn in
  prev <- tail(sim$outputs,1)[3]/sim$nDogs
  
  if (prev < 0.11 & prev > 0.09){
    sim55 <- sim$clone(deep=T)
    sim55$runScenario(doBurn = F, nYears = 5, doVac = T, Cov = .55)
    sim70 <- sim$clone(deep = T)
    sim70$runScenario(doBurn = F, nYears = 5, doVac = T, Cov = .7)
    sim30 <- sim$clone(deep = T)
    sim30$runScenario(doBurn = F, nYears = 5, doVac = T, Cov = .3)
    
    return(matrix(c(sim30$outputs[,3]/sim$nDogs,
                    sim55$outputs[,3]/sim$nDogs,
                    sim70$outputs[,3]/sim$nDogs),ncol = 3))
  } else {
    return(NA)
  }
}

stopCluster(cl)

save.image("ResultsMed_V6.RData")
load("ResultsMed_V6.RData")
