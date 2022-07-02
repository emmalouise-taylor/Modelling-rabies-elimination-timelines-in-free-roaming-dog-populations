library(R6)

Parameters <- R6Class(
  "Parameters",
  portable = F,
  public = list(
    IncubationPeriod = 22.3,    
    durationInf = 14,           
    durationVac = 1095,         
    lifespan = 5*365,           
    mu =  0.006638                
  )
)

Simulation <- R6Class(
  "Simulation",
  portable = F,
  public = list(
    p = Parameters$new(),       
    nDogs = 1000,               
    maxMov = 0.151,                
    inf.Radius = 0.1,           
    sizeM = 100,                
    beta = 0.07,                
    lengthBurn = 365*5,          
    vacCov = 0.3,               
    status = array(),
    exp.time = array(),
    inf.time = array(),
    vac.time = array(),
    dog.age = array(),
    mortality = array(),
    Xcoord = array(),
    Ycoord = array(),
    distMatrix = matrix(),
    outputs = matrix(),
    
    initialize = function(){    
      status <<- rep("S",nDogs)
      inf.time <<- rep(-1,nDogs)  
      exp.time <<- rep(-1,nDogs) 
      vac.time <<- rep(-1,nDogs)  
      dog.age <<- sample.int(p$lifespan,nDogs,replace = T) 
      mortality <<- sample.int(p$mu,nDogs,replace = T)    
      
      initialExp()
      
      Xcoord <<- runif(nDogs,0,sizeM) 
      Ycoord <<- runif(nDogs,0,sizeM) 
      distMatrix <<- sapply(1:length(Xcoord),function(i){sqrt((Xcoord-Xcoord[i])^2+(Ycoord-Ycoord[i])^2)})
      diag(distMatrix) <<- 10000 
      
      outputs <<- table(factor(status, levels = c("S", "E", "I", "V")))
    },
    
    initialExp = function(){  
      expIDs <- sample.int(nDogs,2) 
      status[expIDs] <<- "E"
      exp.time[expIDs] <<- sample.int(p$IncubationPeriod,length(expIDs),replace = T)  

      infIDs <- sample.int(nDogs,11) 
      status[infIDs] <<- "I"
      exp.time[infIDs] <<- -1
      inf.time[infIDs] <<- sample.int(p$durationInf,length(infIDs),replace = T)
    },
    
    oneTimeStep = function(){   
      moveDogs()
      vacDynamics()
      infDynamics()
      expDynamics()
      susDynamics()
      ageDogs()
      outputs <<- rbind(outputs, table(factor(status, levels = c("S", "E", "I", "V"))))
    },
    
    moveDogs = function(){
      angle <- runif(nDogs,0,2*pi)  
      movDist <- rpert(nDogs,min=0.022,mode=0.0315,max=maxMov) 
      Xcoord <<- pmin(sizeM,pmax(0,Xcoord + movDist * cos(angle)))
      Ycoord <<- pmin(sizeM,pmax(0,Ycoord + movDist * sin(angle)))
      distMatrix <<- sapply(1:length(Xcoord),function(i){sqrt((Xcoord-Xcoord[i])^2+(Ycoord-Ycoord[i])^2)})
      diag(distMatrix) <<- 10000 
    },
    
    vacDynamics = function(){
      vacIDs <- which(status == "V")
      vac.time[vacIDs] <<- vac.time[vacIDs]-1
      newSusIDs <- which(vac.time == 0)
      status[newSusIDs] <<- "S"
      vac.time[newSusIDs] <<- -1
    },

    infDynamics = function(){
      infIDs <- which(status == "I")
      inf.time[infIDs] <<- inf.time[infIDs]-1
      deathID <- which(inf.time==0) 
      reborn(deathID)  
    },
    
    expDynamics = function(){
      expIDs <- which(status == "E")
      exp.time[expIDs] <<- exp.time[expIDs]-1
      newInfIDs <- which(exp.time == 0) 
      status[newInfIDs] <<- "I"
      inf.time[newInfIDs] <<- p$durationInf
      exp.time[newInfIDs] <<- -1
    },

    susDynamics = function(){
      susIDs <- which(status == "S")
      nInfDogs <- sapply(susIDs,proximityInf) 
      newExp <- rbinom(length(susIDs),1,1-exp(-beta*nInfDogs)) 
      status[susIDs[which(newExp==1)]] <<- "E"
      exp.time[susIDs[which(newExp==1)]] <<- rpois(length(which(newExp==1)),p$IncubationPeriod) 
    },

    ageDogs = function (){
      dog.age <<- dog.age + 1
      deathID <- which(dog.age == p$lifespan)
      reborn(deathID)
    },

    reborn = function(IDs){ 
      dog.age[IDs] <<- 1 
      status[IDs] <<- "S" 
      inf.time[IDs] <<- -1 
      exp.time[IDs] <<- -1 
      vac.time[IDs] <<- -1 
    },

    proximityInf = function(ID){ 
      closeIDs <- which(distMatrix[ID,]<inf.Radius) 
      return(length(which(status[closeIDs]=="I")))
    },
    
    burnin = function(lb=lengthBurn){ 
      i=0
      repeat{
        oneTimeStep()
        i = i + 1
        if (i==lb) {break}
      }
    },
    
    doVaccination = function(Cov = vacCov){  
      susIDs <- which(status == "S")
      newVac <- rbinom(length(susIDs),1,Cov) 
      status[susIDs[which(newVac==1)]] <<- "V"
      vac.time[susIDs[which(newVac==1)]] <<- p$durationVac
    },
    
    runScenario = function(doBurn = F, lb = lengthBurn, nYears = 5, doVac = F, Cov = vacCov){
      if (doBurn == T) {burnin(lb)} 
      for(i in 1:(nYears*365)){
        oneTimeStep()            
        if (doVac == T) {
          if (i %% 365 == 0){       
            doVaccination(Cov)
          }
        }
      }
    }
    
  ))



