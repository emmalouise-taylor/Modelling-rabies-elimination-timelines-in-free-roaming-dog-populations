library(R6)

Parameters <- R6Class(
  "Parameters",
  portable = F,
  public = list(
    IncubationPeriod = 22.3,    ## duration that dogs are in an exposed state (days) (Beyene 2019 and Fitzpatrick 2014)  
    durationInf = 14,           ## duration that dogs are in an infected state (days)
    durationVac = 1095,         ## duration that dog vaccnation lasts (days)
    lifespan = 5*365,           ## lifespan of dogs (days) Costa 2020, Hampson 2007, Leung and Davis 2017.  
    mu =  0.006638              ## mortality rate of dogs/week (Zinsstag 2009)    
  )
)

Simulation <- R6Class(
  "Simulation",
  portable = F,
  public = list(
    p = Parameters$new(),       ## load parameters    
    nDogs = 1000,               ## number of total dogs     
    maxMov = 0.151,             ## max distnance dog travels in 1 day       
    inf.Radius = 0.1,           ## how far can a dog infect (in km) - calucaltes a pair wise distance    
    sizeM = 100,                ## map size (in km) area estimated in the coverage assessment of the 2012 vaccination campaign    
    beta = 0.07,                ## probability that biting dog is rabid    
    lengthBurn = 365*5,         ## length burn in before convergence      
    vacCov = 0.3,               ## vaccine coverage    
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
    
    initialize = function(){      ## initialize population
      status <<- rep("S",nDogs)
      inf.time <<- rep(-1,nDogs)  ## infection time reset  
      exp.time <<- rep(-1,nDogs)  ## exposed time reset
      vac.time <<- rep(-1,nDogs)  ## vaccination time reset
      dog.age <<- sample.int(p$lifespan,nDogs,replace = T) ## all dogs in population are assigned age, can have multiple dogs with same age
      mortality <<- sample.int(p$mu,nDogs,replace = T)     ## all dogs assigned a mortality rate
      
      initialExp()
      
      Xcoord <<- runif(nDogs,0,sizeM) ## movement of all dogs is random from 0 - size of map (km)
      Ycoord <<- runif(nDogs,0,sizeM) ## movement of all dogs is random from 0 - size of map (km)
      distMatrix <<- sapply(1:length(Xcoord),function(i){sqrt((Xcoord-Xcoord[i])^2+(Ycoord-Ycoord[i])^2)})
      diag(distMatrix) <<- 10000 ## hack the diagonal
      
      outputs <<- table(factor(status, levels = c("S", "E", "I", "V")))
    },
    
    initialExp = function(){  ## introduces a number of exposed dogs to initial population
      expIDs <- sample.int(nDogs,2) 
      status[expIDs] <<- "E"
      exp.time[expIDs] <<- sample.int(p$IncubationPeriod,length(expIDs),replace = T)  ## sample.int introduces a random number of days between 1 and incubation period 

      infIDs <- sample.int(nDogs,11) ## introduces a number of infected dogs to initial population
      status[infIDs] <<- "I"
      exp.time[infIDs] <<- -1
      inf.time[infIDs] <<- sample.int(p$durationInf,length(infIDs),replace = T)
    },
    
    oneTimeStep = function(){   ## executes the list of functions at each timestep 
      moveDogs()
      vacDynamics()
      infDynamics()
      expDynamics()
      susDynamics()
      ageDogs()
      outputs <<- rbind(outputs, table(factor(status, levels = c("S", "E", "I", "V"))))
    },
    
    moveDogs = function(){
      angle <- runif(nDogs,0,2*pi)  ## defines angle all dogs move randomly (angle of radians) 
      movDist <- rpert(nDogs,min=0.022,mode=0.0315,max=maxMov) 
      Xcoord <<- pmin(sizeM,pmax(0,Xcoord + movDist * cos(angle)))
      Ycoord <<- pmin(sizeM,pmax(0,Ycoord + movDist * sin(angle)))
      distMatrix <<- sapply(1:length(Xcoord),function(i){sqrt((Xcoord-Xcoord[i])^2+(Ycoord-Ycoord[i])^2)}) ## #distMatrix - assumption that dogs are close enough to interact an bite
      diag(distMatrix) <<- 10000 
    },
    
    vacDynamics = function(){ ## Assigning characteristics to "V" dogs and returning to "S" state
      vacIDs <- which(status == "V")
      vac.time[vacIDs] <<- vac.time[vacIDs]-1
      newSusIDs <- which(vac.time == 0)
      status[newSusIDs] <<- "S"
      vac.time[newSusIDs] <<- -1
    },

    infDynamics = function(){ ## Assigning characteristics to "I" dogs and death from disease
      infIDs <- which(status == "I")
      inf.time[infIDs] <<- inf.time[infIDs]-1
      deathID <- which(inf.time==0)  ## infected dogs will die once inf.time is 0
      reborn(deathID)  ## all dogs with deathID will move to reborn
    },
    
    expDynamics = function(){ ## Assigning characteristics to "E" dogs and progression to "I" state
      expIDs <- which(status == "E")
      exp.time[expIDs] <<- exp.time[expIDs]-1 
      newInfIDs <- which(exp.time == 0) ## once exp.time is 0 dogs transition to a new "I" state
      status[newInfIDs] <<- "I"
      inf.time[newInfIDs] <<- p$durationInf
      exp.time[newInfIDs] <<- -1
    },

    susDynamics = function(){  ## Assigning characteristics to "S" dogs and progression to "E" state
      susIDs <- which(status == "S")
      nInfDogs <- sapply(susIDs,proximityInf) ## those dogs which are in proximityinf will infect susIDs
      newExp <- rbinom(length(susIDs),1,1-exp(-beta*nInfDogs)) ## beta is infection rate, if no "I" dogs probability = 0
      status[susIDs[which(newExp==1)]] <<- "E"
      status[susIDs[which(newExp==1)]] <<- "E"
      exp.time[susIDs[which(newExp==1)]] <<- rpois(length(which(newExp==1)),p$IncubationPeriod) 
    },

    ## Accounting for age related death
    ## dog.age is a vector
    ageDogs = function (){
      dog.age <<- dog.age + 1
      deathID <- which(dog.age == p$lifespan)
      reborn(deathID)
    },

    reborn = function(IDs){   ## Accounting for replacement of each dog that dies to, they are reborn in to S pop - reset of individuals
      dog.age[IDs] <<- 1      ## Age reset to 1
      status[IDs] <<- "S"     ## status reset to susceptible
      inf.time[IDs] <<- -1    ## infection time reset
      exp.time[IDs] <<- -1    ## exposed time reset
      vac.time[IDs] <<- -1    ## vaccinated time reset
    },

    proximityInf = function(ID){ ## Infected dogs will infect other dogs if within inf.Radius - #checks each dog, for ones that are close by and infected
      closeIDs <- which(distMatrix[ID,]<inf.Radius) ## for certain number of dogs[ID], will look for those which are close enough to infect
      return(length(which(status[closeIDs]=="I")))  ## will return the number of "I" close by
    },
    
    burnin = function(lb=lengthBurn){ ## Burn-in the initial iterations in a Markov chain prior to its convergence to the target distribution
      i=0
      repeat{
        oneTimeStep()
        i = i + 1
        if (i==lb) {break}
      }
    },
    
    doVaccination = function(Cov = vacCov){  ## only dogs in "S" state will be vaccinated 
      susIDs <- which(status == "S")
      newVac <- rbinom(length(susIDs),1,Cov) 
      status[susIDs[which(newVac==1)]] <<- "V"
      vac.time[susIDs[which(newVac==1)]] <<- p$durationVac
    },
    
    runScenario = function(doBurn = F, lb = lengthBurn, nYears = 5, doVac = F, Cov = vacCov){
      if (doBurn == T) {burnin(lb)} ## doburn = toggles burnin as we might not want to burn in everytime
      for(i in 1:(nYears*365)){
        oneTimeStep()            
        if (doVac == T) {
          if (i %% 365 == 0){       ## check timestep i, to see if i/365 results in an intege.  test to see if we are in a tmestep where we do vacciantion - if T doVaccination is called
            doVaccination(Cov)
          }
        }
      }
    }
    
  ))



