
## This code simulates a very simple SIR model
## The things you pass it are:
## 1. seed - a random number generator seed so you can reproduce the results exactly if you want
## 2. tmax - how many time steps do you want the simulation to run
## 3. dt - how often do you want to record the system state?
## 4. initialstate - a vector containing the initial numbers of susceptible, infected, and recovered individuals
## 4. params - a vector containing the parameter values for the model
sir_GEM <- function(seed, tmax, dt, initialstate, params) {
  ## Extract the initial state
  S <- unname(initialstate["S"])
  I <- unname(initialstate["I"])
  R <- unname(initialstate["R"])
  
  ## Extract the parameters
  b <- unname(params["b"]) ## birth rate
  m <- unname(params["m"]) ## natural death rate
  B0 <- unname(params["B0"]) ## transmission rate
  h <- unname(params["h"]) ## half-saturation constant
  v <- unname(params["v"]) ## virulence
  g <- unname(params["g"]) ## recovery rate
  
  ## set RNG seed
  set.seed(seed)
  
  ## initialize a data frame storing the system state at every moment in time
  times <- seq(0,tmax,by=dt)
  storage <- data.frame(time=times,
                        S=rep(NA,length(times)),
                        I=rep(NA,length(times)),
                        R=rep(NA,length(times)))
  ## store the initial system state
  storage[1,c("S","I","R")] <- c(S,I,R)
  t <- 0 ## current time
  
  i <- 1 ## how many things in the storage file? (this just keeps track of which row to put new data into)
  while (t < tmax & I > 0) { ## keep simulating as long as you haven't hit the final timestep and there are still infected hosts
    
    ## set up rates for each possible event
    birth <- b*(S+I) ## susceptible birth
    infect <- B0*v/(h+v)*S*I ## transmission
    Sdeath <- m*S ## susceptible death
    Ideath <- (m+v)*I ## infected death
    recovery <- g*I ## recovery
    events <- c(birth,Sdeath,infect,Ideath,recovery)
    ## advance time
    t <- t + exp(-1/sum(events))/sum(events)
    
    ## probability of every event
    eventProbs <- cumsum(events)/sum(events)
    ## which event happens? Draw a random uniform and figure out which event occurs
    event <- min(which(runif(1) < eventProbs))

    if (event==1) ## susceptible birth
      S <- S + 1
    else if (event==2) ## susceptible death
      S <- S - 1
    else if (event==3) {## new infection
      S <- S - 1
      I <- I + 1
    }
    else if (event==4) ## infected death
      I <- I - 1
    else {## recovery
      I <- I - 1
      R <- R + 1
    }
    
    ## record in the output if necessary
    if (t > storage[i+1,1]) {## if time has passed the next recording time
      storage[i+1,c("S","I","R")] <- c(S,I,R)          
      i <- i + 1 ## advance to the next timestep
    }
  }
  return(storage)
}



sir_GEM <- function(seed, tmax, dt, initialstate, traitmean, traitcv, h2, params, report=FALSE) {
  ## Extract the parameters
  b <- unname(params["b"]) ## birth rate
  m <- unname(params["m"]) ## natural death rate
  B0 <- unname(params["B0"]) ## transmission rate
  h <- unname(params["h"]) ## half-saturation constant
  
  ## set RNG seed
  set.seed(seed)
  ## initialize
  times <- seq(0,tmax,by=dt)
  t <- 0 ## current time
  ## storage for time series; at each timestep, record the system state and the traits of every individual in the population
  stateoutput <- vector("list", length(times))
  traitoutput <- vector("list", length(times))
  
  ## initial system and trait state
  currentstate <- initialstate
  ## draw initial trait distribution from lognormal with meanlog and sdlog given by traitmean and traitcv
  currenttrait <- pick_individuals(initialstate["I"], traitmean, traitsd=traitmean*traitcv)
  
  ## record the initial system and trait states
  stateoutput[[1]] <- currentstate
  traitoutput[[1]] <- currenttrait
  
  i <- 1 ## how many things in the output file?
  while (t < tmax) {
    if (report) {
      print(t)
      print(mean(currenttrait))
      print(currentstate)
    }
    if (length(currenttrait) > 0) { ## as long as the evolving population size is > 0, pick an individual
      ind <- sample(1:length(currenttrait),1)
      trait <- currenttrait[ind]
      S <- currentstate["S"]
      I <- currentstate["I"]
      ## set up rates for each possible event
      birth <- b*(S+I) ## susceptible birth
      infect <- B0*trait/(h+trait)*S*I ## transmission
      Sdeath <- m*S ## susceptible death
      Ideath <- (m+trait)*I ## infected death
      
      ## wheel of fortune
      events <- c(birth,Sdeath,infect,Ideath)
      wheel <- cumsum(events)/sum(events)
      ## which event happens?
      event <- min(which(runif(1) < wheel))
      
      if (event==1) ## susceptible birth
        currentstate["S"] <- currentstate["S"]+1
      else if (event==2) ## susceptible death
        currentstate["S"] <- currentstate["S"]-1
      else if (event==3) {## new infection
        ## the newly infected individual will have a virulence phenotype that is "inherited" from the infected individual who caused the infection
        ## compute the trait of the offspring
        ## offspring trait distribution mean
        trait_parent <- (1-h2)*mean(currenttrait) + h2*trait
        off_std <- sqrt(1-h2^2)*((1-h2)*sd(traitoutput[[1]])+h2*sd(currenttrait))
        newtrait <- pick_individuals(1, trait_parent, off_std)
        ## add this individual to the current trait and current state
        currenttrait <- c(currenttrait, newtrait)
        currentstate["S"] <- currentstate["S"]-1
        currentstate["I"] <- currentstate["I"]+1
      }
      else { ## infected death
        currenttrait <- currenttrait[-ind]
        currentstate["I"] <- currentstate["I"]-1
      }
      
      ## advance time
      t <- t + exp(-1/sum(events))/sum(events)
      
      ## record in the output if necessary
      if (t > times[i+1]) {## if time has passed the next recording time
        stateoutput[[i+1]] <- currentstate ## record the system state
        traitoutput[[i+1]] <- currenttrait
        i <- i + 1 ## advance to the next timestep
      }
    }
    else break
  }
  return(list(states=stateoutput, traits=traitoutput))
}

