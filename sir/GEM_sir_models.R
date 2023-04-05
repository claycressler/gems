
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
    birth <- b*(S+I+R) ## susceptible birth
    infect <- B0*v/(h+v)*S*I ## transmission
    Sdeath <- m*S ## susceptible death
    Ideath <- (m+v)*I ## infected death
    Rdeath <- m*R ## recovered death
    recovery <- g*I ## recovery
    events <- c(birth,Sdeath,infect,Ideath,Rdeath,recovery)
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
    else if (event==5) ## recovered death
      R <- R - 1
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

parms <- c(b = 2, m = 1.7, B0 = 0.05, h = 0.1, v = 0.4, g=0.1)
initial <- c(S=100, I=50, R=0)
tmax <- 100
dt <- 0.25

out <- sir_GEM(1242, tmax, dt, initial, parms)


## Create a new version of the above function that allows every infected individual to be unique
## In this case, every infected individual will have a unique virulence, and thus a unique transmission rate
## For this to work, we also need to specify a distribution from which virulence is drawn, and a variance
## of that distribution. Assume that virulence is drawn from a lognormal distribution.
sir_GEM2 <- function(seed, tmax, dt, initialstate, params) {
  ## Extract the initial state
  S <- unname(initialstate["S"])
  I <- unname(initialstate["I"])
  R <- unname(initialstate["R"])
  
  ## Extract the parameters - virulence is now the mean of a virulence distribution
  b <- unname(params["b"]) ## birth rate
  m <- unname(params["m"]) ## natural death rate
  B0 <- unname(params["B0"]) ## transmission rate
  h <- unname(params["h"]) ## half-saturation constant
  v <- unname(params["v"]) ## virulence
  vSD <- unname(params["vSD"]) ## SD in virulence distribution
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
  
  ## the infected population is now represented by a vector of virulence values, drawn from a lognormal distribution.
  ## Here is function for drawing individual trait values, given means and SDs on a "normal" scale. The function
  ## converts these to the lognormal scale and draws N0 random values.
  pick_individuals <- function(N0, traitmean, traitsd) {
    mu <- log(traitmean^2 / sqrt(traitsd^2+traitmean^2))
    sigma <- sqrt(log(traitsd^2/traitmean^2 + 1))
    ## record this initial distribution in the output
    return(rlnorm(N0, meanlog=mu, sdlog=sigma))
  }
  Itraits <- pick_individuals(I, traitmean=v, traitsd=vSD)
  
  i <- 1 ## how many things in the storage file? (this just keeps track of which row to put new data into)
  
  while (t < tmax & I > 0) { ## keep simulating as long as you haven't hit the final timestep and there are still infected hosts
    
    ## set up rates for each possible event
    ## we have many, many more events because each I individual has a unique virulence and a unique transmission rate
    Ideath <- (m+Itraits) ## previously, this was (m+v)*I
    infect <- B0*Itraits/(h+Itraits)*S ## previously this was B0*v/(h+v)*S*I; the new version is identical to the old version if every I has the same v value
    birth <- b*(S+I+R) ## susceptible birth
    Sdeath <- m*S ## susceptible death
    Rdeath <- m*R ## recovered death
    recovery <- g*I ## recovery
    events <- c(Ideath,infect,birth,Sdeath,Rdeath,recovery)
    ## advance time
    t <- t + exp(-1/sum(events))/sum(events)
    
    ## probability of every event
    eventProbs <- cumsum(events)/sum(events)
    ## which event happens? Draw a random uniform and figure out which event occurs
    event <- min(which(runif(1) < eventProbs))
    
    if (event <= I) {## because Ideath is a vector of length I, the first 1:I events are all deaths of a particular individual 
      ind <- event ## which individual died?
      I <- I - 1 ## remove one infected individual
      Itraits <- Itraits[-ind] ## remove that individual's traits from the Itrait vector
    }
    else if (event > I & event <= 2*I) {## because infect is a vector of length I, the next (I+1):2*I events are all infections
      S <- S - 1 ## remove one susceptible individual
      I <- I + 1 ## add one infected individual
      Itraits <- c(Itraits, pick_individuals(1, traitmean=v, traitsd=vSD)) ## add a new infected individual to trait vector
    }
    else if (event==(2*I+1)) ## birth of a new susceptible
      S <- S + 1
    else if (event==(2*I+2)) ## death of a susceptible
      S <- S - 1
    else if (event==(2*I+3)) ## recovered death
      R <- R - 1
    else {## recovery
      I <- I - 1 
      Itraits <- Itraits[-sample(1:I,1)] ## remove a random infected individual from the vector of traits
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

parms2 <- c(b = 2, m = 1.7, B0 = 0.05, h = 0.1, v = 0.4, g=0.1, vSD=0.04)

out2 <- sir_GEM2(1242, tmax, dt, initial, parms2)
