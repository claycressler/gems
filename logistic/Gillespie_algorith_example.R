## Gillespie algorith for implementing a very simple logistic model
## The model is the following
## dN/dt = (b - bs*N)*N - (d + ds*N)*N
## Thus there is density-dependence in both births and deaths

## Here is a stochastic Gillespie algorithm version of this model
## tmax = number of timesteps to run the model
## pars = vector of parameters
## NO = initial population size
## To start, let's assume that all individuals have identical parameters (b, bs, d, ds)
gillespie_logistic <- function(tmax, pars, N0) {
  b <- pars["b"]
  bs <- pars["bs"]
  d <- pars["d"]
  ds <- pars["ds"]
  
  ## initialize the time at t = 0
  t <- 0
  
  ## initialize a vector to store the time and population size
  ## Note: I am going to show you the easiest way to store the results, but it is not the fastest, computationally
  sim <- c(t, N0)

  ## current population size
  N <- N0
  
  ## start the algorithm
  while (t < tmax) {
    ## compute the birth and death rates
    brate <- (b - bs*N) * N
    drate <- (d + ds*N) * N
    rates <- c(brate, drate)
    
    ## "wheel of fortune"
    wheel <- cumsum(rates)/sum(rates)
    
    ## which event happens? Draw a random uniform to determine
    rand <- runif(1)
    
    ## if event==1, a birth happens
    ## if event==2, a death happens
    event <- 1 + sum(rand > wheel)
    if (event==1) N <- N + 1
    else N <- N - 1
    
    ## what time did the event happen?
    dt <- rexp(1, rate=sum(rates))
    
    ## update t 
    t <- t + dt
    
    ## update the results
    sim <- rbind(sim, c(t, N))
    
  }
  return(sim)
}

## Deterministic model
logistic <- function(t, y, pars) {
  b <- pars["b"]
  bs <- pars["bs"]
  d <- pars["d"]
  ds <- pars["ds"]
  
  N <- y
  
  dNdt <- (b - bs*N) * N - (d + ds*N) * N
  
  return(list(dNdt))
}

## Parameters
pars <- c(b=2, d=0.4, bs=0.01, ds=0.01)
tmax <- 100
## Initial population size
N0 <- 10

## simulate the deterministic model
library(deSolve)
det_out <- ode(y=N0, times=seq(0,tmax,0.1), func=logistic, parms=pars, method="lsoda")

## simulate the Gillespie algorithm
stoch_out <- gillespie_logistic(tmax, pars, N0)

## plot the results
plot(stoch_out[,1], stoch_out[,2], type='l', lwd=2, col=2)
lines(det_out[,1], det_out[,2], lwd=2)

## in general, you want to simulate MANY stochastic replicates and average across them
