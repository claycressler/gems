## Gillespie algorith for implementing a very simple logistic model
## The model is the following
## dN/dt = (b - bs*N)*N - (d + ds*N)*N
## Thus there is density-dependence in both births and deaths

## Here is a stochastic Gillespie algorithm version of this model
## tmax = number of timesteps to run the model
## pars = vector of parameters
## NO = initial population size
## To start, let's assume that all individuals have identical parameters (b, bs, d, ds)
gillespie_logistic <- function(tmax, tstep, pars, N0) {
  b <- pars["b"] 
  bs <- pars["bs"]
  d <- pars["d"]
  ds <- pars["ds"]
  
  ## initialize the time at t = 0
  t <- 0
  
  ## initialize a dataframe to store the time and population size
  sim <- data.frame(time=seq(0,tmax,tstep), N=rep(0, length(seq(0,tmax,tstep)))) 
  sim[1,"N"] <- N0
  
  ## current population size
  N <- N0
  
  ## start the algorithm
  while (t < tmax) {
    ## compute the birth and death rates
    brate <- (b - bs*N) * N
    drate <- (d + ds*N) * N
    rates <- c(brate, drate)
    
    ## what time does the event happen?
    dt <- rexp(1, rate=sum(rates))
    
    ## if these are not equal, don't record anything
    if(max(which(t >= seq(0,tmax,tstep)))!=max(which(t+dt > seq(0,tmax,tstep))))
      sim[max(which(t+dt > seq(0,tmax,tstep))), "N"] <- N
    
    ## update t 
    t <- t + dt
    
    ## "wheel of fortune"
    wheel <- cumsum(rates)/sum(rates)
    
    ## which event happens? Draw a random uniform to determine
    rand <- runif(1)
    
    ## if event==1, a birth happens
    ## if event==2, a death happens
    event <- 1 + sum(rand > wheel)
    if (event==1) N <- N + 1
    else N <- N - 1

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
tstep <- 0.1
## Initial population size
N0 <- 10

## simulate the deterministic model
library(deSolve)
det_out <- ode(y=N0, times=seq(0,tmax,0.1), func=logistic, parms=pars, method="lsoda")

## simulate the Gillespie algorithm
stoch_out <- gillespie_logistic(tmax, tstep, pars, N0)

## plot the results
plot(stoch_out[,1], stoch_out[,2], type='l', lwd=2, col=2)
lines(det_out[,1], det_out[,2], lwd=2)

## in general, you want to simulate MANY stochastic replicates and average across them
store <- array(NA, dim=c(length(seq(0,tmax,tstep)),21))
store[,1] <- seq(0,tmax,tstep)
for (i in 2:21) 
  store[,i] <- gillespie_logistic(tmax, tstep, pars, N0)[,"N"]

plot(store[,1], store[,2], type='l')
for (i in 3:21) lines(store[,1], store[,i])
lines(store[,1], apply(store[,2:21], 1, mean), col=2, lwd=2)
lines(det_out[,1], det_out[,2], col=3, lwd=2)
