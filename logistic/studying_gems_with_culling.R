pick_individuals <- function(N0, traitmean, traitsd) {
  mu <- log(traitmean^2 / sqrt(traitsd^2+traitmean^2))
  sigma <- sqrt(log(traitsd^2/traitmean^2 + 1))
  ## record this initial distribution in the output
  return(rlnorm(N0, meanlog=mu, sdlog=sigma))
}

## initialize the population at the expected eco-evolutionary equilibrium, given the slope, bs, and ds parameters
seed <- 10110
tmax <- 100
bs <- ds <- 0.01
slope <- 0.1

## evolutionary equilibrum: bmax = 1/2s = 5
## ecological equilibrium: N = (bmax-dmin)/(bs-ds) = (1/2s - s*(1/2s)^2)/(bs+ds) = 1/(4s*(bs+ds)) = 125
N0 <- 1/(4*slope*(bs+ds))
traitmean <- 1/(2*slope)
traitcv <- 0.1
h2 <- 0.75

#######################################################################################################
#######################################################################################################
#######################################################################################################
#
#                           MODEL 1: Culling is done as a rate               
#
#######################################################################################################
#######################################################################################################
#######################################################################################################

logistic_GEM_cull <- function(seed, dt, tmax, N0, traitmean, traitcv, h2, bs, ds, slope, X) {
  
  ## set RNG seed
  set.seed(seed)
  
  ## storage for times and trait distributions
  ## the maximum possible size is tmax/dt+1
  time <- rep(0,tmax/dt+1)
  output <- vector(mode='list', length=tmax/dt+1)
  
  ## initialize counter
  i <- 1
  ## initialize time
  t <- 0
  lastrecordtime <- t ## the last timepoint that was recorded
  time[i] <- t
  ## initialize trait distribution from lognormal with meanlog and sdlog given by traitmean and traitcv
  output[[1]] <- pick_individuals(N0, traitmean, traitsd=traitmean*traitcv)
  ## also record it in a vector called 'currentstate' that records the current state of the system only
  currentstate <- output[[1]]
  
  ## if you want to track aspects of the output (for diagnostic purposes), use this to create a text file named for the seed
  ## outputfile <- paste0("/Users/claycressler/Box Sync/Research/John DeLong/seed",seed,".txt")
  
  while (t < tmax & length(currentstate) > 0) { ## as long as population size is > 0 pick an individual
    ind <- sample(1:length(currentstate),1)
    trait <- currentstate[ind]
    N <- length(currentstate)
    ## set up rates for each possible event
    ## birth
    brate <- (trait - bs*N)*N
    ## death (specified by tradeoff between birth and death)
    d <- slope*trait^2
    drate <- (d + ds*N)*N
    ## culling is just X
    cullrate <- X*N
    events <- c(brate,drate,cullrate)
    
    ## wheel of fortune
    wheel <- cumsum(events)/sum(events)
    ## which event happens?
    event <- min(which(runif(1) < wheel))
    ## when did this event happen?
    timestep <- exp(-1/sum(events))/sum(events)
    t <- t + timestep
    
    if (event==1) { ## birth
      ## compute the trait of the offspring
      trait_parent <- (1-h2)*mean(currentstate) + h2*trait
      ## offspring trait distribution standard deviation
      ## if there is more than one individual, this will work
      off_std <- sqrt(1-h2^2)*((1-h2)*sd(output[[1]])+h2*sd(currentstate))
      ## if there is only one individual left in the population, set off_std = 0
      if (is.na(off_std)) off_std <- 0
      ## generate the new individual
      newtrait <- pick_individuals(1, trait_parent, off_std)
      ## add this individual to the current state vector
      currentstate <- c(currentstate, newtrait)
    }
    else  ## death
      currentstate <- currentstate[-ind]
    
    ## record... but do it smartly, to prevent output from getting so big that it crashes R.
    ## Record every event if it moves time forward more than dt time steps, otherewise only record every dt timesteps
    if (timestep > dt) {
      i <- i + 1
      time[i] <- t
      output[[i]] <- currentstate
      lastrecordtime <- lastrecordtime + dt
    }
    else if ((floor(t/dt)*dt) > lastrecordtime) {
      i <- i + 1
      time[i] <- t
      output[[i]] <- currentstate
      lastrecordtime <- floor(t/dt)*dt
    }
    
    ## if you're doing diagnostics, write to the text file
    #writeLines(paste0("\ntrait=", round(trait,2), " bs=", bs, " ds=", ds, " t=",round(t,2)," pop size=",length(currentstate)," memUsed=",mem_used()), con=outputfile)
    
  }
  ## trim down to only what you need
  output <- output[1:i]
  time <- time[1:i]
  
  return(list(time=time,traits=output))
}


set.seed(123497)
seeds <- floor(runif(40,1,1e6))
library(parallel)
mclapply(seeds,
         function(s) logistic_GEM_cull(s, dt=0.1, tmax, N0, traitmean, traitcv, h2, bs, ds, slope,X=1),
         mc.cores=8                     
) -> out

Ntraj <- lapply(out, function(o) lapply(o[[2]],length) %>% unlist) %>% unlist %>% matrix(., ncol=length(out), byrow=F)
btraj <- lapply(out, function(o) lapply(o[[2]],mean) %>% unlist) %>% unlist %>% matrix(., ncol=length(out), byrow=F)
time.seq <- seq(0,tmax,0.1)

## plot the mean trajectories across simulations:
par(mfrow=c(2,1), mar=c(4,4,0.5,0.5), oma=rep(0.5,4))
plot.new() 
plot.window(xlim=c(0,tmax), ylim=c(30,125))
axis(1);axis(2);box('plot')
mtext(side=1, "Time", line=3)
mtext(side=2, "Population size", line=2.5)
for (i in 1:length(out)) lines(time.seq, Ntraj[,i], col=gray(0.7))
lines(time.seq, apply(Ntraj, 1, mean), lwd=2)

plot.new() 
plot.window(xlim=c(0,tmax), ylim=c(3,7))
axis(1);axis(2);box('plot')
mtext(side=1, "Time", line=3)
mtext(side=2, "bmax", line=2.5)
for (i in 1:length(out)) lines(time.seq, btraj[,i], col=gray(0.7))
lines(time.seq, apply(btraj, 1, mean), lwd=2)

#######################################################################################################
#######################################################################################################
#######################################################################################################
#
#                           MODEL 2: Culling is done as a number/time rather than a rate               
#
#######################################################################################################
#######################################################################################################
#######################################################################################################
## implement a cull with X = 30 ind/time
logistic_GEM_cull <- function(seed, dt, tmax, N0, traitmean, traitcv, h2, bs, ds, slope, X) {
  
  ## set RNG seed
  set.seed(seed)
  
  ## storage for times and trait distributions
  ## the maximum possible size is tmax/dt+1
  time <- rep(0,tmax/dt+1)
  output <- vector(mode='list', length=tmax/dt+1)
  
  ## initialize counter
  i <- 1
  ## initialize time
  t <- 0
  lastrecordtime <- t ## the last timepoint that was recorded
  time[i] <- t
  ## initialize trait distribution from lognormal with meanlog and sdlog given by traitmean and traitcv
  output[[1]] <- pick_individuals(N0, traitmean, traitsd=traitmean*traitcv)
  ## also record it in a vector called 'currentstate' that records the current state of the system only
  currentstate <- output[[1]]
  
  ## if you want to track aspects of the output (for diagnostic purposes), use this to create a text file named for the seed
  ## outputfile <- paste0("/Users/claycressler/Box Sync/Research/John DeLong/seed",seed,".txt")
  
  while (t < tmax & length(currentstate) > 0) { ## as long as population size is > 0 pick an individual
    ind <- sample(1:length(currentstate),1)
    trait <- currentstate[ind]
    N <- length(currentstate)
    ## set up rates for each possible event
    ## birth
    brate <- (trait - bs*N)*N
    ## death (specified by tradeoff between birth and death)
    d <- slope*trait^2
    drate <- (d + ds*N)*N
    ## culling is just X
    cullrate <- X
    events <- c(brate,drate,cullrate)
    
    ## wheel of fortune
    wheel <- cumsum(events)/sum(events)
    ## which event happens?
    event <- min(which(runif(1) < wheel))
    ## when did this event happen?
    timestep <- exp(-1/sum(events))/sum(events)
    t <- t + timestep
    
    if (event==1) { ## birth
      ## compute the trait of the offspring
      trait_parent <- (1-h2)*mean(currentstate) + h2*trait
      ## offspring trait distribution standard deviation
      ## if there is more than one individual, this will work
      off_std <- sqrt(1-h2^2)*((1-h2)*sd(output[[1]])+h2*sd(currentstate))
      ## if there is only one individual left in the population, set off_std = 0
      if (is.na(off_std)) off_std <- 0
      ## generate the new individual
      newtrait <- pick_individuals(1, trait_parent, off_std)
      ## add this individual to the current state vector
      currentstate <- c(currentstate, newtrait)
    }
    else  ## death
      currentstate <- currentstate[-ind]
    
    ## record... but do it smartly, to prevent output from getting so big that it crashes R.
    ## Record every event if it moves time forward more than dt time steps, otherewise only record every dt timesteps
    if (timestep > dt) {
      i <- i + 1
      time[i] <- t
      output[[i]] <- currentstate
      lastrecordtime <- lastrecordtime + dt
    }
    else if ((floor(t/dt)*dt) > lastrecordtime) {
      i <- i + 1
      time[i] <- t
      output[[i]] <- currentstate
      lastrecordtime <- floor(t/dt)*dt
    }
    
    ## if you're doing diagnostics, write to the text file
    #writeLines(paste0("\ntrait=", round(trait,2), " bs=", bs, " ds=", ds, " t=",round(t,2)," pop size=",length(currentstate)," memUsed=",mem_used()), con=outputfile)
    
  }
  ## trim down to only what you need
  output <- output[1:i]
  time <- time[1:i]
  
  return(list(time=time,traits=output))
}

set.seed(123497)
seeds <- floor(runif(40,1,1e6))
library(parallel)
mclapply(seeds,
         function(s) logistic_GEM_cull(s, dt=0.1, tmax, N0, traitmean, traitcv, h2, bs, ds, slope,X=30),
         mc.cores=8                     
) -> out2

Ntraj <- lapply(out2, function(o) lapply(o[[2]],length) %>% unlist) %>% unlist %>% matrix(., ncol=length(out2), byrow=F)
btraj <- lapply(out2, function(o) lapply(o[[2]],mean) %>% unlist) %>% unlist %>% matrix(., ncol=length(out2), byrow=F)

## plot the mean trajectories across simulations:
par(mfrow=c(2,1), mar=c(4,4,0.5,0.5), oma=rep(0.5,4))
plot.new() 
plot.window(xlim=c(0,200), ylim=c(50,150))
axis(1);axis(2);box('plot')
mtext(side=1, "Time", line=3)
mtext(side=2, "Population size", line=2.5)
for (i in 1:length(out2)) lines(time.seq, Ntraj[,i], col=gray(0.7))
lines(time.seq, apply(Ntraj, 1, mean), lwd=2)

plot.new() 
plot.window(xlim=c(0,200), ylim=c(3,7))
axis(1);axis(2);box('plot')
mtext(side=1, "Time", line=3)
mtext(side=2, "bmax", line=2.5)
for (i in 1:length(out2)) lines(time.seq, btraj[,i], col=gray(0.7))
lines(time.seq, apply(btraj, 1, mean), lwd=2)


#######################################################################################################
#######################################################################################################
#######################################################################################################
#
#                           MODEL 3: Culling is done as a rate, but only when N > Ncull               
#
#######################################################################################################
#######################################################################################################
#######################################################################################################
## implement a cull with X = 30 ind/time
logistic_GEM_cull <- function(seed, dt, tmax, N0, traitmean, traitcv, h2, bs, ds, slope, X, Ncull) {
  
  ## set RNG seed
  set.seed(seed)
  
  ## storage for times and trait distributions
  ## the maximum possible size is tmax/dt+1
  time <- rep(0,tmax/dt+1)
  output <- vector(mode='list', length=tmax/dt+1)
  
  ## initialize counter
  i <- 1
  ## initialize time
  t <- 0
  lastrecordtime <- t ## the last timepoint that was recorded
  time[i] <- t
  ## initialize trait distribution from lognormal with meanlog and sdlog given by traitmean and traitcv
  output[[1]] <- pick_individuals(N0, traitmean, traitsd=traitmean*traitcv)
  ## also record it in a vector called 'currentstate' that records the current state of the system only
  currentstate <- output[[1]]
  
  ## if you want to track aspects of the output (for diagnostic purposes), use this to create a text file named for the seed
  ## outputfile <- paste0("/Users/claycressler/Box Sync/Research/John DeLong/seed",seed,".txt")
  
  while (t < tmax & length(currentstate) > 0) { ## as long as population size is > 0 pick an individual
    ind <- sample(1:length(currentstate),1)
    trait <- currentstate[ind]
    N <- length(currentstate)
    ## set up rates for each possible event
    ## birth
    brate <- (trait - bs*N)*N
    ## death (specified by tradeoff between birth and death)
    d <- slope*trait^2
    drate <- (d + ds*N)*N
    ## the cullrate is X*(N-Ncull) if N > Ncull, and 0 otherwise
    cullrate <- max(X*(N-Ncull), 0)
    events <- c(brate,drate,cullrate)
    
    ## wheel of fortune
    wheel <- cumsum(events)/sum(events)
    ## which event happens?
    event <- min(which(runif(1) < wheel))
    ## when did this event happen?
    timestep <- exp(-1/sum(events))/sum(events)
    t <- t + timestep
    
    if (event==1) { ## birth
      ## compute the trait of the offspring
      trait_parent <- (1-h2)*mean(currentstate) + h2*trait
      ## offspring trait distribution standard deviation
      ## if there is more than one individual, this will work
      off_std <- sqrt(1-h2^2)*((1-h2)*sd(output[[1]])+h2*sd(currentstate))
      ## if there is only one individual left in the population, set off_std = 0
      if (is.na(off_std)) off_std <- 0
      ## generate the new individual
      newtrait <- pick_individuals(1, trait_parent, off_std)
      ## add this individual to the current state vector
      currentstate <- c(currentstate, newtrait)
    }
    else  ## death
      currentstate <- currentstate[-ind]
    
    ## record... but do it smartly, to prevent output from getting so big that it crashes R.
    ## Record every event if it moves time forward more than dt time steps, otherewise only record every dt timesteps
    if (timestep > dt) {
      i <- i + 1
      time[i] <- t
      output[[i]] <- currentstate
      lastrecordtime <- lastrecordtime + dt
    }
    else if ((floor(t/dt)*dt) > lastrecordtime) {
      i <- i + 1
      time[i] <- t
      output[[i]] <- currentstate
      lastrecordtime <- floor(t/dt)*dt
    }
    
    ## if you're doing diagnostics, write to the text file
    #writeLines(paste0("\ntrait=", round(trait,2), " bs=", bs, " ds=", ds, " t=",round(t,2)," pop size=",length(currentstate)," memUsed=",mem_used()), con=outputfile)
    
  }
  ## trim down to only what you need
  output <- output[1:i]
  time <- time[1:i]
  
  return(list(time=time,traits=output))
}

## This is working! But it is a bit slow. We can speed it up by improving the heritability a bit and increasing the initial variation
bs <- ds <- 0.01
slope <- 0.1
traitcv <- 0.3
h2 <- 0.75
tmax <- 100

set.seed(123210)
seeds <- floor(runif(40,1,1e6))
library(parallel)
mclapply(seeds,
         function(s) logistic_GEM_cull(s, dt=0.1, tmax, N0, traitmean, traitcv, h2, bs, ds, slope,X=1000,Ncull=40),
         mc.cores=8                     
) -> out3

Ntraj <- lapply(out3, function(o) (lapply(o[[2]],length) %>% unlist)[1:1001]) %>% unlist %>% matrix(., ncol=length(out3), byrow=F)
btraj <- lapply(out3, function(o) (lapply(o[[2]],mean) %>% unlist)[1:1001]) %>% unlist %>% matrix(., ncol=length(out3), byrow=F)
time.seq <- seq(0, tmax, 0.1)

## plot the mean trajectories across simulations:
par(mfrow=c(2,1), mar=c(4,4,0.5,0.5), oma=rep(0.5,4))
plot.new() 
plot.window(xlim=c(0,tmax), ylim=c(0,150))
axis(1);axis(2);box('plot')
mtext(side=1, "Time", line=3)
mtext(side=2, "Population size", line=2.5)
for (i in 1:length(out3)) lines(time.seq, Ntraj[,i], col=gray(0.7))
lines(time.seq, apply(Ntraj, 1, mean), lwd=2)

plot.new() 
plot.window(xlim=c(0,tmax), ylim=c(1,6))
axis(1);axis(2);box('plot')
mtext(side=1, "Time", line=3)
mtext(side=2, "Mean(bmax)", line=2.5)
for (i in 1:length(out3)) lines(time.seq, btraj[,i], col=gray(0.7))
lines(time.seq, apply(btraj, 1, mean), lwd=2)
abline(h=2.43961, lty=2, col=2)

