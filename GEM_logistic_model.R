## initialization
tmax <- 70 ## length of time series
dt <- 0.1 ## size of time steps
seed <- 1234320 ## RNG seed
N0 <- 10 ## initial population size
traitmean <- 1 ## initial mean of trait distribution
traitcv <- 0.1 ## initial coefficient of variation of trait distribution
h2 <- 0.8 ## heritability

## logistic equation parameters
bs <- 0.01
ds <- 0.01
slope <- 1/6

pick_individuals <- function(N0, traitmean, traitsd) {
    mu <- log(traitmean^2 / sqrt((traitsd)^2+traitmean^2))
    sigma <- sqrt(log(traitsd^2 + 1))
    ## record this initial distribution in the output
    return(rlnorm(N0, meanlog=mu, sdlog=sigma))
}


logistic_GEM <- function(seed, tmax, dt, N0, traitmean, traitcv, h2, bs, ds, slope) {

    ## set RNG seed
    set.seed(seed)
    ## initialize
    times <- seq(0,tmax,by=dt)
    t <- 0 ## current time
    ## storage for time series; at each timestep, record the traits of every individual in the population
    output <- vector("list", length(times))

    ## draw initial trait distribution from lognormal with meanlog and sdlog given by traitmean and traitcv
    output[[1]] <- pick_individuals(N0, traitmean, traitsd=traitmean*traitcv)
    ## also record it in a vector called 'currentstate' that records the current state of the system only
    currentstate <- output[[1]]

    i <- 1 ## how many things in the output file?
    while (t < tmax) {
        if (length(currentstate) > 0) { ## as long as population size is > 0 pick an individual
            ind <- sample(1:length(currentstate),1)
            trait <- currentstate[ind]
            N <- length(currentstate)
            ## set up rates for each possible event
            ## birth
            brate <- (trait - bs*N)*N
            ## death (specified by tradeoff between birth and death)
            d <- slope*trait^2
            drate <- (d + ds*N)*N
            events <- c(brate,drate)

            ## wheel of fortune
            wheel <- cumsum(events)/sum(events)
            ## which event happens?
            event <- min(which(runif(1) < wheel))

            if (event==1) { ## birth
                ## compute the trait of the offspring
                ## offspring trait distribution mean
                trait_parent <- (1-h2)*mean(currentstate) + h2*trait
                off_std <- sqrt(1-h2^2)*((1-h2)*sd(output[[1]])+h2*sd(currentstate))
                newtrait <- pick_individuals(1, trait_parent, off_std)
                ## add this individual to the current state vector
                currentstate <- c(currentstate, newtrait)
            }
            else  ## death
                currentstate <- currentstate[-ind]

            ## advance time
            t <- t + exp(-1/sum(events))/sum(events)

            ## record in the output if necessary
            if (t > times[i+1]) {## if time has passed the next recording time
                output[[i+1]] <- currentstate ## record the system state
                i <- i + 1 ## advance to the next timestep
            }
        }
    }
    return(output)
}

## here is an alternative version of the model where, instead of randomly choosing the individual and then creating the wheel of fortune based on its traits, instead computes a wheel of fortune for all possible events that occur to *any* individual and then uses the RNG to choose both the individual and the event
logistic_GEM2 <- function(seed, tmax, dt, N0, traitmean, traitcv, h2, bs, ds, slope) {

    ## set RNG seed
    set.seed(seed)
    ## initialize
    times <- seq(0,tmax,by=dt)
    t <- 0 ## current time
    ## storage for time series; at each timestep, record the traits of every individual in the population
    output <- vector("list", length(times))

    ## draw initial trait distribution from lognormal with meanlog and sdlog given by traitmean and traitcv
    output[[1]] <- pick_individuals(N0, traitmean, traitsd=traitmean*traitcv)
    ## also record it in a vector called 'currentstate' that records the current state of the system only
    currentstate <- output[[1]]

    i <- 1 ## how many things in the output file?
    while (t < tmax) {
        if (length(currentstate) > 0) { ## as long as population size is > 0 pick an individual
            ## wheel of fortune
            N <- length(currentstate)
            brates <- (currentstate - bs*N)*N
            d <- slope*currentstate^2
            drates <- (d + ds*N)*N
            ## wheel is set up with births first, then deaths
            events <- c(brates,drates)
            wheel <- cumsum(events)/sum(events)
            rand <- runif(1)
            ## who does the event happen to?
            ind <- ifelse(min(which(rand < wheel)) <= N,
                          min(which(rand < wheel)),
                          min(which(rand < wheel))-N)
            ## was the event a birth or death?
            event <- ifelse(rand <= sum(brates)/sum(events),
                            1, ## birth
                            2) ## death

            if (event==1) { ## birth
                ## compute the trait of the offspring
                ## offspring trait distribution mean
                trait_parent <- (1-h2)*mean(currentstate) + h2*currentstate[ind]
                off_std <- sqrt(1-h2^2)*((1-h2)*sd(output[[1]])+h2*sd(currentstate))
                newtrait <- pick_individuals(1, trait_parent, off_std)
                ## add this individual to the current state vector
                currentstate <- c(currentstate, newtrait)
            }
            else  ## death
                currentstate <- currentstate[-ind]

            ## advance time
            t <- t + exp(-1/sum(events))/sum(events)

            ## record in the output if necessary
            if (t > times[i+1]) {## if time has passed the next recording time
                output[[i+1]] <- currentstate ## record the system state
                i <- i + 1 ## advance to the next timestep
            }
        }
    }
    return(output)
}


## quantitative genetics model (here assuming an accelerating trade-off that will produce an ESS
qg_model <- function(t, y, pars) {
    bs <- pars["bs"]
    ds <- pars["ds"]
    slope <- pars["slope"]
    V <- pars["V"]

    N <- y[1]
    b <- y[2]
    d <- slope*b^2
    dNdt <- (b-bs*N)*N-(d+ds*N)*N
    dbdt <- V*(1-2*b*slope)
    list(c(dNdt,dbdt))
}


## low variation
## initialization
dt <- 0.1 ## size of time steps
N0 <- 10 ## initial population size
traitmean <- 1 ## initial mean of trait distribution
traitcv <- 0.1 ## initial coefficient of variation of trait distribution
h2 <- 0.8 ## heritability

## logistic equation parameters
bs <- 0.01
ds <- 0.01
slope <- 1/6

out1 <- logistic_GEM(23412, 100, dt, N0, traitmean, traitcv, h2, bs, ds, slope)
out2 <- logistic_GEM2(23412, 100, dt, N0, traitmean, traitcv, h2, bs, ds, slope) ## much, much slower because time advances much more slowly

par(mfrow=c(2,2))
plot(seq(0,100,dt), lapply(out1, length) %>% unlist, type='l')
plot(seq(0,100,dt), lapply(out2, length) %>% unlist, col=2, type='l')

plot(seq(0,100,dt), lapply(out1, mean) %>% unlist, type='l')
plot(seq(0,100,dt), lapply(out2, mean) %>% unlist, col=2, type='l')

traitcv <- 0.2 ## initial coefficient of variation of trait distribution
h2 <- 0.9 ## heritability
out3 <- logistic_GEM(23412, 100, dt, N0, traitmean, traitcv, h2, bs, ds, slope)
out4 <- logistic_GEM2(23412, 100, dt, N0, traitmean, traitcv, h2, bs, ds, slope) ## much, much slower because time advances much more slowly

par(mfrow=c(2,2))
plot(seq(0,100,dt), lapply(out3, length) %>% unlist, type='l')
plot(seq(0,100,dt), lapply(out4, length) %>% unlist, col=2, type='l')

plot(seq(0,100,dt), lapply(out3, mean) %>% unlist, type='l')
plot(seq(0,100,dt), lapply(out4, mean) %>% unlist, col=2, type='l')

plot(seq(0,100,dt), lapply(out2, mean) %>% unlist, col=1, type='l')
lines(seq(0,100,dt), lapply(out4, mean) %>% unlist, col=2, type='l')










results <- vector(mode='list', length=20)
for (i in 1:20) {
    print(i)
    results[[i]] <- logistic_GEM(floor(runif(1,1,1e7)), 100, dt, N0, traitmean, traitcv, h2, bs, ds, slope, tradeoff="slope*trait^2")
}
library(deSolve)
out2 <- lsoda(y=c(10,1), times=seq(0,100,0.1), func=qg_model, parms=c(bs=bs, ds=ds, slope=slope, V=traitcv))

par(mfrow=c(1,2))
## calculate mean population dynamics across the simulations
lapply(results, function(res) lapply(res, length) %>% unlist) %>% unlist %>% matrix(., ncol=20, byrow=FALSE) %>% as.data.frame -> out
out$mean <- apply(out, 1, mean)
plot.new()
plot.window(xlim=c(1,nrow(out)), ylim=range(out))
axis(1); axis(2); box('plot')
mtext(side=1, line=3, "time")
mtext(side=2, line=3, "pop size")
for (i in 1:20)
    lines(out[,i], col=gray(0.7))
lines(out$mean, col='blue', lwd=2)
lines(out2[,2], col='red', lwd=2)

## calculate mean trait dynamics across the simulations
lapply(results, function(res) lapply(res, mean) %>% unlist) %>% unlist %>% matrix(., ncol=20, byrow=FALSE) %>% as.data.frame -> out
out$mean <- apply(out, 1, mean)
plot.new()
plot.window(xlim=c(1,nrow(out)), ylim=c(0,10))
axis(1); axis(2); box('plot')
mtext(side=1, line=3, "time")
mtext(side=2, line=3, "trait mean")
for (i in 1:20)
    lines(out[,i], col=gray(0.7))
lines(out$mean, col='blue', lwd=2)
lines(out2[,3], col='red', lwd=2)
