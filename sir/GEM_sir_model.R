
pick_individuals <- function(N0, traitmean, traitsd) {
    mu <- log(traitmean^2 / sqrt((traitsd)^2+traitmean^2))
    sigma <- sqrt(log(traitsd^2/traitmean^2 + 1))
    ## record this initial distribution in the output
    return(rlnorm(N0, meanlog=mu, sdlog=sigma))
}


## need to specify the initial system state (initialstate), the initial trait mean and CV, the "heritability", and a vector of parameter values
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


qg_sir_model <- function(t, y, pars) {
    b <- pars["b"]
    m <- pars["m"]
    B0 <- pars["B0"]
    h <- pars["h"]
    V <- pars["V"]

    S <- y[1]
    I <- y[2]
    v <- y[3]
    beta <- B0*v/(h+v)
    dSdt <- b*(S+I) - beta*S*I - m*S
    dIdt <- beta*S*I - (m+v)*I
    dvdt <- V * (B0*h*S/((h+v)^2) - 1)
    list(c(dSdt, dIdt, dvdt))
}
library(deSolve)

## Compare the SIR GEM to the QG model
sirout <- sir_GEM(1287349876, tmax=200, dt=0.1, initialstate=c(S=60,I=20), traitmean=0.2, traitcv=0.2, h2=0.85, params=c(b = 2, m = 1.7, B0 = 0.05, h = 0.1))

out0 <- ode(y = c(S=100, I=50, v=0.3), times=seq(0, 200, 0.1), func=qg_sir_model, parms=c(b = 2, m = 1.7, B0 = 0.05, h = 0.1, V=0.001))

par(mfrow=c(2,2))
lapply(sirout$states, function(x) x[1]) %>% unlist %>% plot(., type='l', ylab="Susc abund")
lines(out0[,2], col=2)
lapply(sirout$states, function(x) x[2]) %>% unlist %>% plot(., type='l', ylab="Inf abund")
lines(out0[,3], col=2)
lapply(sirout$traits, median) %>% unlist %>% plot(., type='l', ylab="Median trait")
lines(out0[,4], col=2)
lapply(sirout$traits, var) %>% unlist %>% plot(., type='l', ylab="Var trait")


sirout <- sir_GEM(1287349876, tmax=200, dt=0.1, initialstate=c(S=60,I=20), traitmean=0.6, traitcv=0.2, h2=0.85, params=c(b = 2, m = 1.7, B0 = 0.08, h = 0.1), report=TRUE)

out0 <- ode(y = c(S=60, I=20, v=0.6), times=seq(0, 200, 0.1), func=qg_sir_model, parms=c(b = 2, m = 1.7, B0 = 0.12, h = 0.1, V=0.001))
plot(out0[,c(1,3)], type='l')
