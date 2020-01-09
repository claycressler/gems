library(magrittr)
library(parallel)

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

#sirout <- sir_GEM(1287349876, tmax=200, dt=0.1, initialstate=c(S=100,I=50), traitmean=0.3, traitcv=0.2, h2=0.85, params=c(b = 2, m = 1.7, B0 = 0.05, h = 0.1))


for (i in 1:5) {
    print(i)
    params <- c(b = 2, m = 1.7, B0 = seq(0.02,0.1,0.02)[i], h = 0.1)
    set.seed(122134512334)
    seeds <- runif(100, 1, 100000) %>% floor
    mclapply(seeds,
             function(x)
                 sir_GEM(seed=x, tmax=200, dt=0.1, initialstate=c(S=100,I=50), traitmean=0.2, traitcv=0.2, h2=0.85, params=params),
             mc.cores=15
             ) -> out
    saveRDS(out, file=paste0("sir_GEM_b=2_m=1.7_B0=",unname(params["B0"]),"_h=0.1.RDS"))
}





