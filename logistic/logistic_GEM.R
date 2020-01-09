pick_individuals <- function(N0, traitmean, traitsd) {
    mu <- log(traitmean^2 / sqrt(traitsd^2+traitmean^2))
    sigma <- sqrt(log(traitsd^2/traitmean^2 + 1))
    ## record this initial distribution in the output
    return(rlnorm(N0, meanlog=mu, sdlog=sigma))
}

## the original GEM
logistic_GEM <- function(seed, dt, tmax, N0, traitmean, traitcv, h2, bs, ds, slope) {

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
        events <- c(brate,drate)

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

## the many-slice GEM
logistic_GEM2 <- function(seed, dt, tmax, N0, traitmean, traitcv, h2, bs, ds, slope,preallocation=1e6) {

    ## set RNG seed
    set.seed(seed)

    ## storage for times and trait distributions
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

    while (t < tmax & length(currentstate) > 0) { ## as long as population size is > 0 pick an individual
        N <- length(currentstate)
        ## set up per-capita rates for each individual
        ## birth
        brates <- (currentstate - bs*N)
        ## death (specified by tradeoff between birth and death)
        d <- slope*currentstate^2
        drates <- (d + ds*N)
        events <- c(brates,drates)

        ## wheel of fortune
        wheel <- cumsum(events)/sum(events)
        rand <- runif(1)
        ## who does the event happen to?
        ind <- ifelse(min(which(rand < wheel)) <= N,
                      min(which(rand < wheel)),
                      min(which(rand < wheel))-N)
        ## which event happens?
        event <- ifelse(rand <= sum(brates)/sum(events),
                        1, ## birth
                        2) ## death
        ## when did this event happen?
        timestep <-  exp(-1/sum(events))/sum(events)
        t <- t + timestep

        if (event==1) { ## birth
            ## compute the trait of the offspring
            ## offspring trait distribution mean
            trait_parent <- (1-h2)*mean(currentstate) + h2*currentstate[ind]
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
        }
        else if ((floor(t/dt)*dt) > lastrecordtime) {
            i <- i + 1
            time[i] <- t
            output[[i]] <- currentstate
            lastrecordtime <- floor(t/dt)*dt
        }
    }
    ## trim down to only what you need
    output <- output[1:i]
    time <- time[1:i]

    return(list(time=time,traits=output))
}

## an another alternate version of the GEM
## for this one, the per-capita rates are summed to determine which event happens
## then an individual is chosen based on its trait using a second wheel of fortune across all individuals
logistic_GEM3 <- function(seed, dt, tmax, N0, traitmean, traitcv, h2, bs, ds, slope,preallocation=1e6) {

    ## set RNG seed
    set.seed(seed)

    ## storage for times and trait distributions
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

    while (t < tmax & length(currentstate) > 0) { ## as long as population size is > 0 pick an individual
        N <- length(currentstate)
        ## set up per-capita rates for each individual
        ## birth
        brates <- (currentstate - bs*N)
        ## death (specified by tradeoff between birth and death)
        d <- slope*currentstate^2
        drates <- (d + ds*N)
        events <- c(brates,drates)

        ## 1st wheel of fortune to choose an event
        wheel <- cumsum(events)/sum(events)
        rand <- runif(1)
        event <- ifelse(rand <= sum(brates)/sum(events),
                        1, ## birth
                        2) ## death

        ## when did this event happen?
        timestep <-  exp(-1/sum(events))/sum(events)
        t <- t + timestep

        if (event==1) { ## birth
            ## 2nd wheel of fortune to determine who gives birth
            wheel2 <- cumsum(brates)/sum(brates)
            rand2 <- runif(1)
            ind <- min(which(rand2 < wheel2))
            ## compute the trait of the offspring
            ## offspring trait distribution mean
            trait_parent <- (1-h2)*mean(currentstate) + h2*currentstate[ind]
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
        else  {## death
            ## 2nd wheel of fortune to determine who dies
            wheel2 <- cumsum(drates)/sum(drates)
            rand2 <- runif(1)
            ind <- min(which(rand2 < wheel2))
            currentstate <- currentstate[-ind]
        }

        ## record... but do it smartly, to prevent output from getting so big that it crashes R.
        ## Record every event if it moves time forward more than dt time steps, otherewise only record every dt timesteps
        if (timestep > dt) {
            i <- i + 1
            time[i] <- t
            output[[i]] <- currentstate
        }
        else if ((floor(t/dt)*dt) > lastrecordtime) {
            i <- i + 1
            time[i] <- t
            output[[i]] <- currentstate
            lastrecordtime <- floor(t/dt)*dt
        }
    }
    ## trim down to only what you need
    output <- output[1:i]
    time <- time[1:i]

    return(list(time=time,traits=output))
}


