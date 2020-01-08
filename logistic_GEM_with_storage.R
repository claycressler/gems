library(plyr)
pick_individuals <- function(N0, traitmean, traitsd) {
    mu <- log(traitmean^2 / sqrt((traitsd)^2+traitmean^2))
    sigma <- sqrt(log(traitsd^2/traitmean^2 + 1))
    ## record this initial distribution in the output
    return(rlnorm(N0, meanlog=mu, sdlog=sigma))
}

## Note - this can only be used for situations where the population size is not too large, or where you are not simulating many, many replicates.
## otherwise the storage alone will be enough to crash R
logistic_GEM_storage <- function(seed, tmax, N0, traitmean, traitcv, h2, bs, ds, slope) {

    ## set RNG seed
    set.seed(seed)
    ## initialize
    t <- 0 ## current time
    ## storage for time series
    trait=rep(0,5e5)
    tBirth=rep(0,5e5)
    tDeath=rep(0,5e5)
    nOff=rep(0,5e5)
    ancTrait=rep(0,5e5)

    ## draw initial trait distribution from lognormal with meanlog and sdlog given by traitmean and traitcv
    currentstate <- pick_individuals(N0, traitmean, traitsd=traitmean*traitcv)
    currentID <- 1:N0

    ## add these individuals to output
    trait[1:N0] <- currentstate

    ## set the next id to use
    ID <- N0+1
    while (t < tmax & length(currentstate) > 0) {
        ## as long as population size is > 0 pick an individual
        ## pick the id of the individual
        this.ind <- sample(1:length(currentstate), 1)
        ## get the trait of the individual with that id
        this.trait <- currentstate[this.ind]
        ## figure out which index this is in the output vectors
        this.parent <- currentID[this.ind]
        ## current population size
        N <- length(currentstate)
        ## set up rates for each possible event
        ## birth
        brate <- (this.trait - bs*N)*N
        ## death (specified by tradeoff between birth and death)
        d <- slope*this.trait^2
        drate <- (d + ds*N)*N
        events <- c(brate,drate)

        ## wheel of fortune
        wheel <- cumsum(events)/sum(events)
        ## which event happens?
        event <- min(which(runif(1) < wheel))
        ## what time did the event occur?
        t <- t + exp(-1/sum(events))/sum(events)

        ## if there is a birth
        if (event==1) { ## birth
            ## compute the trait of the offspring
            trait_parent <- (1-h2)*mean(currentstate) + h2*this.trait
            ## offspring trait distribution standard deviation
            ## if there is more than one individual, this will work
            off_std <- sqrt(1-h2^2)*((1-h2)*sd(trait[1:N0])+h2*sd(currentstate))
            ## if there is only one individual left in the population, set off_std = 0
            if (is.na(off_std)) off_std <- 0
            ## generate the new individual
            newtrait <- pick_individuals(1, trait_parent, off_std)
            ## add this individual to the current state
            currentstate <- c(currentstate, newtrait)
            currentID <- c(currentID, ID)
            ## add this individual to the overall output
            trait[ID] <- newtrait
            tBirth[ID] <- t
            ancTrait[ID] <- this.trait
            ## increase the nOff for the parent
            nOff[this.parent] <- nOff[this.parent] + 1
            ## update the id
            ID <- ID+1
        }
        if (event==2) {  ## death
            ## remove individual from currentstate by its id
            currentstate <- currentstate[-this.ind]
            currentID <- currentID[-this.ind]
            ## set the time of death for the individual
            tDeath[this.parent] <- t
        }

        ## if ID > allocation, allocate more storage to every vector
        if (ID > length(trait)) {
            print("adding more storage")
            trait <- c(trait, rep(0,1000))
            tBirth <- c(tBirth, rep(0,1000))
            tDeath <- c(tDeath, rep(0,1000))
            nOff <- c(nOff, rep(0,1000))
            ancTrait <- rep(ancTrait, rep(0,1000))
        }

    }
    ## prune the output dataframe to the number of inds actually present in the data
    output <- data.frame(trait=trait[1:(ID)],
                         ancTrait=ancTrait[1:(ID)],
                         nOff=nOff[1:(ID)],
                         tBirth=tBirth[1:(ID)],
                         tDeath=tDeath[1:(ID)])
    ## if the last event was a birth, then there is an additional position that is all 0 - remove it
    if(event==1) output <- output[-ID,]
    ## add some extra information
    mutate(output,
           lifespan=tDeath-tBirth,
           R0=nOff/lifespan) -> output

    return(output)
}

