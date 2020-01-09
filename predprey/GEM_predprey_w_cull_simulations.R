## This script defines a predator prey GEM very similar to the one in GEM_predprey_simulations.R with one addition: a step that culls the predator population to a fixed value. This introduces strong environmental stochasticity that allows us to look at TEAs (for example).
library(magrittr)
library(parallel)

pick_individuals <- function(N0, traitmean, traitsd) {
    mu <- log(traitmean^2 / sqrt((traitsd)^2+traitmean^2))
    sigma <- sqrt(log(traitsd^2/traitmean^2 + 1))
    ## record this initial distribution in the output
    return(rlnorm(N0, meanlog=mu, sdlog=sigma))
}


## need to specify the initial system state (initialstate), the initial trait mean and CV, the "heritability", and a vector of parameter values
predprey_GEM_cull <- function(seed, tmax, dt, initialstate, traitmean, traitcv, h2, params, report=FALSE) {
    ## Extract the parameters
    K <- unname(params["K"]) ## prey carrying capacity
    d <- unname(params["d"]) ## prey minimum death rate
    a0 <- unname(params["a0"]) ## predator attack rate parameter
    e <- unname(params["e"]) ## predator conversion efficiency
    m <- unname(params["m"]) ## predator natural death rate

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
    currenttrait <- pick_individuals(initialstate["C"], traitmean, traitsd=traitmean*traitcv)

    ## record the initial system and trait states
    stateoutput[[1]] <- currentstate
    traitoutput[[1]] <- currenttrait

    i <- 1 ## how many things in the output file?
    while (t < tmax) {
        if (length(currenttrait) > 0) { ## as long as the evolving population size is > 0, pick an individual
            ind <- sample(1:length(currenttrait),1)
            trait <- currenttrait[ind] ## the trait is the prey birth rate
            C <- currentstate["C"]
            P <- currentstate["P"]
            ## set up rates for each possible event
            Cbirth <- trait*C ## prey birth
            Cdeath <- d*C+(trait-d)/K*C^2 + (a0*trait^2)*C*P ## prey death
            Pbirth <- e*(a0*trait^2)*C*P ## predator birth
            Pdeath <- m*P ## preator death

            ## if Cdeath < 0, kill this individual
            if (Cdeath < 0)
                event <- 2
            else {
                ## wheel of fortune
                events <- c(Cbirth,Cdeath,Pbirth,Pdeath)
                wheel <- cumsum(events)/sum(events)
                ## which event happens?
                event <- min(which(runif(1) < wheel))
                if (report) {
                    print(t)
                    print(length(currenttrait))
                    print(currentstate)

                }
            }

            if (event==1) { ## prey birth
                currentstate["C"] <- currentstate["C"]+1
                ## give this new individual a trait
                trait_parent <- (1-h2)*mean(currenttrait) + h2*trait
                if (length(currenttrait) > 1) ## if there is more than one individual in the current population, you can compute the standard deviation of the current trait distribution
                    off_std <- sqrt(1-h2^2)*((1-h2)*sd(traitoutput[[1]])+h2*sd(currenttrait))
                ## otherwise, just set this sd=0
                else off_std <- sqrt(1-h2^2)*((1-h2)*sd(traitoutput[[1]]))
                newtrait <- pick_individuals(1, trait_parent, off_std)
                currenttrait <- c(currenttrait, newtrait)
            }
            else if (event==2) { ## prey death
                currentstate["C"] <- currentstate["C"]-1
                currenttrait <- currenttrait[-ind]
            }
            else if (event==3) ## predator birth
                currentstate["P"] <- currentstate["P"] + 1

            else if (event==4) ## predator death
                currentstate["P"] <- currentstate["P"] - 1

            else if (P==0) ## predator has gone extinct and prey will eventually reach the equilibrium
                break

            else
                print("something has gone awry")

            ## advance time
            t <- t + exp(-1/sum(events))/sum(events)

            ## cull the predator population to only 50 if it's larger than that
            ## this should create a strong TEA
            if (P > 50)
                currentstate["P"] <- 50

            ## record in the output if necessary
            if (t > times[i+1]) {## if time has passed the next recording time
                stateoutput[[i+1]] <- currentstate ## record the system state
                traitoutput[[i+1]] <- currenttrait
                i <- i + 1 ## advance to the next timestep
            }
        }
        else {
            if (report) print("Evolving population is extinct")
            break
        }
        ## if the predator has gone extinct, break the whole thing off before weird numerical simulation errors occur
        if (P==0)
            break
    }

    return(list(states=stateoutput[1:i], traits=traitoutput[1:i]))
}


## Okay, here are the parameter sets I want to run.

paramsets <- vector(mode='list', length=5)
paramsets[[1]] <- c(K=250, a0=0.001, e=3.5, m=2, d=1)
paramsets[[2]] <- c(K=200, a0=0.001, e=2.5, m=1, d=1)
paramsets[[3]] <- c(K=150, a0=0.001, e=1.5, m=0.5, d=1)
paramsets[[4]] <- c(K=100, a0=0.001, e=0.5, m=0.1, d=1)
paramsets[[5]] <- c(K=50, a0=0.001, e=0.2, m=0.02, d=1)

for (i in 1:5) {
    print(i)
    params <- paramsets[[i]]
    Peq <- with(as.data.frame(t(as.data.frame(params))), (4*a0*d^2*e*K-m)/(16*a0^2*d^3*e*K))
    set.seed(1234)
    seeds <- runif(20, 1, 100000) %>% floor
    output <- vector(mode='list', length=20)
    for (j in 12:20) {
        print(j)
        output[[j]] <- predprey_GEM_cull(seed=seeds[j], tmax=500, dt=0.1, initialstate=c(C=unname(params["K"]),P=floor(0.2*Peq)), traitmean=4, traitcv=0.2, h2=0.85, params=params, report=FALSE)
    }
    saveRDS(output, file=paste0("predprey_GEM_cull_K=",unname(params["K"]),"_a0=0.001_e=",unname(params["e"]),"_m=",unname(params["m"]),"_d=1.RDS"))
}
