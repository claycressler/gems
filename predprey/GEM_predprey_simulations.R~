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

## Ignore all of the text between the braces in favor of running the content outside the braces
run <- FALSE
if (run) {
## All of the parameters here produce the same trait equilibrium (b=2) but the prey and predator equilibria get smaller as a0 increases.
set.seed(1234)
seeds <- runif(40, 1, 1e6) %>% floor
## at a0=0.0005, the C equilibrium is 400
## at a0=0.005, the C equilibrium is 40
params <- c(K=500, a0=0.001, e=5, m=4, d=1)
mclapply(seeds,
         function(x)
             predprey_GEM(seed=x, tmax=50, dt=0.1, initialstate=c(C=300,P=100), traitmean=4, traitcv=0.2, h2=0.85, params=params, report=FALSE),
         mc.cores=15
         ) -> out
saveRDS(out, file="predprey_GEM_K=500_a0=0.001_e=5_m=4_d=1.RDS")


## here the final equilibrium is much smaller - need to start pretty close in order to *not* get everything crashing to extinction
for (i in 1:8) {
    print(i)
    this.a0 <- 0.0005*i
    set.seed(1234)
    seeds <- runif(max(40, 40*i/2), 1, 1e6) %>% floor
    params <- c(K=500, a0=this.a0, e=5, m=4, d=1)
    mclapply(seeds,
             function(x)
                 predprey_GEM(seed=x, tmax=200, dt=0.1, initialstate=c(C=300/i,P=100/i), traitmean=4, traitcv=0.2, h2=0.85, params=params, report=FALSE),
             mc.cores=15
             ) -> out
    saveRDS(out, file=paste0("predprey_GEM_K=500_a0=",this.a0,"_e=5_m=4_d=1.RDS"))
}



## I want to choose some parameter values that should promote a bit less tendency to oscillate.
## Given the model I am using, with a = a0*b^2
## params <- c(K=500, a0=0.001, e=5, m=4, d=1)
## Oscillations should be smaller the closer e*a0*b^2*K is to m.
## Given an ESS b of 2*d = 2, and an initial b of 4:

## This is why we get big-ass oscillations: these are huge
with(as.data.frame(t(as.data.frame(params))), e*a0*2^2*K-m) ## = 6
with(as.data.frame(t(as.data.frame(params))), e*a0*4^2*K-m) ## = 36


## Try this:
params <- c(K=250, a0=0.001, e=5, m=4, d=1)
with(as.data.frame(t(as.data.frame(params))), e*a0*2^2*K-m) ## = 1
with(as.data.frame(t(as.data.frame(params))), e*a0*4^2*K-m) ## = 16

## here's a nice set that get close to the equilibrium
## the higher the prey or the lower the predators initially, the more likely the system is to go extinct due to stochasticity. The demographic stochasticity definitely has the effect of exciting the oscillations, though, in a way that may be difficult to prevent or control.
out <- predprey_GEM(seed=1357, tmax=50, dt=0.1, initialstate=c(C=100,P=25), traitmean=4, traitcv=0.2, h2=0.85, params=params, report=FALSE)
par(mfrow=c(1,2))
plot(seq(0,50,0.1), matrix(unlist(out$states), ncol=2, byrow=TRUE)[,1], type='l')
plot(seq(0,50,0.1), matrix(unlist(out$states), ncol=2, byrow=TRUE)[,2], type='l')



## let's systematically try to figure out how to minimize the amplitude of fluctuations
## do this based on the equilibrium predator and prey at the ESS, and knowing that the predator should start below this value and the prey should start above it
## obviously the closer everything is to the equilibrium, the lower the amplitude of fluctuations
for (i in seq(0.1,0.8,0.1)) {
    ## current parameter set
    params <- c(K=250, a0=0.001, e=5, m=4, d=1)
    ## equilibrium prey given this parameter set
    Ceq <- with(as.data.frame(t(as.data.frame(params))), m/(4*a0*d^2*e))
    Peq <- with(as.data.frame(t(as.data.frame(params))), (4*a0*d^2*e*K-m)/(16*a0^2*d^3*e*K))
    for (j in 1:10) {
        out <- predprey_GEM(seed=floor(runif(1,0,100000)), tmax=50, dt=0.1, initialstate=c(C=(1+i)*Ceq,P=i*Peq), traitmean=4, traitcv=0.2, h2=0.85, params=params, report=FALSE)
        print(paste("When C(0) =", 1+i, "*Ceq and P(0) =", i, "*Peq, prey variability is",
                    var(matrix(unlist(out$states), ncol=2, byrow=TRUE)[,1]),
                    "and predator variability is",
                    var(matrix(unlist(out$states), ncol=2, byrow=TRUE)[,2])))
    }
}

## try again just assuming that the predator abundance is constant and low (0.1*Peq)
for (i in seq(0.1,0.8,0.1)) {
    ## current parameter set
    params <- c(K=250, a0=0.001, e=5, m=4, d=1)
    ## equilibrium prey given this parameter set
    Ceq <- with(as.data.frame(t(as.data.frame(params))), m/(4*a0*d^2*e))
    Peq <- with(as.data.frame(t(as.data.frame(params))), (4*a0*d^2*e*K-m)/(16*a0^2*d^3*e*K))
    Cvars <- rep(0,6)
    Pvars <- rep(0,6)
    for (j in 1:6) {
        seed <- floor(runif(1,0,100000))
        print(seed)
        out <- predprey_GEM(seed=seed, tmax=50, dt=0.1, initialstate=c(C=(1+i)*Ceq,P=0.1*Peq), traitmean=4, traitcv=0.2, h2=0.85, params=params, report=FALSE)
        Cvars[j] <- var(matrix(unlist(out$states), ncol=2, byrow=TRUE)[,1])
        Pvars[j] <- var(matrix(unlist(out$states), ncol=2, byrow=TRUE)[,2])
    }
    print(paste0("When C(0) = ", 1+i, "*Ceq and P(0) = ", 0.1, "*Peq, prey variability is ",
                 round(mean(Cvars)), " and predator variability is ", round(mean(Pvars))))
}


## try again just assuming that the prey abundance is at the predator-free carrying capacity, varying the initial predator abundance to see how that affects the strength of oscillations. Presumably, the smaller the predator population, the larger the magnitude of oscillations.
for (i in seq(0.2,0.8,0.2)) {
    ## current parameter set
    params <- c(K=250, a0=0.001, e=5, m=4, d=1)
    ## equilibrium prey given this parameter set
    Ceq <- with(as.data.frame(t(as.data.frame(params))), m/(4*a0*d^2*e))
    Peq <- with(as.data.frame(t(as.data.frame(params))), (4*a0*d^2*e*K-m)/(16*a0^2*d^3*e*K))
    Cvars <- rep(0,6)
    Pvars <- rep(0,6)
    noExtinct <- 0
    for (j in 1:10) {
        seed <- floor(runif(1,0,100000))
        print(seed)
        out <- predprey_GEM(seed=seed, tmax=50, dt=0.1, initialstate=c(C=250,P=i*Peq), traitmean=4, traitcv=0.2, h2=0.85, params=params, report=FALSE)
        if (length(out$states) < 501) noExtinct <- noExtinct+1
        Cvars[j] <- var(matrix(unlist(out$states), ncol=2, byrow=TRUE)[,1])
        Pvars[j] <- var(matrix(unlist(out$states), ncol=2, byrow=TRUE)[,2])
    }
    print(paste0("When C(0) = 250 and P(0) = ", i, "*Peq, prey variability is ",
                 round(mean(Cvars)), ".  Predator variability is ", round(mean(Pvars)),
                 ". The number of predator extinctions was ", noExtinct))
}


## Okay, now what I need to do is vary the carrying capacity of the prey (to reduce the initial prey starting sizes) and other parameters such that the final prey equilibrium gets smaller, but the tendency to oscillate remains manageable.
## With the current parameter set, the tendency to oscillate is given by the magnitude of the following quantity (assuming b = 4, the initial b value)
params <- c(K=250, a0=0.001, e=3.5, m=2, d=1)
with(as.data.frame(t(as.data.frame(params))), e*a0*4^2*K-m) ## = 12
## at the equilibrium, the tendency to oscillate is much smaller
with(as.data.frame(t(as.data.frame(params))), e*a0*2^2*K-m) ## = 1.5
## the prey equilibrium is
with(as.data.frame(t(as.data.frame(params))), m/(4*a0*d^2*e)) ## 142.8
## Try it out!
Ceq <- with(as.data.frame(t(as.data.frame(params))), m/(4*a0*d^2*e))
Peq <- with(as.data.frame(t(as.data.frame(params))), (4*a0*d^2*e*K-m)/(16*a0^2*d^3*e*K))
Cvars <- rep(0,6)
Pvars <- rep(0,6)
noExtinct <- 0
for (j in 1:10) {
    seed <- floor(runif(1,0,100000))
    print(seed)
    out <- predprey_GEM(seed=seed, tmax=50, dt=0.1, initialstate=c(C=unname(params["K"]),P=floor(0.2*Peq)), traitmean=4, traitcv=0.2, h2=0.85, params=params, report=FALSE)
    if (length(out$states) < 501) noExtinct <- noExtinct+1
    Cvars[j] <- var(matrix(unlist(out$states), ncol=2, byrow=TRUE)[,1])
    Pvars[j] <- var(matrix(unlist(out$states), ncol=2, byrow=TRUE)[,2])
}
round(mean(Cvars))
round(mean(Pvars))
noExtinct


params <- c(K=200, a0=0.001, e=2.5, m=1, d=1)
with(as.data.frame(t(as.data.frame(params))), e*a0*4^2*K-m) ## = 7
## at the equilibrium, the tendency to oscillate is much smaller
with(as.data.frame(t(as.data.frame(params))), e*a0*2^2*K-m) ## = 1
## the prey equilibrium is
with(as.data.frame(t(as.data.frame(params))), m/(4*a0*d^2*e)) ## 100
Ceq <- with(as.data.frame(t(as.data.frame(params))), m/(4*a0*d^2*e))
Peq <- with(as.data.frame(t(as.data.frame(params))), (4*a0*d^2*e*K-m)/(16*a0^2*d^3*e*K))
Cvars <- rep(0,6)
Pvars <- rep(0,6)
noExtinct <- 0
for (j in 1:10) {
    seed <- floor(runif(1,0,100000))
    print(seed)
    out <- predprey_GEM(seed=seed, tmax=50, dt=0.1, initialstate=c(C=unname(params["K"]),P=0.2*Peq), traitmean=4, traitcv=0.2, h2=0.85, params=params, report=FALSE)
    if (length(out$states) < 501) noExtinct <- noExtinct+1
    Cvars[j] <- var(matrix(unlist(out$states), ncol=2, byrow=TRUE)[,1])
    Pvars[j] <- var(matrix(unlist(out$states), ncol=2, byrow=TRUE)[,2])
}
round(mean(Cvars))
round(mean(Pvars))
noExtinct


params <- c(K=150, a0=0.001, e=1.5, m=0.5, d=1)
with(as.data.frame(t(as.data.frame(params))), e*a0*4^2*K-m) ## = 3.1
## at the equilibrium, the tendency to oscillate is much smaller
with(as.data.frame(t(as.data.frame(params))), e*a0*2^2*K-m) ## = 0.4
## the prey equilibrium is
with(as.data.frame(t(as.data.frame(params))), m/(4*a0*d^2*e)) ## 83
Ceq <- with(as.data.frame(t(as.data.frame(params))), m/(4*a0*d^2*e))
Peq <- with(as.data.frame(t(as.data.frame(params))), (4*a0*d^2*e*K-m)/(16*a0^2*d^3*e*K))
noExtinct <- 0
for (j in 1:20) {
    seed <- floor(runif(1,0,100000))
    print(seed)
    out <- predprey_GEM(seed=seed, tmax=50, dt=0.1, initialstate=c(C=unname(params["K"]),P=0.2*Peq), traitmean=4, traitcv=0.2, h2=0.85, params=params, report=FALSE)
    if (length(out$states) < 501) noExtinct <- noExtinct+1
}
noExtinct ## 9, so assume about 50% of simulations will go extinct


params <- c(K=100, a0=0.001, e=0.5, m=0.1, d=1)
with(as.data.frame(t(as.data.frame(params))), e*a0*4^2*K-m) ## = 0.7
## at the equilibrium, the tendency to oscillate is much smaller
with(as.data.frame(t(as.data.frame(params))), e*a0*2^2*K-m) ## = 0.1
## the prey equilibrium is
with(as.data.frame(t(as.data.frame(params))), m/(4*a0*d^2*e)) ## 50
Ceq <- with(as.data.frame(t(as.data.frame(params))), m/(4*a0*d^2*e))
Peq <- with(as.data.frame(t(as.data.frame(params))), (4*a0*d^2*e*K-m)/(16*a0^2*d^3*e*K))
noExtinct <- 0
for (j in 1:20) {
    seed <- floor(runif(1,0,100000))
    print(seed)
    out <- predprey_GEM(seed=seed, tmax=50, dt=0.1, initialstate=c(C=unname(params["K"]),P=0.2*Peq), traitmean=4, traitcv=0.2, h2=0.85, params=params, report=FALSE)
    if (length(out$states) < 501) noExtinct <- noExtinct+1
    else {
        par(mfrow=c(2,2))
        plot(seq(0,50,0.1), matrix(unlist(out$states), ncol=2, byrow=TRUE)[,1], type='l')
        plot(seq(0,50,0.1), matrix(unlist(out$states), ncol=2, byrow=TRUE)[,2], type='l')
        plot(seq(0,50,0.1), unlist(lapply(out$traits, median)), type='l')
        plot(seq(0,50,0.1), unlist(lapply(out$traits, var)), type='l')
        readline("Press [enter] to continue")
    }
}
noExtinct ## 9-18, so we need a TON of replicates here


params <- c(K=50, a0=0.001, e=0.2, m=0.02, d=1)
with(as.data.frame(t(as.data.frame(params))), e*a0*4^2*K-m) ## = 0.14
## at the equilibrium, the tendency to oscillate is much smaller
with(as.data.frame(t(as.data.frame(params))), e*a0*2^2*K-m) ## = 0.02
## the prey equilibrium is
with(as.data.frame(t(as.data.frame(params))), m/(4*a0*d^2*e)) ## 25
Ceq <- with(as.data.frame(t(as.data.frame(params))), m/(4*a0*d^2*e))
Peq <- with(as.data.frame(t(as.data.frame(params))), (4*a0*d^2*e*K-m)/(16*a0^2*d^3*e*K))
noExtinct <- 0
for (j in 1:20) {
    seed <- floor(runif(1,0,100000))
    print(seed)
    out <- predprey_GEM(seed=seed, tmax=50, dt=0.1, initialstate=c(C=unname(params["K"]),P=floor(0.2*Peq)), traitmean=4, traitcv=0.2, h2=0.85, params=params, report=FALSE)
    if (length(out$states) < 501) noExtinct <- noExtinct+1
    else {
        par(mfrow=c(2,2))
        plot(seq(0,50,0.1), matrix(unlist(out$states), ncol=2, byrow=TRUE)[,1], type='l')
        plot(seq(0,50,0.1), matrix(unlist(out$states), ncol=2, byrow=TRUE)[,2], type='l')
        plot(seq(0,50,0.1), unlist(lapply(out$traits, median)), type='l')
        plot(seq(0,50,0.1), unlist(lapply(out$traits, var)), type='l')
        readline("Press [enter] to continue")
    }

}
noExtinct ## 9 extinctions
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
    seeds <- runif(100, 1, 100000) %>% floor
    mclapply(seeds,
             function(x)
                 predprey_GEM(seed=x, tmax=200, dt=0.1, initialstate=c(C=unname(params["K"]),P=floor(0.2*Peq)), traitmean=4, traitcv=0.2, h2=0.85, params=params, report=FALSE),
             mc.cores=15
             ) -> out
    saveRDS(out, file=paste0("predprey_GEM_K=",unname(params["K"]),"_a0=0.001_e=",unname(params["e"]),"_m=",unname(params["m"]),"_d=1.RDS"))
}
