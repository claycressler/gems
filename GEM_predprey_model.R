pick_individuals <- function(N0, traitmean, traitsd) {
    mu <- log(traitmean^2 / sqrt((traitsd)^2+traitmean^2))
    sigma <- sqrt(log(traitsd^2 + 1))
    ## record this initial distribution in the output
    return(rlnorm(N0, meanlog=mu, sdlog=sigma))
}


## need to specify the initial system state (initialstate), the initial trait mean and CV, the "heritability", and a vector of parameter values
predprey_GEM <- function(seed, tmax, dt, initialstate, traitmean, traitcv, h2, params, report=FALSE) {
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
        if (report) {
            print(t)
            print(length(currenttrait))
        }
        if (length(currenttrait) > 0) { ## as long as the evolving population size is > 0, pick an individual
            ind <- sample(1:length(currenttrait),1)
            trait <- currenttrait[ind] ## the trait is the prey birth rate
            C <- currentstate["C"]
            P <- currentstate["P"]
            ## set up rates for each possible event
            Cbirth <- trait*C ## prey birth
            Cdeath <- d*C + (trait-d)/K*C*C + (a0*trait^2)*C*P ## prey death
            Pbirth <- e*(a0*trait^2)*C*P ## predator birth
            Pdeath <- m*P ## preator death

            ## wheel of fortune
            events <- c(Cbirth,Cdeath,Pbirth,Pdeath)
            wheel <- cumsum(events)/sum(events)
            ## which event happens?
            event <- min(which(runif(1) < wheel))
            if (report) {
                print(wheel)
                print(paste("event is", event))
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
                if (newtrait > 10) {
                    print(paste("new trait is", newtrait))
                    print(paste("parent trait is", trait_parent))
                    print(paste("currenttrait sd is", sd(currenttrait)))
                    print(paste("offpsring std is", off_std))
                }
                currenttrait <- c(currenttrait, newtrait)
            }
            else if (event==2) { ## prey death
                currentstate["C"] <- currentstate["C"]-1
                currenttrait <- currenttrait[-ind]
            }
            else if (event==3) ## predator birth
                currentstate["P"] <- currentstate["P"] + 1

            else ## predator death
                currentstate["P"] <- currentstate["P"] - 1

            ## advance time
            t <- t + exp(-1/sum(events))/sum(events)

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
    }
    return(list(states=stateoutput, traits=traitoutput))
}

## First, do some simulations that are just a Gillespie algorithm for a non-evolving predator-prey system.
## For b=15, the equilibrium of this non-evolving system is C = 69, P = 18
params <- c(K=2000, a0=1.25e-4, d=1, e=5, m=2.5)
seed <- 1234
tmax <- 10
dt <- 0.1
initialstate <- c(C=1000, P=1000)
traitmean <- 2
traitcv <- 0
h2 <- 1

## this should have only mild oscillations, as it is initiated at its equilibrium
output <- predprey_GEM(seed=1234, tmax=10, dt=0.1, initialstate=c(C=1000, P=1000), traitmean=2, traitcv=0, h2=1, params=params, report=FALSE)
par(mfrow=c(1,2))
plot(unlist(lapply(output$states, function(x) x[1])), type='l', ylab="Prey abundance")
plot(unlist(lapply(output$states, function(x) x[2])), type='l', ylab="Predator abundance")


## Compare against a GillespieSSA simulated using that package
library(GillespieSSA)
parms <- c(b=2, d=1, g=2.5, K=2000, c=5, a=5e-04 )
x0   <- c(N=1000, P=1000)
nu   <- matrix(c(+1, -1, 0,  0,
                  0,  0, 1, -1),nrow=2,byrow=T)
a    <- c("b*N", "d*N+(b-d)/K*N*N + a*P*N","c*a*P*N", "g*P")
tmax <- 10
out <- ssa(x0,a,nu,parms,tmax)
quartz()
par(mfrow=c(1,2))
plot(out$data[,c(1,2)], type='l', ylab="Prey abundance")
plot(out$data[,c(1,3)], type='l', ylab="Predator abundance")


## Comparison looks good!
## Now simulate trait evolution by moving away from the equilibrium. Let's imagine starting with a prey population that hasn't encountered a high predator population before, so there are a lot of prey, not as many predators, and the prey has a birth rate that is too high
output <- predprey_GEM(seed=1234, tmax=40, dt=0.1, initialstate=c(C=1200, P=800), traitmean=2.5, traitcv=0.1, h2=0.9, params=params, report=FALSE)
## where it should be going is C=1000, P=1000, b=2
par(mfrow=c(2,2))
plot(unlist(lapply(output$states, function(x) x[1])), type='l', ylab="Prey abundance")
plot(unlist(lapply(output$states, function(x) x[2])), type='l', ylab="Predator abundance")
plot(unlist(lapply(output$traits, median)), type='l', ylab="Trait median", ylim=c(0, 100))
plot(unlist(lapply(output$traits, function(x) log10(var(x)))), type='l', ylab="log10(Trait variance)", ylim=c(0, 10))



## Try a new version. In this version, I've drastically simplifed the choosing of traits by assuming that CV stays fixed and h2 has no effect on anything.
predprey_GEM2 <- function(seed, tmax, dt, initialstate, traitmean, traitcv, h2, params, report=FALSE) {
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
        if (report) {
            print(t)
            print(length(currenttrait))
        }
        if (length(currenttrait) > 0) { ## as long as the evolving population size is > 0, pick an individual
            ind <- sample(1:length(currenttrait),1)
            trait <- currenttrait[ind] ## the trait is the prey birth rate
            C <- currentstate["C"]
            P <- currentstate["P"]
            ## set up rates for each possible event
            Cbirth <- trait*C ## prey birth
            Cdeath <- d*C + (trait-d)/K*C*C + (a0*trait^2)*C*P ## prey death
            Pbirth <- e*(a0*trait^2)*C*P ## predator birth
            Pdeath <- m*P ## preator death

            ## wheel of fortune
            events <- c(Cbirth,Cdeath,Pbirth,Pdeath)
            wheel <- cumsum(events)/sum(events)
            ## which event happens?
            event <- min(which(runif(1) < wheel))
            if (event==1) { ## prey birth
                currentstate["C"] <- currentstate["C"]+1
                ## give this new individual a trait
                trait_parent <- trait
                off_std <- trait*traitcv
                newtrait <- pick_individuals(1, trait_parent, off_std)
                currenttrait <- c(currenttrait, newtrait)
            }
            else if (event==2) { ## prey death
                currentstate["C"] <- currentstate["C"]-1
                currenttrait <- currenttrait[-ind]
            }
            else if (event==3) ## predator birth
                currentstate["P"] <- currentstate["P"] + 1

            else ## predator death
                currentstate["P"] <- currentstate["P"] - 1

            ## advance time
            t <- t + exp(-1/sum(events))/sum(events)

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
    }
    return(list(states=stateoutput, traits=traitoutput))
}

## simulate the qg model
qg_model <- function(t, y, pars) {
    K <- unname(pars["K"]) ## prey carrying capacity
    d <- unname(pars["d"]) ## prey minimum death rate
    a0 <- unname(pars["a0"]) ## predator attack rate parameter
    e <- unname(pars["e"]) ## predator conversion efficiency
    m <- unname(pars["m"]) ## predator natural death rate
    V <- 0.2

    C <- y[1]
    P <- y[2]
    trait <- y[3]
    Cbirth <- trait*C ## prey birth
    Cdeath <- d*C + (trait-d)/K*C*C + (a0*trait^2)*C*P ## prey death
    Pbirth <- e*(a0*trait^2)*C*P ## predator birth
    Pdeath <- m*P ## preator death

    dCdt <- Cbirth-Cdeath
    dPdt <- Pbirth-Pdeath
    dbdt <- V*(1-C/K - 2 * a0 * trait * P)

    list(c(dCdt,dPdt,dbdt))
}
library(deSolve)
qg <- ode(y=c(C=1400,P=700,b=3.5), time=seq(0,40,0.1), func=qg_model, parms=params)

output2 <- predprey_GEM2(seed=1234, tmax=40, dt=0.1, initialstate=c(C=1400, P=700), traitmean=3.5, traitcv=0.01, h2=0.9, params=params, report=FALSE)
par(mfrow=c(2,2))
plot(unlist(lapply(output2$states, function(x) x[1])), type='l', ylab="Prey abundance")
lines(qg[,2], col=2)
plot(unlist(lapply(output2$states, function(x) x[2])), type='l', ylab="Predator abundance")
lines(qg[,3], col=2)
plot(unlist(lapply(output2$traits, median)), type='l', ylab="Trait median", ylim=c(0, 10))
lines(qg[,4], col=2)
plot(unlist(lapply(output2$traits, var)), type='l', ylab="Trait variance", ylim=c(0, 10))
## An almost perfect match...



## okay so I'm wondering about the effects of the pick_individuals function, and whether it has weird pathologies


## If I look at the Wikipedia article on the lognormal distribution, it says that the two parameters of the lognormal distribution, mu and sigma, can be related to the mean and variance of the non-logarithmized sample values, m and v, according to the following relationships:
## mu = log(m/sqrt(1+v/m^2)) = log(m^2/(sqrt(m^2+v))
## sigma = sqrt(log(1+v/m^2))

## If I modify the pick_individuals function, what will I get?
## Here's the original version
pick_individuals <- function(N0, traitmean, traitsd) {
    mu <- log(traitmean^2 / sqrt((traitsd)^2+traitmean^2))
    sigma <- sqrt(log(traitsd^2 + 1))
    ## record this initial distribution in the output
    return(rlnorm(N0, meanlog=mu, sdlog=sigma))
}
## Here's a modified version
pick_individuals2 <- function(N0, traitmean, traitsd) {
    mu <- log(traitmean^2 / sqrt((traitsd)^2+traitmean^2))
    sigma <- sqrt(log(traitsd^2/traitmean^2 + 1))
    ## record this initial distribution in the output
    return(rlnorm(N0, meanlog=mu, sdlog=sigma))
}

dist1 <- pick_individuals(1000, 1, 0.2)
dist2 <- pick_individuals2(1000, 1, 0.2)
## These are exactly the same, which makes sense because the traitmean is 1
mean(dist1)
sd(dist1)
mean(dist2)
sd(dist2)


## What if it's not
dist1 <- pick_individuals(1000, 2,  0.2)
dist2 <- pick_individuals2(1000, 2, 0.2)
## These are exactly the same, which makes sense because the traitmean is 1
mean(dist1) ## the mean is close to 2
sd(dist1) ## but the sd is 0.4, not 0.2
mean(dist2) ## here the mean is close to 2
sd(dist2) ## and the sd is 0.2, not 0.4

## Check this new version
par(mfrow=c(1,2))
plot.new()
plot.window(xlim=c(0, 2), ylim=c(0, 3))
axis(1); axis(2); box('plot')
for (i in 1:12)
    lines(seq(0.1,2,0.1), sapply(seq(0.1, 2, 0.1), function(m) pick_individuals(1000, m, m*0.05*i) %>% mean), lwd=2, col=i)
legend(x='topleft', paste("CV =", seq(0.05,0.6,0.05)), col=1:6, lwd=1.5, bty='n')

plot.new()
plot.window(xlim=c(0, 2), ylim=c(0, 3))
axis(1); axis(2); box('plot')
for (i in 1:12)
    lines(seq(0.1,2,0.1), sapply(seq(0.1, 2, 0.1), function(m) pick_individuals2(1000, m, m*0.05*i) %>% mean), lwd=2, col=i)
legend(x='topleft', paste("CV =", seq(0.05,0.6,0.05)), col=1:6, lwd=1.5, bty='n')


## Redo the above simulation but with the corrected version of pick_individuals
pick_individuals <- function(N0, traitmean, traitsd) {
    mu <- log(traitmean^2 / sqrt((traitsd)^2+traitmean^2))
    sigma <- sqrt(log(traitsd^2/traitmean^2 + 1))
    ## record this initial distribution in the output
    return(rlnorm(N0, meanlog=mu, sdlog=sigma))
}
## This fixed the problem! Thank god!
output <- predprey_GEM(seed=1234, tmax=40, dt=0.1, initialstate=c(C=1200, P=800), traitmean=2.5, traitcv=0.1, h2=0.9, params=params, report=FALSE)
## where it should be going is C=1000, P=1000, b=2
par(mfrow=c(2,2))
plot(unlist(lapply(output$states, function(x) x[1])), type='l', ylab="Prey abundance")
plot(unlist(lapply(output$states, function(x) x[2])), type='l', ylab="Predator abundance")
plot(unlist(lapply(output$traits, median)), type='l', ylab="Trait median")
plot(unlist(lapply(output$traits, function(x) var(x))), type='l', ylab="Trait variance")

