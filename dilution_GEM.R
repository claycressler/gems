library(deSolve)

dilution.model <- function(t, y, p) {
    Z = y[1]
    I = y[2]
    S = y[3]
    R = y[4]

    sigma=unname(p["sigma"])
    mz=unname(p["mz"])
    a=unname(p["a"])
    epsilon=unname(p["epsilon"])
    v=unname(p["v"])
    e=unname(p["e"])
    ms=unname(p["ms"])
    r=unname(p["r"])
    K=unname(p["K"])

    dZdt = sigma*I - mz*Z - a*Z*(S + I)
    dIdt = epsilon*a*S*Z - (ms+v)*I
    dSdt = e*a*R*(S + I) - epsilon*a*S*Z - ms*S
    dRdt = r*R*(1 - R/K) - a*R*(S + I)

    return(list(c(dZdt,dIdt,dSdt,dRdt)))
}

## This set of parameters does not produce a stochastic simulation that is anything like the analytical simulation
pars <- c(sigma=100,
          mz=1,
          a=0.1,
          epsilon=1,
          v=10,
          e=0.1,
          ms=5,
          r=100,
          K=5000)
## This set of parameters produces a stochastic simulation that matches with the ODE (although the dynamics are pretty trivial - everything goes extinct except resources)
pars <- c(sigma=0.5,
          mz=0.8,
          a=0.005,
          epsilon=0.7156,
          v=0.05,
          e=0.02,
          ms=0.05,
          r=2,
          K=1500)
times <- seq(0,300,0.1)
y0 <- c(Z=100, I=0, S=400, R=500)
out <- ode(y0, times, dilution.model, pars)
tail(out)

## define key evolution values and traits
IS_resource_1 = 9 ## value of species 1 in immunological space
IS_par_A <- 10 ## area under the epsilon curve for parasite generalism
IS_par_p <- 1  ## peak epsilon for parasite
IS_host_range_mean <- 5 ## location of peak epsilon in host immunological space
cv <- 0.2

## Gillespie algorithm
num_replicates <- 10 ## number of simulations
stand_times <- seq(0, tmax)
num_time_steps <- length(stand_times)

## storage vectors
R_stand <- array(NA, dim=c(num_replicates, num_time_steps))
S_stand <- array(NA, dim=c(num_replicates, num_time_steps))
I_stand <- array(NA, dim=c(num_replicates, num_time_steps))
Z_stand <- array(NA, dim=c(num_replicates, num_time_steps))

for (i in 1:num_replicates) { ## start Gillespie algorithm
    ## proximate storage for each Gillespie step
    R <- rep(0, 1e7)
    S <- rep(0, 1e7)
    I <- rep(0, 1e7)
    Z <- rep(0, 1e7)
    t <- rep(0, 1e7)

    ## initial states
    t[1] <- 0
    R[1] <- unname(y0["R"])
    S[1] <- unname(y0["S"])
    I[1] <- unname(y0["I"])
    Z[1] <- unname(y0["Z"])

    ## create initial distribution for evolving parameter
    ## choose a random number generator seed and save it
    #seed <- as.integer(round(runif(1, 1, 1e7)))
    #saveRDS(seed, file=paste("rng_seed",
                             strsplit(as.character(Sys.time()), " ")[[1]][1],
                             "rep", i, sep="_"))
    #set.seed(seed)

    ## Trait distributions - define as many as needed for the number of evolving traits (ncol) and for the initial number of individuals in the evolving population (nrow)
    #x_dist_init <- array(NA, dim=c(Z[1], 2))
    ## create initial distribution for parameter
    #MU <- log(IS_host_p^2 / sqrt((cv*IS_host_p)^2+IS_host_p^2)) ## mean of lognormal
    #SIGMA <- sqrt(log((cv*IS_host_p)^2/IS_host_p^2 + 1)) ## stdev for lognormal
    ## set the initial trait distribution
    #x_dist_init[,1] <- rlnorm(Z[1], MU, SIGMA)

    ## distribution of trait mean
    #MU = log(IS_host_range_mean^2 / sqrt((cv*IS_host_range_mean)^2+IS_host_range_mean^2)) ## mean for lognormal
    #SIGMA = sqrt(log((cv*IS_host_range_mean)^2/IS_host_range_mean^2 + 1)) # std for lognormal
    #x_dist_init[,2] <- rlnorm(Z[1], MU, SIGMA)

    #x_dist <- x_dist_init # reset trait distribution at the start of each simulation
    #x_mean <- apply(x_dist, 2, mean) # initial mean trait
    #x_var <- apply(x_dist, 2, var) # initial mean variance

count <- 1; ## index steps while in loop
while(t[count] < tmax) { ## run until reach tmax
print(t[count])
##if (nrow(x_dist) > 0) {## as long as there is more than one individual in the evolving population, pick another individual
        #    whosnext <- sample(1:nrow(x_dist), 1)
        #    p_next <- x_dist[whosnext,1]
        #    range_mean_next <- x_dist[whosnext,2]
        #    epsilon_next <- p_next - (16*p_next^3*(IS_resource_1 - range_mean_next)^2)/(9*IS_host_A^2)
        #}

        ## Rates of each possible event given by ODEs
        ## birth rate of prey
        b_R <- unname(pars["r"]*R[count])
        ## death rate of prey
        d_R_nat <- unname(pars["r"]*R[count]^2/pars["K"])
        ## consumption rate of resource
        d_R_pred <- unname(pars["a"] * R[count] * (S[count] + I[count]))
        ## birth rate of new susceptibles
        b_S <- unname(pars["e"]*pars["a"]*R[count]*(S[count] + I[count]))
        ## death rate of susceptibles
        d_S <- unname(pars["ms"]*S[count])
        ## death rate of infecteds
        d_I <- unname((pars["ms"]+pars["v"])*I[count])
        ## transmission rate
        beta <- unname(pars["epsilon"]*pars["a"]*S[count]*Z[count]) ## replace w/epsilon_next
        ## shedding new parasite
        b_Z <- unname(pars["sigma"]*I[count])
        ## death of parasite
        d_Z_nat <- unname(pars["mz"]*Z[count])
        ## loss due to consumption
        d_Z_pred <- unname(pars["a"]*Z[count]*(S[count]+I[count]))

        events <- c(b_R, d_R_nat, d_R_pred, b_S, d_S, d_I, beta, b_Z, d_Z_nat, d_Z_pred)
        sc_events <- cumsum(events)/sum(events)
        event_index <- min(which(runif(1) < sc_events))
        if (event_index==1) {## prey birth
            R[count+1] <- R[count]+1
            S[count+1] <- S[count]
            I[count+1] <- I[count]
            Z[count+1] <- Z[count]
        }
        else if (event_index==2) { ## prey death
            R[count+1] <- R[count]-1
            S[count+1] <- S[count]
            I[count+1] <- I[count]
            Z[count+1] <- Z[count]
        }
        else if (event_index==3) { ## prey natural death
            R[count+1] <- R[count]-1
            S[count+1] <- S[count]
            I[count+1] <- I[count]
            Z[count+1] <- Z[count]
        }
        else if (event_index==4) { ## susceptible birth
            R[count+1] <- R[count]
            S[count+1] <- S[count]+1
            I[count+1] <- I[count]
            Z[count+1] <- Z[count]
        }
        else if (event_index==5) { ## susceptible death
            R[count+1] <- R[count]
            S[count+1] <- S[count]-1
            I[count+1] <- I[count]
            Z[count+1] <- Z[count]
        }
        else if (event_index==6) { ## infected death
            R[count+1] <- R[count]
            S[count+1] <- S[count]
            I[count+1] <- I[count]-1
            Z[count+1] <- Z[count]
        }
        else if (event_index==7) { ## transmission
            R[count+1] <- R[count]
            S[count+1] <- S[count]-1
            I[count+1] <- I[count]+1
            Z[count+1] <- Z[count]
        }
        else if (event_index==8) { ## birth of parasite
            R[count+1] <- R[count]
            S[count+1] <- S[count]
            I[count+1] <- I[count]
            Z[count+1] <- Z[count]+1
        }
        else if (event_index==9) { ## natural death of parasite
            R[count+1] <- R[count]
            S[count+1] <- S[count]
            I[count+1] <- I[count]
            Z[count+1] <- Z[count]-1
        }
        else { ## consumption of parasite
            R[count+1] <- R[count]
            S[count+1] <- S[count]
            I[count+1] <- I[count]
            Z[count+1] <- Z[count]-1
        }

        ## when did this event happen?
        t[count+1] <- t[count] + exp(-1/sum(events))/sum(events)
        count <- count+1
    }

plot(t[1:count], R[1:count], type='l', lwd=2, ylim=c(0,1000))
lines(t[1:count], S[1:count], type='l', lwd=2, col=2)
lines(t[1:count], I[1:count], type='l', lwd=2, col=3)
lines(t[1:count], Z[1:count], type='l', lwd=2, col=4)

lines(out[,c(1,5)])
lines(out[,c(1,4)], col=2)
lines(out[,c(1,3)], col=3)
lines(out[,c(1,2)], col=4)
