library(magrittr)
library(parallel)
library(deSolve)

source("GEM_sir_model.R")

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

## Compare the SIR GEM to the QG model
qg <- ode(y = c(S=100, I=50, v=0.3), times=seq(0, 200, 0.1), func=qg_sir_model, parms=c(b = 2, m = 1.7, B0 = 0.05, h = 0.1, V=0.001))

sirout <- sir_GEM(1287349876, tmax=200, dt=0.1, initialstate=c(S=60,I=20), traitmean=0.2, traitcv=0.2, h2=0.85, params=c(b = 2, m = 1.7, B0 = 0.05, h = 0.1))

## Not bad! There is a huge peak in the abundance of infected individuals, but that is really to be expected, given the large number of susceptible hosts.
par(mfrow=c(2,2))
lapply(sirout$states, function(x) x[1]) %>% unlist %>% plot(., type='l', ylab="Susc abund")
lines(qg[,2], col=2)
lapply(sirout$states, function(x) x[2]) %>% unlist %>% plot(., type='l', ylab="Inf abund")
lines(qg[,3], col=2)
lapply(sirout$traits, median) %>% unlist %>% plot(., type='l', ylab="Median trait")
lines(qg[,4], col=2)
lapply(sirout$traits, var) %>% unlist %>% plot(., type='l', ylab="Var trait")

## Realistically, the best way to handle this would be with 
## Analytical calculation of the equilibria of the SI model
si_equil <- function(pars) {
    b <- pars["b"]
    m <- pars["m"]
    B0 <- pars["B0"]
    h <- pars["h"]
    V <- pars["V"]
    c("S"=(sqrt(h)+sqrt(m))^2/B0,
      "I"=((sqrt(h)+sqrt(m))^2*(b-m))/(B0*(-b+sqrt(h*m)+m)),
      "v"=sqrt(h*m))
}

## That's not bad! 

sirout <- sir_GEM(1287349876, tmax=200, dt=0.1, initialstate=c(S=60,I=20), traitmean=0.6, traitcv=0.2, h2=0.85, params=c(b = 2, m = 1.7, B0 = 0.08, h = 0.1), report=TRUE)

out0 <- ode(y = c(S=60, I=20, v=0.6), times=seq(0, 200, 0.1), func=qg_sir_model, parms=c(b = 2, m = 1.7, B0 = 0.12, h = 0.1, V=0.001))
plot(out0[,c(1,3)], type='l')

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





