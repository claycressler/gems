qg_model <- function(t, y, pars) {
    C <- y[1]
    P <- y[2]
    b <- y[3]

    K <- pars["K"]
    e <- pars["e"]
    m <- pars["m"]
    s <- pars["s"]
    amin <- pars["amin"]
    V <- pars["V"]
    a <- amin + s*b^2

    dCdt <- b*C*(1-C/K) - a*C*P
    dPdt <- e*a*C*P - m*P
    dbdt <- V*(1-C/K - 2*s*b*P)
    list(c(dCdt,dPdt,dbdt))
}

pars <- c(K=100, e=0.1, m=0.1, s=0.1, amin=0.01, V=0.01)
y0 <- c(C=99, P=1, b=0.1)
times <- seq(0,1000,0.1)

library(deSolve)
out <- ode(y=y0, times=times, parms=pars, func=qg_model)

par(mfrow=c(1,3), mar=c(5,5,0.5,0.5), oma=rep(0,4))
with(as.data.frame(out), plot(time, C, type='l', lwd=2))
## how far away from the ecological equilibrium, given the current trait value, are the system dynamics?
with(as.data.frame(out), lines(time, pars["m"]/(pars["e"]*(pars["amin"]+pars["s"]*b^2)), type='l', col=2, lwd=2, lty=2))
legend(x='topright', c("Transient ecological dynamics", "Ecological equilibrium given current trait"), bty='n', lty=c(1,2), col=c(1,2))

with(as.data.frame(out), plot(time, P, type='l', lwd=2))
with(as.data.frame(out), lines(time, b*(pars["K"]-pars["m"]/(pars["e"]*(pars["amin"]+pars["s"]*b^2)))/(pars["K"]*(pars["amin"]+pars["s"]*b^2)), type='l', col=2, lty=2, lwd=2))

with(as.data.frame(out), plot(time, b, type='l', ylim=c(0,0.5), lwd=2))
## evolutionary equilibrium given C and K
with(as.data.frame(out), lines(time, (pars["K"]-C)/(2*pars["s"]*pars["K"]*P), col="blue", lty=2, lwd=2))
legend(x='topright', c("Transient trait dynamics", "Transient evolutionary equilibrium"), bty='n', lty=c(1,2), col=c(1,"blue"))



## what if you "trap" the system in a transient, for example by harvesting the prey population.
qg_model_2 <- function(t, y, pars) {
    C <- y[1]
    P <- y[2]
    b <- y[3]

    K <- pars["K"]
    e <- pars["e"]
    m <- pars["m"]
    s <- pars["s"]
    amin <- pars["amin"]
    V <- pars["V"]
    a <- amin + s*b^2

    dCdt <- b*C*(1-C/K) - a*C*P
    dPdt <- e*a*C*P - m*P
    dbdt <- V*(1-C/K - 2*s*b*P)
    list(c(0,dPdt,dbdt))
}

## start the system with P and b and the eco-evolutionary equilibrium, but drop C to its harvested value
y0 <- c(C=30, P=7.9, b=0.316)
out2 <- ode(y=y0, times=times, parms=pars, func=qg_model_2)

par(mfrow=c(1,2), mar=c(5,5,0.5,0.5), oma=rep(0,4))
with(as.data.frame(out2), plot(time, P, type='l', lwd=2))

with(as.data.frame(out2), plot(time, (pars["K"]-C)/(2*pars["s"]*pars["K"]*P), type='l', col="blue", lty=2, lwd=2, ylab="b"))
with(as.data.frame(out2), lines(time, b, type='l', lwd=2))
## evolutionary equilibrium given C and K
legend(x='topright', c("Transient trait dynamics", "Transient evolutionary equilibrium"), bty='n', lty=c(1,2), col=c(1,"blue"))


