library(deSolve)
qg_model <- function(t, y, pars) {
  bs <- pars["bs"]
  ds <- pars["ds"]
  slope <- pars["slope"]
  V <- pars["V"]
  
  N <- y[1]
  b <- y[2]
  d <- slope*b^2
  dNdt <- (b-bs*N)*N-(d+ds*N)*N
  dbdt <- V*(1-2*b*slope)
  list(c(dNdt,dbdt))
}

parms=c(bs=0.04, ds=0.04, slope=0.3/1.8^2, V=0.75*(0.3*1.8)^2)
times <- seq(0,200,0.001)
y0 <- c(N=5, b=1.8)

out <- ode(y=y0, times=times, func=qg_model, parms=parms)

## given any possible age at birth, compute the expected lifetime of that individual, given the dynamics of N

## LRS accumulate, and can be expressed as the integral of the reproductive rate and the probability of survival to any given age

## so for an individual born at time 0 with trait b=1.8, d=0.3
b = 1.8
d = unname(parms["slope"])*b^2
## probability of survival to any given age, for an individual born at time 0
psurv <- exp(-cumsum((d+unname(parms["ds"])*out[1:nrow(out),"N"])*0.001))
## birth rate of that individual
brate <- (b-unname(parms["bs"])*out[,"N"])
## lifetime reproductive success
LRS <- sum(psurv*brate*0.001)

## Apply this idea across the whole N trajectory, comparing the LRS that would be achieved by individuals with different trait values against this background of population size change
dt <- 0.001
out <- ode(y=y0, times=seq(0,200,dt), func=qg_model, parms=parms)
bseq <- seq(1,5,0.01)
lrs <- sapply(bseq, 
              function(b) 
                sum(
                  (
                    (b-unname(parms["bs"])*out[,"N"])*
                      exp(-cumsum((unname(parms["slope"])*b^2+unname(parms["ds"])*out[1:nrow(out),"N"])*dt))
                  )*dt
                )
)
bseq[which.max(lrs)]
## Okay just to check: go much later when dynamics are equilibrial, e.g., starting at time t=100
b = unname(out[100/dt+1,"b"])
d = unname(parms["slope"])*b^2
## it should also be the case that any b that is different from this should have a lower fitness (since this is the ESS)
bseq <- seq(b-2,b+2,0.01)
lrs <- sapply(bseq, 
              function(b) 
                sum(
                  (
                    (b-unname(parms["bs"])*out[(100/dt+1):nrow(out),"N"]) *
                      exp(-cumsum((unname(parms["slope"])*b^2+unname(parms["ds"])*out[(100/dt+1):nrow(out),"N"])*dt))
                  )*dt
                )
)
bseq[which.max(lrs)]
b

## Now what I want to do is calculate the trajectory of the strategy that maximizes the LRS for the first 100 time steps
LRSmax <- rep(0, length(seq(1,140001,100)))
for (j in 1:length(LRSmax)) {
  i <- seq(1,140001,100)[j]
  bseq <- seq(max(c(0.01,out[i,"b"]-2)),out[i,"b"]+2,0.01)
  lrs <- sapply(bseq, 
                function(b) 
                  sum(
                    (
                      (b-unname(parms["bs"])*out[i:nrow(out),"N"]) *
                        exp(-cumsum((unname(parms["slope"])*b^2+unname(parms["ds"])*out[i:nrow(out),"N"])*0.001))
                    )*0.001
                  )
  )
  LRSmax[j] <- bseq[which.max(lrs)]
  print(i)
}
## By t = 140, the trait that maximizes the LRS is also the same as the ESS
par(mfrow=c(1,1), mar=c(5,5,1,1), oma=c(0,0,0,0))
plot(out[,c(1,3)], type='l', lwd=2)
lines(seq(0,140,0.001)[seq(1,140001,100)], LRSmax, lwd=2, lty=2, col=2)
abline(h=1/(2*parms['slope']), col=4)
  
  
  

## Okay, now do it for realsies
## the function for simulating the GEM with storage
source("logistic_GEM_with_storage.R")
library(parallel)

## Here is the set of GEM parameters 
traitmean <- 1.8
traitcv <- 0.3
h2 <- 0.75
slope <- 0.3/1.8^2 ## from dmin = slope*bmax^2 => 0.3 = slope*1.8^2
tmax <- 400
N0 <- 5
bs <- ds <- 0.05
## RNG seeds
set.seed(101)
seeds <- floor(runif(50, 1, 1e7))
mclapply(seeds,
         function(s) logistic_GEM_storage(seed=s, tmax=tmax, N0=N0, traitmean=traitmean, traitcv=traitcv, h2=h2, bs=bs, ds=ds, slope=slope),
         mc.cores=15
) -> out

## Two ways that you could compute stuff: one is to look at the mean trajectory across the stochastic simulations, and ask how well the mean phenotype across the simulations is tracking the predicted ESS across the simulations. The other is to look at this tracking within each simulation, and then average. We'll start with the former.
makeTrajectories <- function(output) {
  time.seq <- seq(0,400,0.1)
  sapply(1:length(time.seq),
         function(t) {
           ## who is alive?
           if (t < length(time.seq))
             alive <- subset(output, tBirth <= time.seq[t] & tDeath > time.seq[t])
           else alive <- subset(output, is.na(tDeath))
           ## at the current moment, how big is the population size and what are the mean and SD of the current trait distribution?
           c(nrow(alive), mean(alive$trait), sd(alive$trait))
         }) %>% t %>% as.data.frame %>% mutate(., time=time.seq) -> traj
  return(traj)
}
mclapply(out, 
         function(o) makeTrajectories(o),
         mc.cores=15) -> trajs
for (i in 1:length(trajs)) 
  mutate(trajs[[i]], rep=i) -> trajs[[i]]

## Compute the LRS maximizing strategy for each trajectory individually (up to timestep 300)
LRSmaximization <- function(N) {
  slope <- 0.3/1.8^2
  bs <- ds <- 0.05
  LRSmax <- rep(0, length(seq(0,300,0.1))) ## set up storage
  for (j in 1:length(LRSmax)) {
    t <- seq(0,300,0.1)[j] ## choose a timepoint
    bseq <- seq(2,6,0.01) ## set up a series of b strategies to test
    lrs <- sapply(bseq, 
                  function(b) 
                    sum(((b-bs*N[j:length(N)]) * exp(-cumsum((slope*b^2+ds*N[j:length(N)])*0.01)))*0.01)) ## calculate LRS for each b strategy
    LRSmax[j] <- bseq[which.max(lrs)] ## find the maximum 
  }
  return(LRSmax)
}
mclapply(trajs, 
         function(tr) LRSmaximization(tr$V1),
         mc.cores=15) -> LRSmaxes
for (i in 1:length(trajs)) 
  mutate(trajs[[i]], LRSmax=c(LRSmaxes[[i]],rep(NA,1000))) -> trajs[[i]]

plot.new()
plot.window(xlim=c(0,3001), ylim=c(1.4,5.5))
axis(1); axis(2); box('plot')
for (i in 1:50) {
  lines(trajs[[i]]$V2, col=gray(0.5,0.5))
  lines(trajs[[i]]$LRSmax, col='pink')
}

## Find the summary trajectories and the LRS maximization strategies for the median N trajectory
trajs %>% do.call(rbind.data.frame, .) -> trajs2
colnames(trajs2)[1:3] <- c("N","mean","sd")
## compute the median trajectories
trajs2 %>% group_by(time) %>% summarise(med.N=median(N,na.rm=T),med.b=median(mean,na.rm=T),med.LRSmax=median(LRSmax,na.rm=T)) -> summary.trajs

## Can also compute the LRS maximizing strategy at the median N
med.LRSmax <- LRSmaximization(summary.trajs$med.N)

## Compare the two maximization strategies
plot(summary.trajs$med.LRSmax, type='l', lwd=2)
lines(med.LRSmax, lwd=2, lty=2, col=2)
lines(summary.trajs$med.b, lwd=2, col=4)


## Plot the deterministic QG expectation and the observed stochastic trajectory of the mean trait
qg_model <- function(t, y, pars) {
  bs <- pars["bs"]
  ds <- pars["ds"]
  slope <- pars["slope"]
  V <- pars["V"]
  
  N <- y[1]
  b <- y[2]
  d <- slope*b^2
  dNdt <- (b-bs*N)*N-(d+ds*N)*N
  dbdt <- V*(1-2*b*slope)
  list(c(dNdt,dbdt))
}
parms=c(bs=0.05, ds=0.05, slope=0.3/1.8^2, V=0.75*(0.3*1.8)^2)
times <- seq(0,400,0.1)
y0 <- c(N=5, b=1.8)
qgsim <- ode(y=y0, times=times, func=qg_model, parms=parms)

## Finally, pull out the ESS and the "mean LRS maximizing strategy"
ESS <- 1/(2*slope)
## The mean LRSmax over the last 100 timesteps
meanLRSmax1 <- mean(tail(LRSmax,1000))
## The b that maximizes LRS at the mean N from t=200 to t=300
meanN <- filter(summary.trajs, time >=300, time <=400)$med.N %>% mean
meanLRSmax2 <- (bs*slope*meanN + sqrt(ds*slope*meanN+(bs*slope*meanN)^2))/slope

## Plot it all together
## Deterministic expectation
par(mar=c(4,4,1,1), oma=rep(0,4))
plot(qgsim[1:3001,c(1,3)], type='l', col=4, lwd=2)
## Stochastic realization
lines(summary.trajs$time, summary.trajs$med.b, lwd=2)
## LRS maximizing strategy, given N trajectory
lines(seq(0,300,0.1), LRSmax, col=2, lwd=2)
## The mean LRSmax over the last 100 timesteps
abline(h=meanLRSmax2, lwd=1.5, lty=2, col=2)
abline(h=ESS, lwd=1.5, lty=2, col=4)