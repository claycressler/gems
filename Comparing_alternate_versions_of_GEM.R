library(magrittr)
library(deSolve)
library(parallel)
library(pryr)

source("logistic_GEM.R")

dt <- 0.1 ## size of time steps
N0 <- 20 ## initial population size
traitmean <- 1.8 ## initial mean of trait distribution
traitcv <- 0.3 ## initial coefficient of variation of trait distribution
h2 <- 0.75 ## heritability
## logistic equation parameters
bs <- 0.005
ds <- 0.005
slope <- 1/6
tmax <- 100

## start with a straight comparison
sims <- vector(mode='list', length=3)
set.seed(1234897)
nreps <- 40
seeds <- floor(runif(nreps, 1, 100000))
mclapply(seeds,
         function(seed)
             logistic_GEM(seed, dt, tmax, N0, traitmean, traitcv, h2, bs, ds, slope),
         mc.cores=min(nreps,12)) -> sims[[1]]
mclapply(seeds,
         function(seed)
             logistic_GEM2(seed, dt, tmax, N0, traitmean, traitcv, h2, bs, ds, slope),
         mc.cores=min(nreps,12)) -> sims[[2]]
mclapply(seeds,
         function(seed)
             logistic_GEM3(seed, dt, tmax, N0, traitmean, traitcv, h2, bs, ds, slope),
         mc.cores=min(nreps,12)) -> sims[[3]]

saveRDS(sims, file="~/methodcomp.RDS")

sims <- readRDS("methodcomp.RDS")

lapply(1:length(sims),
       function(i)
           lapply(1:length(sims[[i]]),
                  function(j)
                      data.frame(time=sims[[i]][[j]]$time,
                                 N=unlist(lapply(sims[[i]][[j]]$traits, length)),
                                 med_b=lapply(sims[[i]][[j]]$traits, median) %>% unlist,
                                 var_b=lapply(sims[[i]][[j]]$traits, var) %>% unlist,
                                 rep=j,
                                 method=i
                                 )
                  ) %>% do.call(rbind.data.frame,.)
       ) %>% do.call(rbind.data.frame,.) -> out2

lapply(1:3,
       function(n)
           sapply(seq(0, 99.9, 0.1),
                  function(t)
                      subset(out2, method==n & time>=t & time<(t+0.1)) %>% apply(., 2, mean)
                  ) %>% t
       ) %>% do.call(rbind.data.frame,.) -> out3




png(file="Method_comparsion.png", width=10, height=6, units='in', res=400)
par(mfrow=c(3,3), mar=c(1,1,0.5,0.5), oma=c(4,4,0,0))
for (n in 1:3) {
    plot.new()
    plot.window(xlim=c(0,100), ylim=range(out2$N))
    axis(1, tick=T, labels=F)
    if (n==1) {
        axis(2)
        mtext(side=2, line=3, "Population size")
    }
    else axis(2, tick=T, labels=F)
    box(bty='L')
    for (i in 1:20) with(subset(out2, method==n & rep==i), lines(time, N, col=gray(0.7)))
    with(subset(out3, method==n), lines(time, N, lwd=2))
}

for (n in 1:3) {
    plot.new()
    plot.window(xlim=c(0,100), ylim=range(out2$med_b))
    axis(1, tick=T, labels=F)
    if (n==1) {
        axis(2)
        mtext(side=2, line=3, "Median trait")
    }
    else axis(2, tick=T, labels=F)
    box(bty='L')
    for (i in 1:20) with(subset(out2, method==n & rep==i), lines(time, med_b, col=gray(0.7)))
    with(subset(out3, method==n), lines(time, med_b, lwd=2))
}

for (n in 1:3) {
    plot.new()
    plot.window(xlim=c(0,100), ylim=range(out2$var_b))
    axis(1)
    if (n==1) {
        axis(2)
        mtext(side=2, line=3, "Trait variance")
    }
    else axis(2, tick=T, labels=F)
    box(bty='L')
    for (i in 1:20) with(subset(out2, method==n & rep==i), lines(time, var_b, col=gray(0.7)))
    with(subset(out3, method==n), lines(time, var_b, lwd=2))
}
dev.off()

## Why is evolution so much slower in the alternative versions of the GEM?

