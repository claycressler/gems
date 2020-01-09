source("logistic_GEM.R")
library(tidyverse)
library(magrittr)
library(deSolve)
library(parallel)

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

qg_model_2 <- function(t, y, pars) {
    bs <- pars["bs"]
    ds <- pars["ds"]
    slope <- pars["slope"]
    V <- pars["V"]

    N <- y[1]
    b <- y[2]
    d <- slope*b^2
    dNdt <- (b-bs*N)*N-(d+ds*N)*N
    dbdt <- V*(1/(d+ds*N)*(1-(2*slope*b*(b-bs*N))/(d+ds*N)))
    list(c(dNdt,dbdt))
}


## Recreating Fig. 2 of the manuscript
dt <- 0.1
N0 <- 5
bs_seq <- ds_seq <- c(0.1, 0.01)
slope <- 1/10
traitmean <- 2 ## initial mean of trait distribution
tmax <- 100
traitcv <- 0.3
h2 <- 0

## Run this on a multicore Mac
run <- FALSE
if (run) {
    set.seed(123431234)
    seeds <- floor(runif(50, 1, 1e5))
    out <- vector(mode='list', length=length(bs_seq))
    for (i in 1:length(bs_seq)) {
        print(i)
        bs <- bs_seq[i]
        ds <- ds_seq[i]
        mclapply(seeds,
                 function(s) logistic_GEM(s, dt, tmax, N0, traitmean, traitcv, h2, bs, ds, slope),
                 mc.cores=12) -> out[[i]]
    }
    saveRDS(out, file="simulation_results_h2=0_6-7-19.RDS")
}
out <- readRDS("simulation_results_h2=0_6-7-19.RDS")


## Process this data, calculating the population size and mean traits for each replicate
lapply(1:length(out),
         function(i) {
             o <- out[[i]]
             lapply(1:length(o),
                    function(j) {
                        o1 <- o[[j]]
                        print(j)
                        data.frame(time=o1$time,
                                   N=lapply(o1$traits, length) %>% unlist,
                                   med_b=lapply(o1$traits, function(t) mean(t, na.rm=T)) %>% unlist,
                                   var_b=lapply(o1$traits, function(t) var(t, na.rm=T)) %>% unlist,
                                   rep=j,
                                   bs=bs_seq[i]
                                   )
                    }
                    ) %>% do.call(rbind.data.frame,.)
         }
         ) %>% do.call(rbind.data.frame,.) -> out2


## calculate means
lapply(bs_seq,
       function(n)
           sapply(seq(0, 99.9, 0.1),
                  function(t)
                      with(subset(out2, bs==n & time>=t & time<(t+0.1)),
                           c(time=mean(time),
                             N=mean(c(N,rep(0,max(length(N),50)-length(N)))),
                             Nmin=sort(c(N, rep(0,max(length(N),50)-length(N))))[12],
                             Nmax=sort(c(N, rep(0,max(length(N),50)-length(N))))[37],
                             meanb=mean(med_b,na.rm=TRUE),
                             meanbmin=sort(med_b)[floor(0.25*length(med_b))],
                             meanbmax=sort(med_b)[floor(0.75*length(med_b))],
                             varb=mean(var_b,na.rm=TRUE),
                             varbmin=sort(var_b)[floor(0.25*length(var_b))],
                             varbmax=sort(var_b)[floor(0.75*length(var_b))],
                             bs=n
                             )
                             )
                  ) %>% t
       ) %>% do.call(rbind.data.frame,.) -> out4

lapply(bs_seq,
       function(n)
           lsoda(y=c(N=N0, b=traitmean), times=seq(0,tmax,dt), func=qg_model, parms=c(bs=n, ds=n, slope=slope, V=(traitcv*traitmean)^2*h2)) %>% as.data.frame %>% mutate(., bs=n)
       ) %>% do.call(rbind.data.frame,.) -> deter

lapply(bs_seq,
       function(n)
           lsoda(y=c(N=N0, b=traitmean), times=seq(0,tmax,dt), func=qg_model_2, parms=c(bs=n, ds=n, slope=slope, V=(traitcv*traitmean)^2*h2)) %>% as.data.frame %>% mutate(., bs=n)
       ) %>% do.call(rbind.data.frame,.) -> deter2





## Plot using out4 - should be more similar to John's figure
jpeg(file="Transient_attractor_Fig2_version_1b.jpg", height=6, width=6, units='in', res=500)
par(mfrow=c(3,length(bs_seq)), mar=c(4,4,0.5,0.5), oma=c(0.5,0.5,2,0.5))
for (n in bs_seq) {
    plot.new()
    plot.window(xlim=c(0,100), ylim=range(c(unlist(subset(out4, bs==n)[,c("N","Nmin","Nmax")]),
                                   subset(deter, bs==n)$N,
                                   subset(deter2, bs==n)$N)))
    axis(1)
    axis(2)
    box('plot')
    if (n==bs_seq[1]) mtext(side=2, line=3, "Population size")
    with(subset(out4, bs==n), polygon(c(time, rev(time)), c(Nmax, rev(Nmin)), col=gray(0.8), border=NA) )
    with(subset(out4, bs==n), lines(time, N, lwd=2))
    with(subset(deter, bs==n), lines(time, N, lwd=2, lty=2, col='red'))
    with(subset(deter2, bs==n), lines(time, N, lwd=2, lty=2, col='blue'))
    if(n==0.1) mtext(side=3, line=0.5, expression(b[s]=="0.1"), bty='n')
    if(n==0.01) mtext(side=3, line=0.5, expression(b[s]=="0.01"), bty='n')
    if(n==0.001) mtext(side=3, line=0.5, expression(b[s]=="0.001"), bty='n')
}

for (n in bs_seq) {
    plot.new()
    plot.window(xlim=c(0,100), ylim=range(c(unlist(subset(out4, bs==n)[,c("meanb","meanbmin","meanbmax")]),
                                   subset(deter, bs==n)$b,
                                   subset(deter2, bs==n)$b)))
    axis(1)
    axis(2)
    box('plot')
    if (n==bs_seq[1]) mtext(side=2, line=3, "Trait mean")
    with(subset(out4, bs==n), polygon(c(time, rev(time)), c(meanbmax, rev(meanbmin)), col=gray(0.8), border=NA) )
    with(subset(out4, bs==n), lines(time, meanb, lwd=2))
    with(subset(deter, bs==n), lines(time, b, lwd=2, lty=2, col='red'))
    with(subset(deter2, bs==n), lines(time, b, lwd=2, lty=2, col='blue'))
}

for (n in bs_seq) {
    plot.new()
    plot.window(xlim=c(0,100), ylim=range(unlist(subset(out4, bs==n)[,c("varb","varbmin","varbmax")])))
    axis(1)
    axis(2)
    box('plot')
    if (n==bs_seq[1]) mtext(side=2, line=3, "Trait variance")
    with(subset(out4, bs==n), polygon(c(time, rev(time)), c(varbmax, rev(varbmin)), col=gray(0.8), border=NA) )
    with(subset(out4, bs==n), lines(time, varb, lwd=2))
    abline(h=(traitcv*traitmean)^2*h2, lwd=2, col='red')
    abline(h=(traitcv*traitmean)^2*h2, lwd=2, lty=2, col='blue')
}
dev.off()




##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################



## high heritability version, just to see what happen
h2 <- 0.99
run <- FALSE
if (run) {
    set.seed(123431234)
    seeds <- floor(runif(50, 1, 1e5))
    out <- vector(mode='list', length=3)
    for (i in 1:3) {
        print(i)
        bs <- bs_seq[i]
        ds <- ds_seq[i]
        mclapply(seeds,
                 function(s) logistic_GEM(s, dt, tmax, N0, traitmean, traitcv, h2, bs, ds, slope),
                 mc.cores=12) -> out[[i]]
    }
    saveRDS(out, file="~/simulation_results_2_6-7-19.RDS")
}
out <- readRDS("simulation_results_2_6-7-19.RDS")

## Process this data, calculating the population size and mean traits for each replicate
lapply(1:length(out),
         function(i) {
             o <- out[[i]]
             lapply(1:length(o),
                    function(j) {
                        o1 <- o[[j]]
                        data.frame(time=o1$time,
                                   N=lapply(o1$traits, length) %>% unlist,
                                   med_b=lapply(o1$traits, function(t) mean(t, na.rm=T)) %>% unlist,
                                   var_b=lapply(o1$traits, function(t) var(t, na.rm=T)) %>% unlist,
                                   rep=j,
                                   bs=bs_seq[i]
                                   )
                    }
                    ) %>% do.call(rbind.data.frame,.)
         }
         ) %>% do.call(rbind.data.frame,.) -> out2


## calculate means
lapply(bs_seq,
       function(n)
           sapply(seq(0, 99.9, 0.1),
                  function(t)
                      with(subset(out2, bs==n & time>=t & time<(t+0.1)),
                           c(time=mean(time),
                             N=mean(c(N,rep(0,max(length(N),50)-length(N)))),
                             Nmin=sort(c(N, rep(0,max(length(N),50)-length(N))))[12],
                             Nmax=sort(c(N, rep(0,max(length(N),50)-length(N))))[37],
                             meanb=mean(med_b,na.rm=TRUE),
                             meanbmin=sort(med_b)[floor(0.25*length(med_b))],
                             meanbmax=sort(med_b)[floor(0.75*length(med_b))],
                             varb=mean(var_b,na.rm=TRUE),
                             varbmin=sort(var_b)[floor(0.25*length(var_b))],
                             varbmax=sort(var_b)[floor(0.75*length(var_b))],
                             bs=n
                             )
                             )
                  ) %>% t
       ) %>% do.call(rbind.data.frame,.) -> out4

lapply(bs_seq,
       function(n)
           lsoda(y=c(N=N0, b=traitmean), times=seq(0,tmax,dt), func=qg_model, parms=c(bs=n, ds=n, slope=slope, V=(traitcv*traitmean)^2*h2)) %>% as.data.frame %>% mutate(., bs=n)
       ) %>% do.call(rbind.data.frame,.) -> deter

lapply(bs_seq,
       function(n)
           lsoda(y=c(N=N0, b=traitmean), times=seq(0,tmax,dt), func=qg_model_2, parms=c(bs=n, ds=n, slope=slope, V=(traitcv*traitmean)^2*h2)) %>% as.data.frame %>% mutate(., bs=n)
       ) %>% do.call(rbind.data.frame,.) -> deter2





## Plot using out4 - should be more similar to John's figure
jpeg(file="Transient_attractor_Fig2_version_1b.jpg", height=6, width=6, units='in', res=500)
par(mfrow=c(3,3), mar=c(4,4,0.5,0.5), oma=c(0.5,0.5,2,0.5))
for (n in bs_seq) {
    plot.new()
    plot.window(xlim=c(0,100), ylim=range(c(unlist(subset(out4, bs==n)[,c("N","Nmin","Nmax")]),
                                   subset(deter, bs==n)$N,
                                   subset(deter2, bs==n)$N)))
    axis(1)
    axis(2)
    box('plot')
    if (n==bs_seq[1]) mtext(side=2, line=3, "Population size")
    with(subset(out4, bs==n), polygon(c(time, rev(time)), c(Nmax, rev(Nmin)), col=gray(0.8), border=NA) )
    with(subset(out4, bs==n), lines(time, N, lwd=2))
    with(subset(deter, bs==n), lines(time, N, lwd=2, lty=2, col='red'))
    with(subset(deter2, bs==n), lines(time, N, lwd=2, lty=2, col='blue'))
    if(n==0.1) mtext(side=3, line=0.5, expression(b[s]=="0.1"), bty='n')
    if(n==0.01) mtext(side=3, line=0.5, expression(b[s]=="0.01"), bty='n')
    if(n==0.001) mtext(side=3, line=0.5, expression(b[s]=="0.001"), bty='n')
}

for (n in bs_seq) {
    plot.new()
    plot.window(xlim=c(0,100), ylim=range(c(unlist(subset(out4, bs==n)[,c("meanb","meanbmin","meanbmax")]),
                                   subset(deter, bs==n)$b,
                                   subset(deter2, bs==n)$b)))
    axis(1)
    axis(2)
    box('plot')
    if (n==bs_seq[1]) mtext(side=2, line=3, "Trait mean")
    with(subset(out4, bs==n), polygon(c(time, rev(time)), c(meanbmax, rev(meanbmin)), col=gray(0.8), border=NA) )
    with(subset(out4, bs==n), lines(time, meanb, lwd=2))
    with(subset(deter, bs==n), lines(time, b, lwd=2, lty=2, col='red'))
    with(subset(deter2, bs==n), lines(time, b, lwd=2, lty=2, col='blue'))
}

for (n in bs_seq) {
    plot.new()
    plot.window(xlim=c(0,100), ylim=range(unlist(subset(out4, bs==n)[,c("varb","varbmin","varbmax")])))
    axis(1)
    axis(2)
    box('plot')
    if (n==bs_seq[1]) mtext(side=2, line=3, "Trait variance")
    with(subset(out4, bs==n), polygon(c(time, rev(time)), c(varbmax, rev(varbmin)), col=gray(0.8), border=NA) )
    with(subset(out4, bs==n), lines(time, varb, lwd=2))
    abline(h=(traitcv*traitmean)^2*h2, lwd=2, col='red')
    abline(h=(traitcv*traitmean)^2*h2, lwd=2, lty=2, col='blue')
}
dev.off()


##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################



## An alternative version with slightly larger starting populations and slightly different carrying capacities
dt <- 0.1
N0 <- 10
bs_seq <- ds_seq <- c(0.05, 0.01, 0.005)
slope <- 1/10
traitmean <- 2 ## initial mean of trait distribution
tmax <- 100
traitcv <- 0.3
h2 <- 0.75

run <- FALSE
if (run) {
    set.seed(123431234)
    seeds <- floor(runif(50, 1, 1e5))
    out <- vector(mode='list', length=3)
    for (i in 1:3) {
        print(i)
        bs <- bs_seq[i]
        ds <- ds_seq[i]
        mclapply(seeds,
                 function(s) logistic_GEM(s, dt, tmax, N0, traitmean, traitcv, h2, bs, ds, slope),
                 mc.cores=12) -> out[[i]]
    }
    saveRDS(out, file="~/simulation_results_2_6-6-19.RDS")
}
out <- readRDS("simulation_results_2_6-6-19.RDS")

lapply(1:length(out),
         function(i) {
             o <- out[[i]]
             lapply(1:length(o),
                    function(j) {
                        o1 <- o[[j]]
                        data.frame(time=o1$time,
                                   N=lapply(o1$traits, length) %>% unlist,
                                   med_b=lapply(o1$traits, function(t) mean(t, na.rm=T)) %>% unlist,
                                   var_b=lapply(o1$traits, function(t) var(t, na.rm=T)) %>% unlist,
                                   rep=j,
                                   bs=bs_seq[i]
                                   )
                    }
                    ) %>% do.call(rbind.data.frame,.)
         }
         ) %>% do.call(rbind.data.frame,.) -> out2

lapply(bs_seq,
       function(n)
           sapply(seq(0, 99.9, 0.1),
                  function(t)
                      subset(out2, bs==n & time>=t & time<(t+0.1)) %>% apply(., 2, function(x) mean(x, na.rm=T))
                  ) %>% t
       ) %>% do.call(rbind.data.frame,.) -> out3

lapply(bs_seq,
       function(n)
           lsoda(y=c(N=N0, b=traitmean), times=seq(0,tmax,dt), func=qg_model, parms=c(bs=n, ds=n, slope=slope, V=(traitcv*traitmean)^2*h2)) %>% as.data.frame %>% mutate(., bs=n)
       ) %>% do.call(rbind.data.frame,.) -> deter

lapply(bs_seq,
       function(n)
           lsoda(y=c(N=N0, b=traitmean), times=seq(0,tmax,dt), func=qg_model_2, parms=c(bs=n, ds=n, slope=slope, V=(traitcv*traitmean)^2*h2)) %>% as.data.frame %>% mutate(., bs=n)
       ) %>% do.call(rbind.data.frame,.) -> deter2


## Plot
jpeg(file="Transient_attractor_Fig2_version_2.jpg", height=6, width=6, units='in', res=500)
par(mfrow=c(3,3), mar=c(4,4,0.5,0.5), oma=rep(0.5,4))
for (n in bs_seq) {
    plot.new()
    plot.window(xlim=c(0,100), ylim=range(subset(out2, bs==n)$N))
    axis(1)
    axis(2)
    box('plot')
    if (n==bs_seq[1]) mtext(side=2, line=3, "Population size")
    for (i in 1:50)
        with(subset(out2, rep==i & bs==n), lines(time, N, col=gray(0.8)))
    with(subset(out3, bs==n), lines(time, N, lwd=2))
    with(subset(deter, bs==n), lines(time, N, lwd=2, lty=2, col='red'))
    with(subset(deter2, bs==n), lines(time, N, lwd=2, lty=2, col='blue'))
}

for (n in bs_seq) {
    plot.new()
    plot.window(xlim=c(0,100), ylim=range(subset(out2, bs==n)$med_b, na.rm=T))
    axis(1)
    axis(2)
    box('plot')
    if (n==bs_seq[1]) mtext(side=2, line=3, "Trait mean")
    for (i in 1:50)
        with(subset(out2, rep==i & bs==n), lines(time, med_b, col=gray(0.8)))
    with(subset(out3, bs==n), lines(time, med_b, lwd=2))
    with(subset(deter, bs==n), lines(time, b, lwd=2, lty=2, col='red'))
    with(subset(deter2, bs==n), lines(time, b, lwd=2, lty=2, col='blue'))
}

for (n in bs_seq) {
    plot.new()
    plot.window(xlim=c(0,100), ylim=c(0,1))
    axis(1)
    axis(2)
    box('plot')
    if (n==bs_seq[1]) mtext(side=2, line=3, "Trait variance")
    for (i in 1:50)
        with(subset(out2, rep==i & bs==n), lines(time, var_b, col=gray(0.8)))
    with(subset(out3, bs==n), lines(time, var_b, lwd=2))
    abline(h=(traitcv*traitmean)^2*h2, lwd=2, col='red')
    abline(h=(traitcv*traitmean)^2*h2, lwd=2, lty=2, col='blue')
}

dev.off()



##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################




## Recreating Fig. 2 of the manuscript
dt <- 0.1
N0 <- 5
bs_seq <- ds_seq <- c(0.1, 0.01, 0.001)
slope <- 1/10
traitmean <- 2 ## initial mean of trait distribution
tmax <- 100
traitcv <- 0.3
h2 <- 0.75

## Run this on a multicore Mac
run <- FALSE
if (run) {
    set.seed(123431234)
    seeds <- floor(runif(50, 1, 1e5))
    out <- vector(mode='list', length=3)
    for (i in 1:3) {
        print(i)
        bs <- bs_seq[i]
        ds <- ds_seq[i]
        mclapply(seeds,
                 function(s) logistic_GEM(s, dt, tmax, N0, traitmean, traitcv, h2, bs, ds, slope),
                 mc.cores=12) -> out[[i]]
    }
    saveRDS(out, file="~/simulation_results_6-7-19.RDS")
}

tmax <- 60
set.seed(123431234)
seeds <- floor(runif(20, 1, 1e5))
i <- 2
print(i)
bs <- bs_seq[i]
ds <- ds_seq[i]

logistic_GEM_storage(seeds[1], tmax, N0, traitmean, traitcv, h2, bs, ds, slope) -> out

## parent-offspring regression
with(out, plot(ancTrait, trait))
with(out, lm(trait~ancTrait)) %>% summary

mclapply(seeds,
         function(s) logistic_GEM_storage(s, tmax, N0, traitmean, traitcv, h2, bs, ds, slope),
         mc.cores=4) -> out


out <- readRDS("simulation_results_6-7-19.RDS")
