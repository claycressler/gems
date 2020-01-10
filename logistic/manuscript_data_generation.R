## This code produces versions of the figures used in the "transient evolutionary attractor" manuscript. Basically, it uses the logistic_GEM to simulate ecological and evolutionary dynamics, varying the strength of density-dependence (and thus the population size). This allows us to see how stochasticity (e.g., drift) can play an important role in both the evolutionary dynamics and the ecological dynamics.
library(magrittr)
library(parallel)

## the function for simulating the GEM with storage
source("logistic_GEM_with_storage.R")

## Here is the set of GEM parameters that will remain constant across the different simulations
traitmean <- 1.8
traitcv <- 0.3
h2 <- 0.75
slope <- 0.3/1.8^2 ## from dmin = slope*bmax^2 => 0.3 = slope*1.8^2
tmax <- 400
N0 <- 5

## RNG seeds
set.seed(101)
seeds <- floor(runif(500, 1, 1e7))

## Parameter sets varying density dependence
## dd <- c(0.05, 0.025, 0.01, 0.005)
## Here are the values were are actually planning to use in the manuscript
dd <- c(0.05, 0.025, 0.01)

for (val in dd) {
    print(val)
    ## set the density dependence
    bs <- ds <- val

    ## simulate stochastic replicates of the GEM
    mclapply(seeds,
             function(s) logistic_GEM_storage(seed=s, tmax=tmax, N0=N0, traitmean=traitmean, traitcv=traitcv, h2=h2, bs=bs, ds=ds, slope=slope),
             mc.cores=15
             ) -> out

    ## Process the data for each of these stochastic simulations into trajectories of population size and mean and sd of the trait distribution
    time.seq <- 0:tmax
    lapply(1:length(out),
           function(o)
               sapply(1:length(time.seq),
                      function(t) {
                          ## who is alive?
                          if (t < length(time.seq))
                              alive <- subset(out[[o]], tBirth <= time.seq[t] & tDeath > time.seq[t])
                          else alive <- subset(out[[o]], is.na(tDeath))
                          ## at the current moment, how big is the population size and what are the mean and SD of the current trait distribution?
                          c(nrow(alive), mean(alive$trait), sd(alive$trait))
                      }) %>% t %>% as.data.frame %>% mutate(., time=time.seq, rep=o)
           ) %>% do.call(rbind.data.frame, .) -> trajs
    colnames(trajs)[1:3] <- c("N","mean","sd")

    ## Save the output dataframes
    saveRDS(out, file=paste0("logistic_GEM_storage_output_bs=ds=",val,".RDS"))
    saveRDS(trajs, file=paste0("logistic_GEM_storage_trajectories_bs=ds=",val,".RDS"))
}

## One extra parameter set, just to match John
bs <- ds <- 0.04
## simulate stochastic replicates of the GEM
mclapply(seeds,
         function(s) logistic_GEM_storage(seed=s, tmax=tmax, N0=N0, traitmean=traitmean, traitcv=traitcv, h2=h2, bs=bs, ds=ds, slope=slope),
         mc.cores=15
         ) -> out

## Process the data for each of these stochastic simulations into trajectories of population size and mean and sd of the trait distribution
time.seq <- 0:tmax
lapply(1:length(out),
       function(o)
           sapply(1:length(time.seq),
                  function(t) {
                      ## who is alive?
                      if (t < length(time.seq))
                          alive <- subset(out[[o]], tBirth <= time.seq[t] & tDeath > time.seq[t])
                      else alive <- subset(out[[o]], is.na(tDeath))
                      ## at the current moment, how big is the population size and what are the mean and SD of the current trait distribution?
                      c(nrow(alive), mean(alive$trait), sd(alive$trait))
                  }) %>% t %>% as.data.frame %>% mutate(., time=time.seq, rep=o)
       ) %>% do.call(rbind.data.frame, .) -> trajs
colnames(trajs)[1:3] <- c("N","mean","sd")

## Save the output dataframes
saveRDS(out, file=paste0("logistic_GEM_storage_output_bs=ds=",val,".RDS"))
saveRDS(trajs, file=paste0("logistic_GEM_storage_trajectories_bs=ds=",val,".RDS"))


## Create a figure showing the mean population size trajectory and the middle 50% of observations at each time point

plotq <- FALSE
if (plotq) {
    dd <- c(0.05, 0.025, 0.01)

    library(magrittr)
    library(deSolve)
    ## Simulate the QG model for comparison
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

    ## This takes a long time to run - easier to just do it outside the figure plotting loop to speed things along
    for (val in dd) {

        ## Pull out the relationship between traits and fitness for the first 20 and last 20 timesteps
        ## Either exclude individuals that lived less than one day (there are individuals that only live for <0.1 time steps but manage to have an offspring, so their R0 values are enormous)
        ## or just plot realized fitness (total number of offspring)
        out <- readRDS(file=paste0("logistic_GEM_storage_output_bs=ds=",val,".RDS")) ## uncomment this if you want to be able to do something like plotting the relationship between traits and fitness
        lapply(out, function(o) subset(o, tBirth <= 20 & tDeath > 20)[c('trait','nOff')]) %>% do.call(rbind.data.frame, .) -> nOff.early
        ## break trait space into small intervals and determine in which interval observed fitness is the highest
        ## important consideration: want to have approximately equal numbers of individuals per bin, otherwise mean across individuals can be strongly skewed by one or two animals
        ## the approach on the next two lines divides into trait bins of each width, which doesn't totally work
        ##tseq <- seq(floor(10*min(nOff.early$trait))/10, ceiling(10*max(nOff.early$trait))/10, 0.2)
        ##sapply(1:(length(tseq)-1), function(tr) subset(nOff.early, trait > tseq[tr] & trait < tseq[tr+1])$nOff %>% mean)
        ##rseq <- seq(1,nrow(nOff.early),1000)
        ##sapply(1:(length(rseq)-1), function(r) with(nOff.early, nOff[order(trait)][rseq[r]:rseq[r+1]]) %>% mean)
        lapply(out, function(o) subset(o, tBirth <= 380 & tDeath > 380)[c('trait','nOff')]) %>% do.call(rbind.data.frame, .) -> nOff.late
        assign(paste0("nOff.early.",val), nOff.early)
        assign(paste0("nOff.late.",val), nOff.late)

        ## then remove out to free up space
        rm(out)
    }

    png(filename="Simulation_trajectories_for_bs=ds=0.05_0.025_0.01.png", width=12, height=8, units='in', res=300)
    par(mfcol=c(4,length(dd)), mar=c(3,3,1,0.5), oma=c(1,3,0,0), family="HersheySans")
    for (val in dd) {

        ## trajs just contains the trajectories themselves
        trajs <- readRDS(file=paste0("logistic_GEM_storage_trajectories_bs=ds=",val,".RDS"))

        nOff.early <- get(paste0("nOff.early.",val))
        nOff.late <- get(paste0("nOff.late.",val))
        
        ## simulate the QG trajectory
        ## note: additive genetic variance is calculated based on the definition of heritability:
        ## h2 = Va/Vp, where Va is the additive genetic variance and Vp is the total phenotypic variance. The initial phenotypic variance is the initial traitsd^2, so Vp=(traitmean*traitcv)^2. Then Va=h2*(traitmean*traitcv)^2.
        tmax <- 400
        qg <- ode(y=c(N=5, b=1.8), times=0:tmax, func=qg_model, parms=c(bs=val, ds=val, slope=0.3/1.8^2, V=0.75*(0.3*1.8)^2))
        time.seq <- 0:(tmax-1)
        sapply(time.seq,
               function(t)
                   with(subset(trajs, time==t & N > 0),
                        c(time=t,
                          med.N=median(N),
                          low50.N=sort(N)[floor(0.25*length(N))],
                          high50.N=sort(N)[floor(0.75*length(N))],
                          med.mean=median(mean,na.rm=T),
                          low50.mean=sort(mean)[floor(0.25*length(N))],
                          high50.mean=sort(mean)[floor(0.75*length(N))],
                          med.sd=median(sd,na.rm=T),
                          low50.sd=sort(sd)[floor(0.25*length(N))],
                          high50.sd=sort(sd)[floor(0.75*length(N))]
                          )
                        )
               ) %>% t %>% as.data.frame -> summary.trajs
        ## compute an early and a late TEA
        s <- 0.3/1.8^2
        N.early <- summary.trajs$med.N[21]
        TEA.early <- val*N.early+sqrt((val*N.early)^2 + val*N.early/s)
        N.late <- mean(summary.trajs$med.N[351:391])
        TEA.late <- val*N.late+sqrt((val*N.late)^2 + val*N.late/s)
        
        plot.new()
        plot.window(xlim=c(0,tmax), ylim=c(min(summary.trajs$low50.N), max(summary.trajs$high50.N)))
        axis(1); axis(2); box('plot')
        with(summary.trajs, polygon(c(time, rev(time)), c(high50.N, rev(low50.N)), col='skyblue', border=NA))
        with(summary.trajs, lines(time, med.N, lwd=2))
        lines(qg[,'time'], qg[,'N'], lwd=2, col='orange')
        if(val==dd[1]) mtext(side=2, 'Abundance', line=2)
        
        plot.new()
        plot.window(xlim=c(0,tmax), ylim=c(min(summary.trajs$low50.mean), max(summary.trajs$high50.mean)))
        axis(1); axis(2); box('plot')
        with(summary.trajs, polygon(c(time, rev(time)), c(high50.mean, rev(low50.mean)), col='skyblue', border=NA))
        with(summary.trajs, lines(time, med.mean, lwd=2))
        lines(qg[,'time'], qg[,'b'], lwd=2, col='orange')
        if(val==dd[1]) mtext(side=2, 'Trait mean', line=2)
        abline(h=TEA.late, lty=3, lwd=2)
        
        plot.new()
        plot.window(xlim=c(0,tmax), ylim=c(min(summary.trajs$low50.sd^2), max(summary.trajs$high50.sd^2)))
        axis(1); axis(2); box('plot')
        with(summary.trajs, polygon(c(time, rev(time)), c(high50.sd^2, rev(low50.sd^2)), col='skyblue', border=NA))
        with(summary.trajs, lines(time, med.sd^2, lwd=2))
        if(val==dd[1]) mtext(side=2, 'Trait variance', line=2)
        abline(h=0.75*(0.3*1.8)^2, lwd=2, col='orange')
        if(val==dd[2]) mtext(side=1, 'Time', line=2)

        with(nOff.early, plot(trait, nOff, col=gray(0.8), xlim=range(c(nOff.early$trait, nOff.late$trait)), ylim=range(c(nOff.early$nOff, nOff.late$nOff))), xlab='', ylab='')
        with(nOff.late, points(trait, nOff, col=gray(0.4)))
        ##abline(v=TEA.early, lty=2, lwd=2, col='purple')
        abline(v=TEA.late, lwd=2, col='blue')
        abline(v=5.4, lwd=2, col='orange')
        if(val==dd[1]) mtext(side=2, 'Lifetime reproductive\nsuccess', line=2)
        if(val==dd[2]) mtext(side=1, 'Birth rate trait', line=2)
            
    }

    dev.off()
    
}    

## Plot the bs=ds=0.04 case
png(filename="Simulation_trajectories_for_bs=ds=0.04.png", width=3, height=6, units='in', res=300)
par(mfcol=c(3,1), mar=c(3,3,1,0.5), oma=c(1,3,0,0), family="HersheySans")
val <- 0.04
## trajs just contains the trajectories themselves
trajs <- readRDS(file=paste0("logistic_GEM_storage_trajectories_bs=ds=",val,".RDS"))

## simulate the QG trajectory
## note: additive genetic variance is calculated based on the definition of heritability:
## h2 = Va/Vp, where Va is the additive genetic variance and Vp is the total phenotypic variance. The initial phenotypic variance is the initial traitsd^2, so Vp=(traitmean*traitcv)^2. Then Va=h2*(traitmean*traitcv)^2.
tmax <- 400
time.seq <- 0:(tmax-1)
sapply(time.seq,
       function(t)
           with(subset(trajs, time==t & N > 0),
                c(time=t,
                  med.N=median(N),
                  low50.N=sort(N)[floor(0.25*length(N))],
                  high50.N=sort(N)[floor(0.75*length(N))],
                  med.mean=median(mean),
                  low50.mean=sort(mean)[floor(0.25*length(N))],
                  high50.mean=sort(mean)[floor(0.75*length(N))],
                  med.sd=median(sd,na.rm=T),
                  low50.sd=sort(sd)[floor(0.25*length(N))],
                  high50.sd=sort(sd)[floor(0.75*length(N))]
                  )
                )
       ) %>% t %>% as.data.frame -> summary.trajs
qg <- ode(y=c(N=5, b=1.8), times=0:tmax, func=qg_model, parms=c(bs=val, ds=val, slope=0.3/1.8^2, V=0.75*(summary.trajs$med.sd[1])^2))



## compute an early and a late TEA
s <- 0.3/1.8^2
N.late <- mean(summary.trajs$med.N[351:391])
TEA.late <- val*N.late+sqrt((val*N.late)^2 + val*N.late/s)
        
plot.new()
plot.window(xlim=c(0,tmax), ylim=c(min(summary.trajs$low50.N), max(summary.trajs$high50.N)))
axis(1); axis(2); box('plot')
with(summary.trajs, polygon(c(time, rev(time)), c(high50.N, rev(low50.N)), col='skyblue', border=NA))
with(summary.trajs, lines(time, med.N, lwd=2))
lines(qg[,'time'], qg[,'N'], lwd=2, col='orange')
mtext(side=2, 'Abundance', line=3)

plot.new()
plot.window(xlim=c(0,tmax), ylim=c(min(summary.trajs$low50.mean), max(summary.trajs$high50.mean)))
axis(1); axis(2); box('plot')
with(summary.trajs, polygon(c(time, rev(time)), c(high50.mean, rev(low50.mean)), col='skyblue', border=NA))
with(summary.trajs, lines(time, med.mean, lwd=2))
lines(qg[,'time'], qg[,'b'], lwd=2, col='orange')
mtext(side=2, 'Trait mean', line=3)
abline(h=TEA.late, lty=3, lwd=2)

plot.new()
plot.window(xlim=c(0,tmax), ylim=c(min(summary.trajs$low50.sd^2), max(summary.trajs$high50.sd^2)))
axis(1); axis(2); box('plot')
with(summary.trajs, polygon(c(time, rev(time)), c(high50.sd^2, rev(low50.sd^2)), col='skyblue', border=NA))
with(summary.trajs, lines(time, med.sd^2, lwd=2))
mtext(side=2, 'Trait variance', line=2)
abline(h=0.75*(0.3*1.8)^2, lwd=2, col='orange')
mtext(side=1, 'Time', line=3)

dev.off()

    
